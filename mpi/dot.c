#include <stdio.h>
#include "algorithm.h"
#include "mark_time.h"
#include "utils.h"
#include "fileManager.h"
#include "dot_matrix.h"
#include "output.h"
#include "filtering.h"
#include "schedules.h"


int main (int argc, char* argv[]) {    
  struct config *cfg;
  MATCH_ARRAY_TYPE **wm;
  struct filemanager *fm;  /*  fasta/fastq file maneger  */
  struct outfile *output=NULL;
  struct sequence_t *seq; /*   contains information on the last loaded seq  */
  struct dot_matrix *dm;
  Dot_input *param;
  int i, res, rank, commsize;
  struct sequences_info * sinfo;  
  struct assigments * assings=NULL;  
  unsigned int myassingments, iter = 0, *num_assigs=NULL;
  unsigned long int *sizes=NULL, mysize;  
  long int myoffset, *offsets=NULL;  
  char tempOuputFile[MAX_PARAM_LEN+50], *path, *filename, *aux;
  #ifdef DEBUG_TIME  
    double tread=0, tcomp=0, twrite=0, twait=0, tcomun=0, start, end;
  #endif

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  
  wm = loadConfig(argc, argv, &cfg);
  
  verbose_output = cfg->verbose;  /*  rely on external variable  */

  /*  Initialize the file manager  */
  fm = filemanager_init (cfg->svalue);
  if (fm == NULL) {
    free_weights_matrix (wm);
    free (cfg);  /* free config params struct  */
    exit (EXIT_FAILURE);
  }
  fm->buffer_size=0;

    
  /*  Output file not specifyed  */
  if ((cfg->flags & O_FLAG) == 0){
    fprintf(stderr, "ERROR: You must specify an output file\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  else {
    aux = strdup(cfg->output_filename);
    path = dirname(aux);
    aux = strdup(cfg->output_filename);
    filename = basename(aux);
    sprintf(tempOuputFile,"%s/tmp-%s-%d.dot",path,filename,rank);

    output = output_create (tempOuputFile);
  }

  
  if(rank == 0){
    #ifdef DEBUG_TIME    
      start = MPI_Wtime();
    #endif    

    sinfo = filemanager_seq_count(fm);

    #ifdef DEBUG_TIME
      end = MPI_Wtime();
      tread += end - start;
      start = MPI_Wtime();
    #endif    
    
    if(cfg->schedule == BLOCK){
      assings = staticSchedule(sinfo, commsize);
      #ifdef DEBUG
        fprintf(stderr, "Using BLOCK Schedule\n");
      #endif
    }else if(cfg->schedule == BALANCED){
      assings = staticBalancedScheduleEnhanced(sinfo, commsize);
      #ifdef DEBUG
        fprintf(stderr, "Using BALANCED Schedule\n");
      #endif
    }else{
      fprintf(stderr, "Incorrect Schedule Mode\n");      
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    num_assigs = assings->num_assigs;
    offsets = assings->offsets;
    sizes = assings->size;
    print_header(output);

    #ifdef DEBUG_TIME
      end = MPI_Wtime();
      tcomp += end - start;
    #endif
    
  }  

  #ifdef DEBUG_TIME
    start = MPI_Wtime();
  #endif  
  MPI_Scatter(num_assigs, 1, MPI_UNSIGNED, &myassingments, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Scatter(offsets, 1, MPI_LONG, &myoffset, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Scatter(sizes, 1, MPI_UNSIGNED_LONG, &mysize, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  #ifdef DEBUG_TIME
    end = MPI_Wtime();
    tcomun += end - start;
  #endif
  
  #ifdef DEBUG
    fprintf(stderr,"Proc: %d | num_assigs: %u | offset: %ld | size: %lu\n", rank, myassingments, myoffset, mysize);    
  #endif  

  if(fseek (fm->pf, myoffset, SEEK_SET)){
    fprintf(stderr,"ERROR: fseek(%ld)\n",myoffset);
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  fm->finish = false;  
  
  /* algorithm run attempt */	
  seq = NULL;
  if ((param = dot_obj_init(cfg, wm, fm, output, 0)) == NULL) { 
      printf("Error in main() for param pointer\n"); 
      free_weights_matrix (wm);
      free (cfg);  /* free config params struct  */
      filemanager_destroy (fm);  /*  Releasefile managemant resources  */
      exit (EXIT_FAILURE); 
  }

  while(iter < myassingments) {
    #ifdef DEBUG_TIME
      start = MPI_Wtime();
    #endif    
    
    seq = filemanager_next_seq (fm, seq);
  
    #ifdef DEBUG_TIME
      end = MPI_Wtime();
      tread += end - start;

      start = MPI_Wtime();
    #endif
    if (verbose_output) {
        if (seq != NULL) {
            fprintf (stderr,"Rank %d Processing %s length: %ld bp (%u of %u)\n",rank, seq->label, seq->sequence_size, iter+1, myassingments);
            fflush (stderr);
        }
    }
    if ( seq == NULL ) {
      #ifdef DEBUG_TIME
        end = MPI_Wtime();
        tcomp += end - start;    
      #endif
        break;
    } 
    else {      
      param->sequence = seq->sequence;
      param->IDSeq = seq->label;

      dm = dot_init (seq, wm);
      if (dm == NULL) {
          perror("Error in creating dot_matrix\n");
          return 1; 
      }
      param->matrix = dm;      
      res = start_TRs_search(param);
      if (res != 0) {
          perror("Something was wrong searching TRs\n");
          return 1; 
      }
      #ifdef DEBUG_TIME
        end = MPI_Wtime();
        tcomp += end - start;
      #endif             

      if (param->TRs_bundle->trs_found_offset > 0) {
        #ifdef DEBUG_TIME
          start = MPI_Wtime();
        #endif        
        /* FINAL LIST FILTERING */
        filter (param->TRs_bundle, param->config_params);
        #ifdef DEBUG_TIME
          end = MPI_Wtime();
          tcomp += end - start;

          start = MPI_Wtime();
        #endif        
        print_TRs_list_toFile (param->output , param->IDSeq, param->sequence, param->TRs_bundle);
        #ifdef DEBUG_TIME
          end = MPI_Wtime();
          twrite += end - start;
        #endif
      } else {
        printf("No Tandem Repeats found in %s\n", param->IDSeq);        
      }      
      reset_dot_obj(param);
      dot_free(dm);      
      iter++;
    }    
  }

  if (iter != myassingments) {
        fprintf(stderr, "Rank %d Error processing sequences. Done %u of %u\n", rank, iter, myassingments);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 1);
  }

  output_destroy (output);  

  #ifdef DEBUG_TIME
    start = MPI_Wtime();    
  #endif  

  MPI_Barrier(MPI_COMM_WORLD);

  #ifdef DEBUG_TIME
    end = MPI_Wtime();
    twait += end - start;
  #endif  
  
  if(rank == 0){
    #ifdef DEBUG_TIME
      start = MPI_Wtime();
    #endif    

    if(merge_files(cfg->output_filename, commsize)){
      perror("ERROR merging temp output files\n");      
    }
    if(remove_tmpfiles(cfg->output_filename, commsize)){
      perror("ERROR deleting tempfiles\n");      
    }
    
    #ifdef DEBUG_TIME
      end = MPI_Wtime();
      twrite += end - start;
    #endif
  }
  
  destroy_dot_obj(&param);      
  filemanager_destroy  (fm);
  free_weights_matrix (wm);
  free (cfg);  /* free config params struct  */

  #ifdef DEBUG_TIME
    int numDebugData=8, idx;
    double debug_data[] = {rank, myassingments, mysize, tread, tcomp, tcomun, twrite, twait};
    double *ddg=NULL;
    if(rank==0){
        ddg = (double *) malloc(sizeof(double)*numDebugData*commsize);
    }

    MPI_Gather(debug_data, numDebugData, MPI_DOUBLE, ddg, numDebugData, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(rank==0){
      fprintf(stderr, "| RANK | #seqs | bps | t_reading | t_computing | t_communicating | t_writing | t_waiting | TOTAL |\n");   
      for(i=0; i<commsize; i++){
        idx = numDebugData * i;
        double total = ddg[idx+3] + ddg[idx+4] + ddg[idx+5] + ddg[idx+6] + ddg[idx+7];
        fprintf(stderr,"%d\t%d\t%ld\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\t%.2lf\n",(int)ddg[idx], (int)ddg[idx+1], (long int)ddg[idx+2], ddg[idx+3], ddg[idx+4], ddg[idx+5], ddg[idx+6], ddg[idx+7], total);
      }
      free(ddg);
    }  
  #endif

  MPI_Finalize();

  return EXIT_SUCCESS;
}
