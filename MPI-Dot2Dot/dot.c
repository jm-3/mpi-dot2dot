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
  Dot_Thread_input *param;
  int res, rank, commsize;
  struct sequences_info * sinfo;  
  struct assigments * assings=NULL;  
  unsigned int myassingments, iter = 0, *num_assigs=NULL;
  unsigned long int *sizes=NULL, mysize;  
  long int myoffset, *offsets=NULL;  
  char tempOuputFile[MAX_PARAM_LEN+50];

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
    perror("ERROR: You must specify an output file\n");
    MPI_Abort(MPI_COMM_WORLD, 1);
    //output = output_create (NULL);
  }
  else {
    sprintf(tempOuputFile,"tmp-%s-%d.dot",cfg->output_filename,rank);
#ifdef DEBUG
    printf("Proc: %d | outputfile: %s\n", rank, tempOuputFile);
#endif
    output = output_create (tempOuputFile);
  }

  
  if(rank == 0){
    double start, end;    

    start = MPI_Wtime();

    sinfo = filemanager_seq_count(fm);
    assings = staticBalancedSchedule(sinfo, commsize);
    num_assigs = assings->num_assigs;
    offsets = assings->offsets;
    sizes = assings->size;
    print_header(output);

    end = MPI_Wtime();

    fprintf(stderr, "Rank 0. Time reading secuences %.3lf seconds\n", (end - start));
  }  

  
  MPI_Scatter(num_assigs, 1, MPI_UNSIGNED, &myassingments, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  MPI_Scatter(offsets, 1, MPI_LONG, &myoffset, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  MPI_Scatter(sizes, 1, MPI_UNSIGNED, &mysize, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
  
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
  if ((param = dot_Thread_obj_init(cfg, wm, fm, output, 0)) == NULL) { 
      printf("Error in main() for param pointer\n"); 
      free_weights_matrix (wm);
      free (cfg);  /* free config params struct  */
      filemanager_destroy (fm);  /*  Releasefile managemant resources  */
      exit (EXIT_FAILURE); 
  }
  while(iter < myassingments) {    
    seq = filemanager_next_seq (fm, seq);
    if (verbose_output) {
        if (seq != NULL) {
            fprintf (stderr,"Rank %d Processing %s length: %ld bp (%u of %u)\n",rank, seq->label, seq->sequence_size, iter+1, myassingments);
            fflush (stderr);
        }
    }
    if ( seq == NULL ) {
        break;
    } 
    else {      
      if ((param->sequence = copy_string(seq->sequence)) == NULL) {
          perror("Error in copying sequence in dot_thread_fn\n");            
          return 1;
      };
      if ((param->IDSeq = copy_string(seq->label)) == NULL) {
          perror("Error in copying label in dot_thread_fn\n");              
          return 1;
      };
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

      if (param->thread_TRs_bundle->trs_found_offset > 0) {
        /* FINAL LIST FILTERING */
        filter (param->thread_TRs_bundle, param->config_params);
        print_TRs_list_toFile (param->output , param->IDSeq, param->sequence, param->thread_TRs_bundle);
      } else {
        printf("No Tandem Repeats found in %s\n", param->IDSeq);        
      }      
      reset_dot_Thread_obj(param);
      dot_free(dm);      
      iter++;
    }    
  }

  if (iter != myassingments) {
        fprintf(stderr, "Rank %d Error processing secuences. Done %u of %u\n", rank, iter, myassingments);
        fflush(stderr);
        MPI_Abort(MPI_COMM_WORLD, 1);
  }

#ifdef DEBUG
  fprintf(stderr,"Proc: %d | Waiting in Barrier\n", rank);    
#endif
  output_destroy (output);  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank == 0){
    double start, end;    

    start = MPI_Wtime();

    if(merge_files(cfg->output_filename, commsize)){
      perror("ERROR merging temp output files\n");      
    }
    if(remove_tmpfiles(cfg->output_filename, commsize)){
      perror("ERROR deleting tempfiles\n");      
    }
    sync();
    
    end = MPI_Wtime();

    fprintf(stderr, "Rank 0. Time reading secuences %.3lf seconds\n", (end - start));
  }
  
  destroy_dot_Thread_obj(&param);      
  filemanager_destroy  (fm);
  free_weights_matrix (wm);
  free (cfg);  /* free config params struct  */

  MPI_Finalize();

  return EXIT_SUCCESS;
}
