#include <stdio.h>
#include "algorithm.h"
#include "mark_time.h"
#include "utils.h"
#include "fileManager.h"
#include "dot_matrix.h"
#include "output.h"
#include "filtering.h"
#include "schedules.h"
#include <omp.h>


int main (int argc, char* argv[]) {    
  struct config *cfg;
  MATCH_ARRAY_TYPE **wm;
  struct filemanager *fm;  /*  fasta/fastq file maneger  */
  struct outfile *output=NULL;
  struct sequence_t *seq; /*   contains information on the last loaded seq  */
  struct dot_matrix *dm;
  Dot_input *param=NULL;
  int res, rank, commsize;
  struct sequences_info * sinfo;  
  struct assigments * assings=NULL;  
  unsigned int myassingments, iter = 0, *num_assigs=NULL;
  long int myoffset, *offsets=NULL, mysize=0, *sizes=NULL;  
  #ifdef DEBUG_TIME  
    double tread=0, tcomp=0, twrite=0, twait=0, tcomun=0, start, end;
  #endif
  MPI_Request request;
  MPI_Status status;
  int provided;
  char nodename[MPI_MAX_PROCESSOR_NAME];
  int nodenamelen, color, key, intranode_commsize, intranode_rank, total_cores, *coresperproc=NULL, mycores;
  long int *intranode_sizes=NULL;
  MPI_Comm intranode_comm;


  MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
  if (provided < MPI_THREAD_FUNNELED) MPI_Abort(MPI_COMM_WORLD, 1);
  MPI_Comm_size(MPI_COMM_WORLD, &commsize);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // setting num threads
  /*if (rank == 0) {
    char *nthreads_str = getenv("OMP_NUM_THREADS");
    if (nthreads_str)
      nthreads = atoi(nthreads_str);
    else
      nthreads = 1;
  }
  MPI_Bcast(&nthreads, 1, MPI_INT, 0, MPI_COMM_WORLD);
  omp_set_num_threads(nthreads);

  #ifdef DEBUG
    if(rank == 0){      
      fprintf(stderr, "Config: %d proceses, %d threads per proc\n", commsize, omp_get_max_threads());
    }
  #endif*/

  
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
    output = output_create_tmp (cfg->output_filename, rank);
  }

  
  if(rank == 0){
    #ifdef DEBUG_TIME    
      start = MPI_Wtime();
    #endif    

    sinfo = filemanager_seq_count(fm, cfg->schedule);

    #ifdef DEBUG_TIME
      end = MPI_Wtime();
      tread += end - start;
      start = MPI_Wtime();
    #endif    
    
    if(cfg->schedule == BLOCK){
      #ifdef DEBUG
        fprintf(stderr, "Using BLOCK Schedule\n");
      #endif      
      assings = staticSchedule(sinfo, commsize);
    }else if(cfg->schedule == BALANCED){
      #ifdef DEBUG
        fprintf(stderr, "Using BALANCED Schedule\n");
      #endif      
      assings = staticBalancedScheduleEnhanced(sinfo, commsize);
    }else{
      fprintf(stderr, "Incorrect Schedule Mode\n");      
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    sequences_info_free(sinfo);

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
  MPI_Scatter(sizes, 1, MPI_LONG, &mysize, 1, MPI_LONG, 0, MPI_COMM_WORLD);
  
  #ifdef DEBUG_TIME
    end = MPI_Wtime();
    tcomun += end - start;
  #endif

  if(rank == 0) assings_free(assings);  

  /* communication along others procs in the same node */
  if(MPI_Get_processor_name(nodename, &nodenamelen) != MPI_SUCCESS) { fprintf(stderr, "Error getting processor name\n"); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); }
  color = hash(nodename);
  key = rank;
  if(MPI_Comm_split(MPI_COMM_WORLD, color, key, &intranode_comm) != MPI_SUCCESS) { fprintf(stderr, "Error creating intranode comm\n"); MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); }
   
  MPI_Comm_size(intranode_comm, &intranode_commsize);
  MPI_Comm_rank(intranode_comm, &intranode_rank);
  
  if(intranode_rank == 0)
    intranode_sizes = (long *) malloc(sizeof(long) * intranode_commsize);

  MPI_Gather(&mysize, 1, MPI_LONG, intranode_sizes, 1, MPI_LONG, 0, intranode_comm);

  /* distribute cores among procs */
  if(intranode_rank == 0){
    total_cores = omp_get_num_procs();

    coresperproc = shareCoresPerProcBySizes(intranode_sizes, intranode_commsize, total_cores);

    #ifdef DEBUG  
      fprintf(stderr, "Node: %s Total cores: %d ", nodename, total_cores);
      for(int i=0; i<intranode_commsize; i++)
        fprintf(stderr, "| Proc %d, size %ld, num_cores %d", i, intranode_sizes[i], coresperproc[i]);    
      fprintf(stderr, "\n");
    #endif
  }
  MPI_Scatter(coresperproc, 1, MPI_INT, &mycores, 1, MPI_INT, 0, intranode_comm);
  omp_set_num_threads(mycores);

  if(myassingments > 0){
    if(fseek (fm->pf, myoffset, SEEK_SET)){
      fprintf(stderr,"ERROR: fseek(%ld)\n",myoffset);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
    seq = NULL;
    if ((param = dot_obj_init(cfg, wm, fm, output, 0)) == NULL) { 
      fprintf(stderr, "Error in main() for param pointer\n"); 
      free_weights_matrix (wm);
      free (cfg);  /* free config params struct  */
      filemanager_destroy (fm);  /*  Releasefile managemant resources  */
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
    }
  }
  
  fm->finish = false;  
  
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
      mysize += seq->sequence_size;

      param->sequence = seq->sequence;
      param->IDSeq = seq->label;

      dm = dot_init (seq, wm);
      if (dm == NULL) {
          perror("Error in creating dot_matrix\n");
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); 
      }
      param->matrix = dm;      
      res = start_TRs_search(param);
      if (res != 0) {
          perror("Something was wrong searching TRs\n");
          MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
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

  #ifdef DEBUG
    fprintf(stderr,"Proc: %d | num_assigs: %u | offset: %ld | size: %ld\n", rank, myassingments, myoffset, mysize);    
  #endif   

  #ifdef DEBUG_TIME
    start = MPI_Wtime();    
  #endif  

  MPI_Ibarrier(MPI_COMM_WORLD, &request);

  if(rank == 0)
    MPI_Wait(&request, &status);
  else{
    int flag;
    do{
      sleep(1);
      MPI_Test(&request, &flag, MPI_STATUS_IGNORE);
    }while(!flag);
  }

  #ifdef DEBUG_TIME
    end = MPI_Wtime();
    twait += end - start;
  #endif  
  
  if(rank == 0){
    #ifdef DEBUG_TIME
      start = MPI_Wtime();
    #endif    

    if(merge_files(cfg->output_filename, commsize)){
      perror("ERROR merging temp output files");      
    }
    if(remove_tmpfiles(cfg->output_filename, commsize)){
      perror("ERROR deleting tempfiles");      
    }
    
    #ifdef DEBUG_TIME
      end = MPI_Wtime();
      twrite += end - start;
    #endif
  }
  
  if(param != NULL) destroy_dot_obj(&param);      
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
      for(int i=0; i<commsize; i++){
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
