#include <stdio.h>
#include "algorithm.h"
#include "mark_time.h"
#include "utils.h"
#include "fileManager.h"
#include "dot_matrix.h"
#include "output.h"
#include "filtering.h"

_Bool verbose_output;


int main (int argc, char* argv[]) {  
  struct paramList *par_list=NULL;
  struct config *cfg;
  MATCH_ARRAY_TYPE **wm;
  struct filemanager *fm;  /*  fasta/fastq file maneger  */
  struct outfile *output;
  struct sequence_t *seq; /*   contains information on the last loaded seq  */
  struct dot_matrix *dm;
  Dot_Thread_input *param;
  int res;
  
  cfg = command_line_parser (argc, argv);

  verbose_output = cfg->verbose;  /*  rely on external variable  */

  if ((cfg->flags & U_FLAG) != 0) {
    print_usage();
    exit(EXIT_SUCCESS);
  }
  if ((cfg->flags & (S_FLAG + C_FLAG)) != S_FLAG + C_FLAG) {
    printf("Missing REQUIRED files: SequenceFile or ConfigurationFile\n");
    print_usage();
    exit (EXIT_FAILURE);
  }
  
  if ((par_list = loadConfigFromFile (cfg->cvalue)) == NULL) {
    perror ("opening configuration file\n");
    exit (EXIT_FAILURE);
  }

  param_list_parser (par_list, cfg); /*  Update cfg  */
  wm = read_weights_matrix (par_list);
  free_paramList (par_list);   /*  From here it is no longer used  */
  if (wm == NULL) exit (EXIT_FAILURE);
  /*
    printf("NVALUE: %s\n", cfg->output_filename);
    printf("XVALUE: %d\n", cfg->xvalue);
    printf("FVALUE: %f\n", cfg->fvalue);
    printf("JVALUE: %d\n", cfg->jvalue);
    printf("GVALUE: %d\n", cfg->gvalue);
  */

  
  /*  Integrity check - must be moved in the appropr. place */
  if ((cfg->xvalue < cfg->nvalue) && (cfg->xvalue > 0)) {
    printf("max_length must be higher than min_length\n");
    print_usage();
    free_weights_matrix (wm);
    free (cfg);  /* free config params struct  */
    exit(EXIT_FAILURE);
  }

  if ((cfg->fvalue <= 0) || (cfg->fvalue > 1)) {
    printf("MinMatch parameter must be included in (0,1]\n");
    print_usage();
    free_weights_matrix (wm);
    free (cfg);  /* free config params struct  */
    exit(EXIT_FAILURE);
  }	
  /*  Initialize the file manager  */
  fm = filemanager_init (cfg->svalue);
  if (fm == NULL) {
    free_weights_matrix (wm);
    free (cfg);  /* free config params struct  */
    exit (EXIT_FAILURE);
  }

    
  /*  Output file not specifyed  */
  if ((cfg->flags & O_FLAG) == 0)
    output = output_create (NULL);
  else 
    output = output_create (cfg->output_filename);

  print_header (output);
  
   
  /* algorithm run attempt */	
  seq = NULL;
  if ((param = dot_Thread_obj_init(cfg, wm, fm, output, 0)) == NULL) { 
      printf("Error in main() for param pointer\n"); 
      free_weights_matrix (wm);
      free (cfg);  /* free config params struct  */
      filemanager_destroy (fm);  /*  Releasefile managemant resources  */
      exit (EXIT_FAILURE); 
  }
  do {       
    seq = filemanager_next_seq (fm, seq);
    if (verbose_output) {
        if (seq != NULL) {
            fprintf (stderr,"Processing %s length: %ld bp\n",seq->label, seq->sequence_size);
            fflush (stderr);
        }
    }
    if ( seq == NULL ) {
        continue;
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

  	/* Print TRs on file  */       
#ifdef DEBUG_THREAD      
      printf("WRITING...\n");        
#endif          
      if (param->thread_TRs_bundle->trs_found_offset > 0) {
        /* FINAL LIST FILTERING */
        filter (param->thread_TRs_bundle, param->config_params);
        print_TRs_list_toFile (param->output , param->IDSeq, param->sequence, param->thread_TRs_bundle);
      } else {
        printf("No Tandem Repeats found in %s\n", param->IDSeq);        
      }      
      reset_dot_Thread_obj(param);
      dot_free(dm);
    }    

  } while (seq != NULL);

  if ((param->file_manager)->finish != true) {
        perror ("Error reading input file\n");
        return 1;
  }
  
  destroy_dot_Thread_obj(&param);    
  output_destroy (output);
  filemanager_destroy  (fm);
  free_weights_matrix (wm);
  free (cfg);  /* free config params struct  */
  exit (EXIT_SUCCESS);
}
