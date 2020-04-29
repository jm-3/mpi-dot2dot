#include <stdio.h>
#include "fileManager.h"
#include "utils.h"
#include "configReader.h"
#include "schedules.h"

int main(int argc, char* argv[]){
	struct config *cfg;
	struct filemanager *fm;  /*  fasta/fastq file maneger  */	
	struct sequences_info * sinfo;
	struct assigments * assings; 

	cfg = command_line_parser (argc, argv);

	  /*  Initialize the file manager  */
  	fm = filemanager_init (cfg->svalue);
  	if (fm == NULL) {    	
    	free (cfg);  /* free config params struct  */
    	exit (EXIT_FAILURE);
  	}

  	sinfo = filemanager_seq_count(fm);

  	printf("Num secuencias: %u\n", sinfo->num_seqs);
  	printf("Offsets: ");
  	for(int i=0; i<sinfo->num_seqs; i++)
  		printf("%lu ",sinfo->offsets[i]);
  	printf("\n");

	assings = staticSchedule(sinfo, 24);

	for(int i=0; i<24; i++)
		printf("%d %u %ld\n",i,assings->num_assigs[i], assings->offsets[i]);



	return 0;
}
