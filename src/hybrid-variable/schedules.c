#include "schedules.h"

struct assigments * staticSchedule(struct sequences_info *sinfo, int commsize){
	struct assigments * assings;	
	unsigned int i, quotient, remainder, counter;
	
	assings = assings_init(commsize);

	quotient = sinfo->num_seqs / commsize;
	if(quotient == 0)
		fprintf(stderr, "WARNING: Number of secuences (%u) is lower than number of processes (%d). Some of them will be idle\n", sinfo->num_seqs, commsize);

	for(i=0; i<commsize; i++){
		assings->num_assigs[i] = quotient;
	}

    remainder = sinfo->num_seqs % commsize;
	for(i=0; i<remainder; i++){
		assings->num_assigs[i]++;	
	}

	assings->offsets[0] = 0;
	counter = assings->num_assigs[0];
	for(i=1; i<commsize; i++){
		if(assings->num_assigs[i] > 0)
			assings->offsets[i] = sinfo->offsets[counter];	
		else
			assings->offsets[i] = -1;
		counter += assings->num_assigs[i];
	}

	#ifdef DSCHEDULER
		fprintf(stderr, "Total sequences: %u\n", sinfo->num_seqs);
		for(i=0; i<commsize; i++)
			fprintf(stderr, "Rank %d NumAssigs %u\n", i, assings->num_assigs[i]);
	#endif

	return assings;

}

struct assigments * staticBalancedSchedule(struct sequences_info *sinfo, int commsize){
	struct assigments * assings;
	unsigned int i, counter=0;
	long int off=0;
	unsigned long int quotient, sumsizes=0, accum;

	assings = assings_init(commsize);

	for(i=0; i<sinfo->num_seqs;i++){
		sumsizes += sinfo->sizes[i];		
	}

	quotient = sumsizes / commsize;

	for(i=0; i<commsize; i++){
		accum = 0;
		assings->num_assigs[i] = 0;	
		while(accum < quotient && counter < sinfo->num_seqs){
			accum += sinfo->sizes[counter];
			assings->num_assigs[i]++;
			counter++;
		};
		if(counter <= sinfo->num_seqs){
			assings->offsets[i] = sinfo->offsets[off];
			assings->size[i] = accum;			
		}
		else{
			assings->offsets[i] = -1;
			assings->size[i] = 0;
		}
		off += assings->num_assigs[i];
	}

	if(verbose_output){
		double stddev;
		unsigned long int min, max;
		stddev = calcStdDevli(assings->size, commsize);
		min = getMinValueli(assings->size, commsize);
		max = getMaxValueli(assings->size, commsize);
		fprintf(stderr, "Static Balanced Scheduler: sum %lu, np %d, average %lu, stddev %.3lf, min size: %lu, max size %lu\n", sumsizes, commsize, quotient, stddev, min, max);		
	}

	return assings;
}

struct assigments * staticBalancedScheduleEnhanced(struct sequences_info *sinfo, int commsize){
	struct assigments * assings;
	unsigned int i, counter=0, maxnumassigsperproc;
	long int off=0;
	long int average, sumsizes=0, accum, totalaccum=0;	

	assings = assings_init(commsize);

	if (sinfo->num_seqs < commsize){
		fprintf(stderr, "WARNING: Number of secuences (%u) is lower than number of processes (%d). Some of them will be idle\n", sinfo->num_seqs, commsize);
	}


	for(i=0; i<sinfo->num_seqs; i++){
		sumsizes += sinfo->sizes[i];
	}
	
	for(i=0; i<commsize; i++){
		accum = 0;
		assings->num_assigs[i] = 0;		
		maxnumassigsperproc = (sinfo->num_seqs - counter) / (commsize - i);
		average = (sumsizes - totalaccum) / (commsize - i);
		while(counter < sinfo->num_seqs){			
			if((((accum + sinfo->sizes[counter]) > (average*1.1)) && (accum > average*0.8)) || assings->num_assigs[i] == maxnumassigsperproc){				
				break;
			}
			accum += sinfo->sizes[counter];
			assings->num_assigs[i]++;
			counter++;
		};
		if(counter <= sinfo->num_seqs){
			assings->offsets[i] = sinfo->offsets[off];
			assings->size[i] = accum;			
		}
		else{
			assings->offsets[i] = -1;
			assings->size[i] = 0;
		}
		off += assings->num_assigs[i];
		totalaccum += accum;
	}

	#ifdef DSCHEDULER
		double stddev;
		unsigned long int min, max;
		average = sumsizes / commsize;
		stddev = calcStdDevli(assings->size, commsize);
		min = getMinValueli(assings->size, commsize);
		max = getMaxValueli(assings->size, commsize);
		fprintf(stderr, "Static Balanced Scheduler Enhanced: sum %lu, np %d, average %lu, stddev %.3e, min size: %lu, max size %lu\n", sumsizes, commsize, average, stddev, min, max);
		fprintf(stderr,"Sizes per proc: ");
		for(i=0; i<commsize; i++)		
			fprintf(stderr, "Rank %d NumAssigs %u Size %ld\n", i, assings->num_assigs[i], assings->size[i]);
		fprintf(stderr,"\n");
	#endif

	return assings;
}


struct assigments * assings_init(int n){
	struct assigments *a;

	a = (struct assigments *) memalloc(sizeof(struct assigments), "Error allocating memory for struct assigments\n");
	a->num_assigs = (unsigned int *) memalloc(sizeof(unsigned int) * n, "Error allocating memory for assings->num_assigs\n");
	a->offsets = (long int *) memalloc(sizeof(long int) * n, "Error allocating memory for assings->offsets\n");
	a->size = (long int *) memalloc(sizeof(long int) * n, "Error allocating memory for assings->size\n");

	return a;
}

void assings_free(struct assigments *a){
	if(a->num_assigs != NULL) free(a->num_assigs);
	if(a->offsets != NULL) free(a->offsets);
	if(a->size != NULL) free(a->size);
	
}


/*******************************************************************************************/
/* hybrid */
int *shareCoresPerProcBySizes(int long *sizes, int commsize, int total_cores){
	int *cores, i, j, remaining_cores, mproc;
	int mincoresperproc = 1;

	if(total_cores < commsize) fprintf(stderr, "Warning: there are less cores per node than processes\n");

	cores = (int *) memalloc(sizeof(int)*commsize, "Error allocating memory for cores\n");
	// initial share
	for(i=0; i<commsize; i++)
		cores[i] = mincoresperproc;

	remaining_cores = total_cores - (mincoresperproc * commsize);

	for(i=0; i<remaining_cores; i++) { 
        mproc = 0; 
    	for(j=1; j<commsize; j++) { 
            if( sizes[j]/cores[j] > sizes[mproc]/ cores[mproc] ) { 
    			mproc = j; 
    		}
    	}
    	cores[mproc]++;
    }
        
	return cores;
}
