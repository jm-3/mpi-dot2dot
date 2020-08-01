#include "schedules.h"

struct assigments * staticSchedule(struct sequences_info *sinfo, int commsize){
	struct assigments * assings;	
	int i, j, quotient, remainder, counter;
	unsigned long int size, index;
	
	assings = (struct assigments *) memalloc(sizeof(struct assigments), "Error allocating memory for struct assigments\n");
	assings->num_assigs = (unsigned int *) memalloc(sizeof(unsigned int) * commsize, "Error allocating memory for assings->num_assigs\n");
	assings->offsets = (long int *) memalloc(sizeof(long int) * commsize, "Error allocating memory for assings->offsets\n");
	assings->size = (unsigned long int *) memalloc(sizeof(unsigned long int) * commsize, "Error allocating memory for assings->size\n");

	quotient = sinfo->num_seqs / commsize;
	if(quotient == 0)
		fprintf(stderr, "WARNING: Number of secuences (%ld) is lower than number of processes (%d). Some of them will be idle\n", sinfo->num_seqs, commsize);

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

	index = 0;
	for(i=0; i<commsize; i++){
		size = 0;
		for(j=0; j<assings->num_assigs[i]; j++){
			size += sinfo->sizes[index];
			index++;
		}
		assings->size[i] = size;

	}

	return assings;

}

struct assigments * staticBalancedSchedule(struct sequences_info *sinfo, int commsize){
	struct assigments * assings;
	unsigned int i, counter=0;
	long int off=0;
	unsigned long int quotient, sumsizes=0, accum;

	assings = (struct assigments *) memalloc(sizeof(struct assigments), "Error allocating memory for struct assigments\n");
	assings->num_assigs = (unsigned int *) memalloc(sizeof(unsigned int) * commsize, "Error allocating memory for assings->num_assigs\n");
	assings->offsets = (long int *) memalloc(sizeof(long int) * commsize, "Error allocating memory for assings->offsets\n");
	assings->size = (unsigned long int *) memalloc(sizeof(unsigned long int) * commsize, "Error allocating memory for assings->size\n");

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
		stddev = calcStdDev(assings->size, commsize);
		min = getMinValue(assings->size, commsize);
		max = getMaxValue(assings->size, commsize);
		fprintf(stderr, "Static Balanced Scheduler: sum %lu, np %d, average %lu, stddev %.3lf, min size: %lu, max size %lu\n", sumsizes, commsize, quotient, stddev, min, max);		
	}

	return assings;
}

struct assigments * staticBalancedScheduleEnhanced(struct sequences_info *sinfo, int commsize){
	struct assigments * assings;
	unsigned int i, counter=0, maxnumassigsperproc;
	long int off=0, *ranges, maxrange;
	unsigned long int average, sumsizes=0, accum, totalaccum=0;	

	assings = (struct assigments *) memalloc(sizeof(struct assigments), "Error allocating memory for struct assigments\n");
	assings->num_assigs = (unsigned int *) memalloc(sizeof(unsigned int) * commsize, "Error allocating memory for assings->num_assigs\n");
	assings->offsets = (long int *) memalloc(sizeof(long int) * commsize, "Error allocating memory for assings->offsets\n");
	assings->size = (unsigned long int *) memalloc(sizeof(unsigned long int) * commsize, "Error allocating memory for assings->size\n");

	for(i=0; i<sinfo->num_seqs;i++){
		sumsizes += sinfo->sizes[i];	
	}
	
	average = sumsizes / commsize;
	maxnumassigsperproc = sinfo->num_seqs / commsize;

	ranges = (long int *) memalloc(sizeof(long int) * sinfo->num_seqs, "Error allocating memory for ranges\n");
	for(i=0; i<sinfo->num_seqs; i++)
		ranges[i] = sinfo->sizes[i] - average;
	maxrange = getMaxValueli(ranges, sinfo->num_seqs);

	for(i=0; i<commsize; i++){
		accum = 0;
		assings->num_assigs[i] = 0;		
		maxnumassigsperproc = (sinfo->num_seqs - counter) / (commsize - i);
		average = (sumsizes - totalaccum) / (commsize - i);
		while(counter < sinfo->num_seqs){			
			if((((accum + sinfo->sizes[counter]) > (average*1.1)) && (accum > average*0.5)) || assings->num_assigs[i] == maxnumassigsperproc){				
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

	if(verbose_output){
		double stddev;
		unsigned long int min, max;
		average = sumsizes / commsize;
		stddev = calcStdDev(assings->size, commsize);
		min = getMinValue(assings->size, commsize);
		max = getMaxValue(assings->size, commsize);
		fprintf(stderr, "Static Balanced Scheduler Enhanced: sum %lu, np %d, average %lu, stddev %.3e, min size: %lu, max size %lu\n", sumsizes, commsize, average, stddev, min, max);
		fprintf(stderr,"Sizes per proc: ");
		for(i=0; i<commsize; i++)		
			fprintf(stderr," %lu", assings->size[i]);
		fprintf(stderr,"\n");
	}

	return assings;
}
