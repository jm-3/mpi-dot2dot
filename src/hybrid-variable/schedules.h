#ifndef _SCHEDULES_H_
#define _SCHEDULES_H_

#include "utils.h"
#include "fileManager.h"


struct assigments
{
	unsigned int * num_assigs;
	long int * offsets;
	long int * size;
};

struct assigments * staticSchedule(struct sequences_info *, int);
struct assigments * staticBalancedSchedule(struct sequences_info *, int);
struct assigments * staticBalancedScheduleEnhanced(struct sequences_info *, int);
struct assigments * assings_init(int n);
void assings_free(struct assigments *a);
int *shareCoresPerProcBySizes(int long *sizes, int commsize, int total_cores, int mincoresperproc, int maxcoresperproc);

#endif