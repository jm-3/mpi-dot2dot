#ifndef _SCHEDULES_H_
#define _SCHEDULES_H_

#include "utils.h"
#include "fileManager.h"


struct assigments
{
	unsigned int * num_assigs;
	long int * offsets;
	unsigned long int * size;
};

struct assigments * staticSchedule(struct sequences_info *, int);
struct assigments * staticBalancedSchedule(struct sequences_info *, int);
struct assigments * staticBalancedScheduleEnhanced(struct sequences_info *, int);

#endif