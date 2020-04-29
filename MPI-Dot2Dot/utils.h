#ifndef __UTILS__
#define __UTILS__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void print_usage (void);
void print_version (void);
void *memalloc(size_t,char*);
void *memrealloc(void*,size_t,char*);
float calcStdDev(unsigned long int *, int );
unsigned long int getMinValue(unsigned long int *, int);
unsigned long int getMaxValue(unsigned long int *, int);
long int getMaxValueli(long int *, int);

#endif /* UTILS */
