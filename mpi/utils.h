#ifndef __UTILS__
#define __UTILS__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <libgen.h>
#include <string.h>

void print_usage (void);
void print_version (void);
void *memalloc(size_t,char*);
void *memrealloc(void*,size_t,char*);
double calcStdDev(unsigned long int *, int );
unsigned long int getMinValue(unsigned long int *, int);
unsigned long int getMaxValue(unsigned long int *, int);
long int getMaxValueli(long int *, int);
char* toLower(char* s);
const char *get_filename_ext(const char *filename);
char *remove_filename_ext(char* myStr);

#endif /* UTILS */
