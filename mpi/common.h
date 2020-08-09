#ifndef __COMMON__
#define __COMMON__

#include <stdbool.h>  /*  Makes _Bool and true/false available  */
#include <mpi.h>
#include "utils.h"

#define MATCH_ARRAY_TYPE char
#define MAXIND_ARRAY (256 * sizeof (char))
#define MAX_PARAM_LEN 1024
#define MAX_LABEL_LENGTH  65
/*#define BUFF_SIZE 1048576 */
#define BUFF_SIZE 16777216

_Bool verbose_output;


/*  Native type for filtering  */
typedef enum {NONE, THRESHOLD, LIGHT, FAIR, HEAVY} strategy_t;
/*  Native typer for file manager  */
typedef enum {UNKNOWN, FASTA, FASTQ} seqfile_t;
/* Types of scheduling */
typedef enum {WRONG, BLOCK, BALANCED} schedule_t;

#endif /* __COMMON__ */
