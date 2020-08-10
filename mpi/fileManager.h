#ifndef FILEMANAGER_H_
#define FILEMANAGER_H_

#include "common.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/*  #define MAX_LABEL_LENGTH  65  defined in common.h  */
/*  #define MAX_PARAM_LEN 1024  defined in common.h  */
/*  typedef enum {UNKNOWN, FASTA, FASTQ} seqfile_t;  Declared in common.h */
typedef enum {H_PRE_SI, H_PRE_LABEL, H_LABEL, H_POST_LABEL, SEQUENCE, FQ_PLUS, FQ_SCORE} parse_status_t;


struct filemanager {
  seqfile_t filetype;
  char buffer[BUFF_SIZE];
  long int buffer_size;   /*  Buffer size for malloc  */
  unsigned int offset;
  char empty_identifier[MAX_LABEL_LENGTH];  /*  Used when a sequence has no identifier  */
  FILE *pf;
  _Bool finish;  /*  True if no more sequence is to be read  */  
};

struct sequence_t {
  char label[MAX_LABEL_LENGTH];/*header*/
  char *sequence;
  long int buffer_size;/*byte size fed to malloc*/
  long int sequence_size;
  int label_size;
};
  
struct filemanager *filemanager_init (char *filename);
void __filemanager_read_id (char *filename, char *id_buff);
seqfile_t __filemanager_get_filetype (struct filemanager *fmobj);
void filemanager_destroy  (struct filemanager *fmobj);

/*  if seq is not null reuse it and return it, otherwise allocate a new seq  */
struct sequence_t *__filemanager_next_seq (struct filemanager *fmobj, struct sequence_t *seq);
struct sequence_t *filemanager_next_seq (struct filemanager *fmobj, struct sequence_t *seq);

#define SINFO_SIZE 5000
struct sequences_info {
  unsigned int num_seqs;
  long int * offsets;
  long int * sizes;
};

struct sequences_info * filemanager_seq_count(struct filemanager *fmobj);
void sequences_info_free(struct sequences_info *);


#endif /* FILEMANAGER_H_ */
