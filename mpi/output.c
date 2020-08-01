#include "output.h"
#include "common.h"

struct outfile *output_create (char *filename) {
  struct outfile *out;
  char buff[MAX_PARAM_LEN];
  int len;
  out = (struct outfile *) malloc (sizeof (struct outfile));
  /* In case of memory faiulure just print the dot format on output  */
  if (out == NULL) return NULL;
  if (filename == NULL) {
    out->pf = NULL;
    out->filetype = dot;
    return out;
  }
  len = strlen (filename);
  if (len < 4) {  /*  cannot have an extension  */
    sprintf (buff, "%s.dot", filename);
    out->pf = fopen (buff, "w");
    out->filetype = dot;
    return out;
  }
  if (strcmp (filename + len - 4, ".bed") == 0) {
    out->pf = fopen (filename, "w");
    out->filetype = bed;
    return out;
  }
  if (strcmp (filename + len - 4, ".dot") == 0) {
    out->pf = fopen (filename, "w");
    out->filetype = dot;
    return out;
  }
  sprintf (buff, "%s.dot", filename);
  out->pf = fopen (buff, "w");
  out->filetype = dot;
  return out;
}

struct outfile *output_create_tmp (char *output_filename, int rank){
  char tempOuputFile[MAX_PARAM_LEN+50], *path, *filename, *aux;
  char ext[20];

  aux = strdup(output_filename);
  path = dirname(aux);
  aux = strdup(output_filename);
  filename = basename(aux);
  
  strcpy (ext, get_filename_ext(filename));
  filename = remove_filename_ext(filename);
  if(strcmp(ext, "") == 0)
    strcpy(ext, ".dot");

  sprintf(tempOuputFile,"%s/tmp-%s-%d%s", path, filename, rank, ext);
  
  return output_create (tempOuputFile);
}

void output_destroy (struct outfile *out) {
  if (out == NULL) return;
  if (out->pf != NULL) fclose (out->pf);
  free (out);
  out = NULL;
}

void print_header (struct outfile *out) {
  if (out == NULL) {
    printf ("#ID\tStart\tEnd\tCN\tLen\tPurity\tSequence\n");
    return;
  }
  if (out->filetype == dot && out->pf != NULL) {
    fprintf (out->pf, "#ID\tStart\tEnd\tCN\tLen\tPurity\tSequence\n");
    return;
  }
  if (out->filetype == dot && out->pf == NULL) {
    printf ("#ID\tStart\tEnd\tCN\tLen\tPurity\tSequence\n");
    return;
  }
  if (out->filetype == bed && out->pf != NULL) {
    fprintf (out->pf, "#ID\tStart\tEnd\n");
    return;
  }
  if (out->filetype == bed && out->pf == NULL) {
    printf ("#ID\tStart\tEnd\n");
    return;
  }
}

void print_TRs_list_toFile (struct outfile *out, char* label, char* sequence, TRs_Result_Bundle* trs_bundle) {
  TRs_result_t* trs_toPrint = trs_bundle->TRs_found;
  unsigned short int* motifs_lengths = trs_bundle->motif_lengths;
  unsigned int i;
  char buff[32769], buff2[32768];
  /*Print the record of tandem repeats to FILENAME*/
  i=0;
  while ( i < trs_bundle->trs_found_offset ) {
    /*  Print only unfiltered TRs */
    if (trs_toPrint[i].valid_TR == true) {
      dump_TR_elem (out, sequence, &trs_toPrint[i], motifs_lengths, buff, buff2, 32769);
      if (out->pf == NULL) printf ("%s\t%s", label ,buff);
      else fprintf (out->pf, "%s\t%s", label ,buff);
    }  
    ++i;
  }
}

void dump_TR_elem (struct outfile *out, char* sequence, TRs_result_t* tr, unsigned short int* motif_lengths, char *outStr, char *tmpbuff, size_t n) {

  unsigned short int j, tmpbuff_index;
  short int motifs_number;
  long int seq_index, i;
  int chars_number;

  if (out->filetype == bed) {
    snprintf (outStr, n, "%ld\t%ld\n", tr->origin_position, tr->origin_position + tr->full_length);
    return;
  }
  chars_number = snprintf (outStr, n, "%ld\t%ld\t%d\t%d\t%.3f\t",
			   tr->origin_position+1, tr->origin_position + tr->full_length,
			   tr->copy_number,tr->period, tr->purity_percentage);

#ifdef DEBUG_OUTPUT
  /*   These are for debugging */
  strcpy (tmpbuff, outStr);
  if (tr->valid_TR == true) 
    snprintf (outStr, n, "%s%.3f\tok\t", tmpbuff, tr->stats);  
  else  
    snprintf (outStr, n, "%s%.3f\tno\t", tmpbuff, tr->stats);  
#endif
	/* rifaccio lo switch per tenere in considerazione il blocco if DEBUG precedente*/

  strcpy (tmpbuff, outStr);

  /* Build the tandem repeat from the sequence  */
  seq_index = tr->origin_position;
  tmpbuff_index = chars_number;
  motifs_number = tr->motifs_number;
  for (i = tr->motif_start_index; motifs_number > 0 ; ++i, --motifs_number) {
    j=0;
    while ( j < motif_lengths[i] ) {
      tmpbuff[tmpbuff_index] = sequence[seq_index];
      ++j;
      ++seq_index;
      ++tmpbuff_index;
    }
    tmpbuff[tmpbuff_index] = ' ';
    ++tmpbuff_index;
  }
  tmpbuff[tmpbuff_index - 1] = '\0';

  snprintf (outStr, n, "%s\n", tmpbuff);
  /*print_usintArray(motif_lengths, tr->motifs_number);*/
  /*strcpy (tmpbuff, outStr);  */
  tmpbuff[0] = '\0';
}

void print_usintArray(unsigned short int* array, int length) {
	int i;
	for(i=0; i<length; i++) {
		printf("%hu ", array[i]);
	}
	printf("\n");
}

int merge_files(const char* fullpath, int n) {
  int reads, writes, len;
  char buffer[BUFF_SIZE];
  char outputFile[MAX_PARAM_LEN+50], tempOuputFile[MAX_PARAM_LEN+50], *path, *filename, *aux;
  char ext[20];

  aux = strdup(fullpath);
  path = dirname(aux);
  aux = strdup(fullpath);
  filename = basename(aux);
  
  strcpy (ext, get_filename_ext(filename));
  filename = remove_filename_ext(filename);
  if(strcmp(ext, "") == 0)
    strcpy(ext, ".dot");

  len = strlen (fullpath);
  if (len < 4) {  /*  cannot have an extension  */
    sprintf (outputFile, "%s.dot", fullpath);    
  }else
    if (strcmp (fullpath + len - 4, ".bed") == 0 || strcmp (fullpath + len - 4, ".dot") == 0) {
      sprintf (outputFile, "%s", fullpath);    
    }else
      sprintf (outputFile, "%s.dot", fullpath);
  

  FILE* fp = fopen(outputFile, "w");
  if (fp == NULL) return 1;

  for (int i = 0; i < n; i++) {
    sprintf(tempOuputFile,"%s/tmp-%s-%d%s",path,filename, i, ext);

    FILE* fp_read = fopen(tempOuputFile, "r");
    if (fp_read == NULL) return 1;

    while ((reads = fread(buffer, sizeof(char), BUFF_SIZE, fp_read))) {
        writes = fwrite(buffer, sizeof(char), reads, fp);
        if (reads != writes) return 1;
    }

    fclose(fp_read);
 }

 fclose(fp);

 return 0;
}

int remove_tmpfiles(const char* fullpath, int n){
  char tempOuputFile[MAX_PARAM_LEN+70], *path, *filename, *aux;
  char ext[20];

  aux = strdup(fullpath);
  path = dirname(aux);
  aux = strdup(fullpath);
  filename = basename(aux);
  
  strcpy (ext, get_filename_ext(filename));
  filename = remove_filename_ext(filename);
  if(strcmp(ext, "") == 0)
    strcpy(ext, ".dot");

  for (int i = 0; i < n; i++) {
    sprintf(tempOuputFile,"%s/tmp-%s-%d%s", path, filename, i, ext);
       
    if(remove(tempOuputFile)) return 1;    
 }
 return 0;
}

