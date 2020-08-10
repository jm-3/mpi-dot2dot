#include "utils.h"

void print_version (void) {
  printf ("Dot2dot-HPC version 0.9 build 23-07-2020\n");
}

void print_usage (void) {
  printf ("dot [options] in () the same option in the config file\n");
  printf("   -s, --sequence (Sequence): \n");
  printf("   sequence file name (accept fasta/multifasta/fastq).\n");
  printf("   The file format is detected inspecting the sequences’ header.\n");
  printf("   -c, --config:   [Command line only]\n");
  printf("   configuration file name (required)\n");
  printf("   -o, --output (Outfile):\n");
  printf("   output file name. (bed/dot)\n");
  printf("   -l, --minmotif (MinMotifLen):   [default=2]\n");
  printf("   minimum size in bp of the motif sequence \n");
  printf("   -L, --maxmotif (MaxMotifLen): [default=30]\n"); 
  printf("   maximum size in bp of the motif sequence\n");
  printf("   -m, --minmatch (MinMatch): [range (0,1) default=1]\n");
  printf("    Minimum overall matching score normalized in the range [0,1]\n");
  printf("   -G, --maxgaps (MaxGaps): [default=0]\n");
  printf("   maximum number of mismatches in a motif (expressed in bp).\n");
  printf("   -I, --maxinsert (MaxInsert): [default=0]\n");
  printf("   maximum insert size expressed in bp.\n");
  printf("   -S, --schedule (Schedule): [BLOCK | BALANCED. Default=BALANCED]\n");
  printf("   scheduling type to use.\n");  
  printf("   -v, --verbose: \n");
  printf("   print on the standard error the id of a sequence once loaded\n");
  printf("   -V, --version: \n");
  printf("   print dot-to-dot version and exit\n");
  printf("   -h, --help: \n");
  printf("   print a help page and exit\n");
  printf("\nConfig file options\n");
  printf("   FilterType: [default=NONE]\n");
  printf("   sets the level of filtering\n");
  printf("   * NONE: disable filtering\n");
  printf("   * THRESHOLD: filters-out only: too small TRs\n");
  printf("   * LIGHT: apply light filtering (see manual for details\n");
  printf("   * FAIR: apply fair filtering (see manual for details)\n");
  printf("   * HEAVY: apply heavy filtering (see manual for details)\n");
  printf("   MinTRLen: [default=12]\n");
  printf("   minimum length of a tandem repeat to be included in the output.\n");
  printf("   MinPurity: [range (0,1) default=0.1]\n");
  printf("   minimum purity of a tandem repeat to be included in the output.\n");
  printf("   Tolerance: [range (0,1) default=0]\n");
  printf("   Set a degree of tolerance to be accepted using heavy filtering\n");
  printf("   AllowOverlap: [range (Y/N) default=Y]\n");
  printf("   Enable/Disable  overlapping results in the final output.\n\n");
}

void * memalloc(size_t size, char *msg){
  void * aux;
  if((aux=malloc(size))==NULL){
    perror(msg);
    return NULL;
  }
  return aux;
}

void * memrealloc(void *ptr, size_t size, char *msg){
  void * aux;
  if((aux=realloc(ptr,size))==NULL){
    perror(msg);
    return NULL;
  }
  return aux;
}

double calcStdDevli(long int *values, int n){
  double sigma=0.0, sum=0.0, average;
  int i;

  for (i=0; i < n; i++)
    sum += values[i];
   
  average =  sum / (double) n;

  for(i=0; i<n; i++)
    sigma += (values[i] - average) * (values[i] - average);

  return sqrt(sigma/n);
}

long int getMinValueli(long int *values, int n){
  unsigned long int min;
  int i;

  min = values[0];

  for(i=1; i<n; i++)
    if(values[i] < min)
      min = values [i];

  return min;
}



long int getMaxValueli(long int *values, int n){
  long int max;
  int i;

  max = values[0];

  for(i=1; i<n; i++)
    if(values[i] > max)
      max = values [i];

  return max;
}

char* toLower(char* s) {
  for(char *p=s; *p; p++) *p=tolower(*p);
  return s;
}

const char *get_filename_ext(const char *filename) {
  const char *dot = strrchr(filename, '.');
  if(!dot || dot == filename) return "";
  return dot;
}

char *remove_filename_ext(char* myStr) {
  char *retStr;
  char *lastExt;
  if (myStr == NULL) return NULL;
  if ((retStr = malloc (strlen (myStr) + 1)) == NULL) return NULL;
  strcpy (retStr, myStr);
  lastExt = strrchr (retStr, '.');
  if (lastExt != NULL)
      *lastExt = '\0';
  return retStr;
}
