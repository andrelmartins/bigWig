#include <stdlib.h>
#include <stdio.h>
#include "bigwiglib.h"

int main(int argc, char ** argv) {
  char * chrom;
  int start;
  int end;
  int step;
  int do_abs = 0;
  char * filename;
  bigWig_t * bw;
  double * result;
  int len;
  int is_blank;
  int i;

  if (argc < 6) {
    printf("Usage: %s <bigWig file> <chrom> <start> <end> <step> [<do abs>=0]\n",
	   argv[0]);
    return EXIT_FAILURE;
  }

  filename = argv[1];
  chrom = argv[2];
  start = atoi(argv[3]);
  end = atoi(argv[4]);
  step = atoi(argv[5]);
  
  if (argc >= 7)
    do_abs = atoi(argv[6]);

  if (!is_bigwig(filename)) {
    fprintf(stderr, "Not a bigwig file: %s\n", filename);
    return EXIT_FAILURE;
  }

  bw = bigwig_load(filename, NULL);

  if (bw == NULL)
    return EXIT_FAILURE;

  if (bigwig_valid_chrom(bw, chrom) == 0) {
    fprintf(stderr, "Error: unknown chrom '%s'.\n", chrom);
    bigwig_free(bw);
    return EXIT_FAILURE;
  }

  result = bigwig_readf(bw, chrom, start, end, step, do_abs, &len, &is_blank);

  if (is_blank == 0) {
    for (i = 0; i < len; ++i)
      printf(" %g", result[i]);
    printf("\n");
  } else
    printf("no reads\n");

  bigwig_free(bw);
  free(result);

  return EXIT_SUCCESS;
}
