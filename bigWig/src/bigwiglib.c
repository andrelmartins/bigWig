#include "bw_query.h"
#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "options.h"
#include "udc.h"
#include "bigWig.h"
#include "obscure.h"
#include "errCatch.h"
#include "hmmstats.h"
#include "localmem.h"

int is_bigwig(char * filename) {
  return isBigWig(filename);
}

bigWig_t * bigwig_load(const char * filename, const char * udc_dir) {
  bigWig_t * bigwig = NULL;
  struct errCatch * err;

  /* set cache */
  if (udc_dir != NULL)
    udcSetDefaultDir((char*) udc_dir);

  /* setup error management & try to open file */
  err = errCatchNew();
  if (errCatchStart(err))
    bigwig = bigWigFileOpen((char*)filename);
  errCatchEnd(err);
  if (err->gotError) {
    fprintf(stderr, "error: %s\n", err->message->string);
    errCatchFree(&err);
    return NULL;
  }
  errCatchFree(&err);

  return bigwig;
}

void bigwig_free(bigWig_t * bw) {
  if (bw != NULL)
    bbiFileClose(&bw);
}


int bigwig_valid_chrom(bigWig_t * bw, const char * chrom) {
  return bw_has_chrom(bw, chrom);
}

long bigwig_chrom_size(bigWig_t * bw, const char * chrom) {
  return bw_chrom_size(bw, chrom);
}

double * bigwig_readf(bigWig_t * bw, const char * chrom, int start, int end, int step, int abs, int * out_length, int * out_is_blank) {
  int nIntervals;
  int size = bw_step_query_size(start, end, step);
  bwStepOp bwOp;
  double * result = (double*) calloc(size, sizeof(double));
  
  bw_select_op(&bwOp, "sum", 1);
  nIntervals = bw_step_query(bw, &bwOp, chrom, start, end, step, 0.0, abs, 0.0, result);

  if (nIntervals > 0)
    *out_is_blank = 0;
  else
    *out_is_blank = 1;

  *out_length = size;
  return result;
}

int * bigwig_readi(bigWig_t * bw, const char * chrom, int start, int end, int step, int abs, int * out_length, int * out_is_blank) {
  double * data = bigwig_readf(bw, chrom, start, end, step, abs, out_length, out_is_blank);
  int * result;
  int i, len;

  if (data == NULL)
    return NULL;
  
  len = *out_length;
  result = (int*) calloc(len, sizeof(int));

  if (*out_is_blank == 0) {
    for (i = 0; i < len; ++i)
      result[i] = (int) round(data[i]);
  }

  /* clean-up */
  free(data);

  return result;
}
