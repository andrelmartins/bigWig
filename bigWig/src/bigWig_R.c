/*
  C implementation of R interface

*/

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

#include "bwgExtra.h"

#include <string.h>
#include <R.h>
#include <Rdefines.h>

typedef struct bbiFile bigWig_t;

static SEXP wrap_int(int value) {
  SEXP result = NEW_INTEGER(1);
  INTEGER(result)[0] = value;
  return result;
}

static SEXP wrap_bool(int value) {
  SEXP result = NEW_LOGICAL(1);
  LOGICAL(result)[0] = value;
  return result;
}

static SEXP wrap_dbl(double value) {
  SEXP result = NEW_NUMERIC(1);
  REAL(result)[0] = value;
  return result;
}

static SEXP chrom_names(bigWig_t * bw) {
  SEXP result;
  struct bbiChromInfo * chrom, * chromList = bbiChromList(bw);
  result = allocVector(STRSXP, slCount(chromList));
  
  for (chrom = chromList; chrom != NULL; chrom = chrom->next)
    SET_STRING_ELT(result, chrom->id, mkChar(chrom->name));

  bbiChromInfoFreeList(&chromList);

  return result;
}

static SEXP chrom_sizes(bigWig_t * bw) {
  SEXP result;
  struct bbiChromInfo * chrom, * chromList = bbiChromList(bw);
  result = NEW_NUMERIC(slCount(chromList));
  
  for (chrom = chromList; chrom != NULL; chrom = chrom->next)
    REAL(result)[chrom->id] = chrom->size;

  bbiChromInfoFreeList(&chromList);

  return result;
}

static void bigwig_finalizer(SEXP ptr) {
  bigWig_t * bw;
  bw = R_ExternalPtrAddr(ptr);
  if (!bw) return;
  bbiFileClose(&bw);
  R_ClearExternalPtr(ptr); /* not really necessary */
}

SEXP bigWig_load(SEXP filename, SEXP udcDir) {
  SEXP ans, ans_names, ptr;
  bigWig_t * bigwig = NULL;
  struct errCatch * err;
  struct bbiSummaryElement sum;
  const char * cache = NULL;

  PROTECT(filename = AS_CHARACTER(filename));
  if (udcDir != R_NilValue) {
    PROTECT(udcDir = AS_CHARACTER(udcDir));
    cache = CHAR(STRING_ELT(udcDir, 0));
    if (cache != NULL)
      udcSetDefaultDir(strdup(cache)); /* this causes a memory leak
					  but it's a small one ... */
    UNPROTECT(1);
  }

  /* load bigWig */
  err = errCatchNew();
  if (errCatchStart(err))
    bigwig = bigWigFileOpen((char*) CHAR(STRING_ELT(filename, 0)));
  errCatchEnd(err);
  if (err->gotError) {
    /* TODO: Add some resource clean-up because Kent code doesn't
             clean up after itself on errors.
    */
    REprintf("error: %s\n", err->message->string);
    errCatchFree(&err);
    UNPROTECT(1);
    return R_NilValue;
  }
  errCatchFree(&err);

  PROTECT(ans = allocVector(VECSXP, 13));
  PROTECT(ans_names = allocVector(STRSXP, 13));

  /* make external pointer */
  ptr = R_MakeExternalPtr(bigwig, install("BIGWIG_struct"), R_NilValue);
  PROTECT(ptr);
  R_RegisterCFinalizerEx(ptr, bigwig_finalizer, TRUE);
  setAttrib(ans, install("handle_ptr"), ptr);

  /* fill info */
  sum = bbiTotalSummary(bigwig);
  SET_VECTOR_ELT(ans, 0, wrap_int(bigwig->version));
  SET_STRING_ELT(ans_names, 0, mkChar("version"));
  SET_VECTOR_ELT(ans, 1, wrap_bool(bigwig->uncompressBufSize > 0));
  SET_STRING_ELT(ans_names, 1, mkChar("isCompressed"));
  SET_VECTOR_ELT(ans, 2, wrap_bool(bigwig->isSwapped));
  SET_STRING_ELT(ans_names, 2, mkChar("isSwapped"));
  SET_VECTOR_ELT(ans, 3, wrap_dbl((double) (bigwig->unzoomedIndexOffset - bigwig->unzoomedDataOffset))); /* R doesn't support long int, so we'll use doubles
													    see: http://r.789695.n4.nabble.com/long-integer-in-R-td1478798.html
													 */
  SET_STRING_ELT(ans_names, 3, mkChar("primaryDataSize"));
  if (bigwig->levelList != NULL) {
    long long indexEnd = bigwig->levelList->dataOffset;
    SET_VECTOR_ELT(ans, 4, wrap_dbl(indexEnd - bigwig->unzoomedIndexOffset));
  } else
    SET_VECTOR_ELT(ans, 4, R_NilValue);
  SET_STRING_ELT(ans_names, 4, mkChar("primaryIndexSize"));

  SET_VECTOR_ELT(ans, 5, wrap_int(bigwig->zoomLevels));
  SET_STRING_ELT(ans_names, 5, mkChar("zoomLevels"));
  SET_VECTOR_ELT(ans, 6, chrom_names(bigwig));
  SET_STRING_ELT(ans_names, 6, mkChar("chroms"));
  SET_VECTOR_ELT(ans, 7, chrom_sizes(bigwig));
  SET_STRING_ELT(ans_names, 7, mkChar("chromSizes"));
  SET_VECTOR_ELT(ans, 8, wrap_dbl(sum.validCount));
  SET_STRING_ELT(ans_names, 8, mkChar("basesCovered"));
  SET_VECTOR_ELT(ans, 9, wrap_dbl(sum.sumData / sum.validCount));
  SET_STRING_ELT(ans_names, 9, mkChar("mean"));
  SET_VECTOR_ELT(ans, 10, wrap_dbl(sum.minVal));
  SET_STRING_ELT(ans_names, 10, mkChar("min"));
  SET_VECTOR_ELT(ans, 11, wrap_dbl(sum.maxVal));
  SET_STRING_ELT(ans_names, 11, mkChar("max"));
  SET_VECTOR_ELT(ans, 12, wrap_dbl(calcStdFromSums(sum.sumData, sum.sumSquares, sum.validCount)));
  SET_STRING_ELT(ans_names, 12, mkChar("std"));

  setAttrib(ans, R_NamesSymbol, ans_names);

  UNPROTECT(4);

  return ans;
}

SEXP bigWig_unload(SEXP obj) {
  SEXP ptr;

  PROTECT(ptr = GET_ATTR(obj, install("handle_ptr")));
  if (ptr == R_NilValue)
    error("invalid bigwig object");

  bigwig_finalizer(ptr);

  UNPROTECT(1);

  return R_NilValue;
}

SEXP bigWig_query(SEXP obj, SEXP chrom, SEXP start, SEXP end, SEXP clip) {
  SEXP ptr, res = R_NilValue;
  bigWig_t * bigwig;
  int do_clip;

  PROTECT(clip = AS_LOGICAL(clip));
  do_clip = LOGICAL(clip)[0] == TRUE;
  PROTECT(chrom = AS_CHARACTER(chrom));
  PROTECT(start = AS_INTEGER(start));
  PROTECT(end = AS_INTEGER(end));
  PROTECT(ptr = GET_ATTR(obj, install("handle_ptr")));
  if (ptr == R_NilValue)
    error("invalid bigWig object");

  bigwig = R_ExternalPtrAddr(ptr);
  if (bigwig == NULL) {
    error("bigWig object has been unloaded");
  } else {
    struct lm * localMem = lmInit(0); /* use default value */
    struct bbiInterval * intervals;
    
    if (do_clip)
      intervals = bigWigIntervalQuery(bigwig,
			    (char*) CHAR(STRING_ELT(chrom, 0)),
			    INTEGER(start)[0],
			    INTEGER(end)[0],
			    localMem);
    else
      intervals = bigWigIntervalQueryNoClip(bigwig,
                                      (char*) CHAR(STRING_ELT(chrom, 0)),
                                      INTEGER(start)[0],
                                      INTEGER(end)[0],
                                      localMem);

    /* convert result into an R matrix */
    int nIntervals = slCount(intervals);
    
    if (nIntervals > 0) {
      int nx, ny;
      double * mptr;
      struct bbiInterval * interval;  
      int i;
      nx = nIntervals;
      ny = 3;
      PROTECT(res = allocMatrix(REALSXP, nx, ny));
      mptr = REAL(res);

      for (i = 0, interval = intervals; 
           interval != NULL;
           interval = interval->next, i++) {
        mptr[i] = (double) interval->start;
        mptr[i + nx] = (double) interval->end;
        mptr[i + nx*2] = interval->val;
      }

      UNPROTECT(1);
    }

    /* clean-up */
    lmCleanup(&localMem);
  }

  UNPROTECT(5);

  return res;
}
