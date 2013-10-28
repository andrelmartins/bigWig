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

/*
  Query version to speed up meta-plots
*/
SEXP bigWig_query_by_step(SEXP obj, SEXP chrom, SEXP start, SEXP end, SEXP step, SEXP doSum, SEXP gapValue) {
  SEXP ptr, res = R_NilValue;
  bigWig_t * bigwig;
  int do_sum = 0;

  PROTECT(doSum = AS_LOGICAL(doSum));
  if (LOGICAL(doSum)[0] == TRUE)
    do_sum = 1;

  PROTECT(chrom = AS_CHARACTER(chrom));
  PROTECT(start = AS_INTEGER(start));
  PROTECT(end = AS_INTEGER(end));
  PROTECT(step = AS_INTEGER(step));
  PROTECT(gapValue = AS_NUMERIC(gapValue));
  PROTECT(ptr = GET_ATTR(obj, install("handle_ptr")));
  if (ptr == R_NilValue)
    error("invalid bigWig object");

  bigwig = R_ExternalPtrAddr(ptr);
  if (bigwig == NULL) {
    error("bigWig object has been unloaded");
  } else {
    int istart = INTEGER(start)[0];
    int iend = INTEGER(end)[0];
    int istep = INTEGER(step)[0];
    struct lm * localMem = lmInit(0); /* use default value */
    struct bbiInterval * intervals
      = bigWigIntervalQuery(bigwig,
                          (char*) CHAR(STRING_ELT(chrom, 0)),
                          istart,
                          iend,
                          localMem);

    /* convert result into an R matrix */
    int nIntervals = slCount(intervals);
    
    if (nIntervals > 0) {
      int size = (iend - istart)/istep;
      struct bbiInterval * interval;
      int count = 0;
      double sum = 0.0;
      double left, right;
      int idx;
      double d_gap_value = REAL(gapValue)[0];

      PROTECT(res = NEW_NUMERIC(size));

      /* init to 'gapValue' */
      for (idx = 0; idx < size; ++idx)
        REAL(res)[idx] = d_gap_value;

      left = istart;
      right = istart + istep;
      idx = 0;
      for (interval = intervals; interval != NULL && idx < size; interval = interval->next) {
        /* interval starts beyond current step */
        if (((double)interval->start) >= right) {
          /* save current value */
          if (count > 0) {
            if (do_sum == 1)
              REAL(res)[idx] = sum;
            else
              REAL(res)[idx] = sum / count;
          }
          
          count = 0;
          sum = 0.0;
          while (idx < size && ((double)interval->start) >= right) {
            ++idx;
            left += istep;
            right += istep;
          }
        }
        
        /* interval starts at or before the current step */
        if (((double)interval->start) < right) {
          sum += interval->val;
          ++count;
        }
        
        /* interval ends beyond or atthe current step */
        if (((double)interval->end) >= right && idx < size) {
          do {
            /* save current step */
            if (count > 0) {
              if (do_sum == 1)
                REAL(res)[idx] = sum;
              else
                REAL(res)[idx] = sum/count;
            }
            
            ++idx;
            left += istep;
            right += istep;
            
            count = 1;
            sum = interval->val;
          } while (((double)interval->end) >= right && idx < size);
          if (left >= ((double)interval->end)) {
            count = 0; /* current interval now ends before left */
            sum = 0.0;
          }
        }
      }
      
      if (count > 0 && idx < size) {
        if (do_sum == 1)
          REAL(res)[idx] = sum;
        else
          REAL(res)[idx] = sum/count;
      }
      
      UNPROTECT(1);
    }
    
    /* clean-up */
    lmCleanup(&localMem);
  }
  
  UNPROTECT(7);
  
  return res;
}

#define MODE_MAX 0
#define MODE_MIN 1
#define MODE_SUM 2
#define MODE_AVG 3

static int aggregator_to_mode(const char * aggregator) {
  if (!strcmp("max", aggregator))
    return MODE_MAX;
  if (!strcmp("min", aggregator))
    return MODE_MIN;
  if (!strcmp("sum", aggregator))
    return MODE_SUM;
  if (!strcmp("mean", aggregator))
    return MODE_AVG;

  error("unknown aggregator value: %s", aggregator);
  return -1;
}

SEXP bigWig_bed_query(SEXP bed, SEXP bwPlus, SEXP bwMinus, SEXP gapValue, SEXP weighted, SEXP aggregator) {
  SEXP ptr, res = R_NilValue;
  bigWig_t * bigwig_plus, * bigwig_minus;

  double gap_value = 0;
  int has_gapvalue = 0;
  int is_weighted = 0;
  int mode = MODE_SUM;
  int protect_count = 0;
  int i, N;

  SEXP chroms, starts, ends, strands = R_NilValue;
  int has_strand = 0;

  struct lm * localMem;
  struct bbiInterval * intervals;

  //PROTECT(bed = AS_LIST(chrom)); /* a data frame is just a list :) */
  PROTECT(bed);
  
  chroms = VECTOR_ELT(bed, 0);
  if (!isFactor(chroms) && TYPEOF(chroms) != STRSXP)
    error("first column of bed file must be a factor or character vector");

  PROTECT(starts = AS_INTEGER(VECTOR_ELT(bed, 1)));
  PROTECT(ends = AS_INTEGER(VECTOR_ELT(bed, 2)));

  if (length(bed) >= 6) {
    has_strand = 1;
    strands = VECTOR_ELT(bed, 5);

    if (!isFactor(strands) && TYPEOF(strands) != STRSXP)
      error("sixth column of bed file must be a factor or character vector");
  }

  /* options */
  if (gapValue != R_NilValue) {
    has_gapvalue = 1;
    PROTECT(gapValue = AS_NUMERIC(gapValue));
    gap_value = REAL(gapValue)[0];
    UNPROTECT(1);
  }

  PROTECT(weighted = AS_LOGICAL(weighted));
  if (LOGICAL(weighted)[0] == TRUE)
    is_weighted = 1;

  PROTECT(aggregator = AS_CHARACTER(aggregator));
  mode = aggregator_to_mode(CHAR(STRING_ELT(aggregator, 0)));


  /* Plus is mandatory */
  PROTECT(ptr = GET_ATTR(bwPlus, install("handle_ptr")));
  if (ptr == R_NilValue)
    error("invalid bigWig object");

  protect_count = 6;

  /* bigwig files */
  bigwig_plus = R_ExternalPtrAddr(ptr);
  if (bigwig_plus == NULL)
    error("bigWig plus object has been unloaded");
  
  /* Minus is optional */
  if (bwMinus == R_NilValue) {
    bigwig_minus = NULL;
    if (has_strand == 1) {
      has_strand = 0;
      warning("bed has strand column but no 'minus' bigWig supplied");
    }
  } else {
    PROTECT(ptr = GET_ATTR(bwMinus, install("handle_ptr")));
    if (ptr == R_NilValue)
      error("invalid bigWig object");

    /* bigwig files */
    bigwig_minus = R_ExternalPtrAddr(ptr);
    if (bigwig_minus == NULL)
      error("bigWig minus object has been unloaded");

    ++protect_count;
  }

  /* parse inputs */
  N = length(chroms);
  PROTECT(res = NEW_NUMERIC(N));
  ++protect_count;
  
  localMem = lmInit(0); /* use default value */

  for (i = 0; i < N; ++i) {
    /* collect chrom, start, end, strand (if any) */
    const char * chrom;
    char strand = '+';
    int start = INTEGER(starts)[i];
    int end = INTEGER(ends)[i];
    bigWig_t * bigwig;
    int nIntervals;

    if (isFactor(chroms)) {
      int idx = INTEGER(chroms)[i] - 1;
      chrom = CHAR(STRING_ELT(GET_LEVELS(chroms), idx));
    } else
      chrom = CHAR(STRING_ELT(chroms, i));

    if (has_strand == 1) {
      if (isFactor(strands)) {
        int idx = INTEGER(strands)[i] - 1;
        strand = CHAR(STRING_ELT(GET_LEVELS(strands), idx))[0];
      } else
        strand = CHAR(STRING_ELT(strands, i))[0];
    }

    if (strand == '+')
      bigwig = bigwig_plus;
    else
      bigwig = bigwig_minus;

    /* do query */
    intervals = bigWigIntervalQuery(bigwig, (char*)chrom, start, end, localMem);
    nIntervals = slCount(intervals);

    /* walk over data (depending on mode, weighted) */
    if (nIntervals > 0) {
      double accum;
      double weight = 1;
      double covered = 0;
      struct bbiInterval * interval = intervals;
      
      /* get first value */
      if (mode == MODE_MIN && mode == MODE_MAX)
        accum = interval->val;
      else {
        if (mode == MODE_AVG && is_weighted == 1)
          weight = (double) (interval->end - interval->start);
        covered = (interval->end - interval->start);
        
        accum = weight * interval->val;
      }
      
      for (interval = interval->next; interval != NULL; interval = interval->next) {
        if (mode == MODE_MIN) {
          if (accum > interval->val)
            accum = interval->val;
        } else if (mode == MODE_MAX) {
          if (accum < interval->val)
            accum = interval->val;
        } else {
          if (mode == MODE_AVG && is_weighted == 1)
            weight = (double) (interval->end - interval->start);
          covered += (interval->end - interval->start);
          
          accum += weight * interval->val;
        }
      }
      
      /* add gap value if weighted */
      if (mode == MODE_AVG && has_gapvalue == 1 && is_weighted == 1)
        accum += ((end - start) - covered) * gap_value;
      
      /* store result */
      if (mode == MODE_AVG) {
        if (is_weighted == 1)
          REAL(res)[i] = accum / ((double) (end - start));
        else
          REAL(res)[i] = accum / nIntervals;
      } else
        REAL(res)[i] = accum;
      
    } else if (has_gapvalue == 1)
      REAL(res)[i] = gap_value;
    else
      REAL(res)[i] = NA_REAL;
    
    /* reclaim memory */
    if (i % 1000) {
      lmCleanup(&localMem);
      localMem = lmInit(0); /* use default value */
    }
  }

  lmCleanup(&localMem);

  /* clean up */
  UNPROTECT(protect_count);

  return res;
}

SEXP bigWig_chrom_step_sum(SEXP obj, SEXP chrom, SEXP step, SEXP defaultValue) {
  SEXP ptr, res = R_NilValue;
  bigWig_t * bigwig;

  PROTECT(chrom = AS_CHARACTER(chrom));
  PROTECT(step = AS_INTEGER(step));
  PROTECT(defaultValue = AS_NUMERIC(defaultValue));
  PROTECT(ptr = GET_ATTR(obj, install("handle_ptr")));
  if (ptr == R_NilValue)
    error("invalid bigWig object");

  bigwig = R_ExternalPtrAddr(ptr);
  if (bigwig == NULL) {
    error("bigWig object has been unloaded");
  } else {
    const char * cchrom = CHAR(STRING_ELT(chrom, 0));
    int istep = INTEGER(step)[0];
    double defval = REAL(defaultValue)[0];
    int n;
    int i, j, k;
    double * valptr;

    struct bigWigValsOnChrom *chromVals = bigWigValsOnChromNew();

    if (!bigWigValsOnChromFetchData(chromVals, (char*)cchrom, bigwig))
      error("could not retrieve information on chrom: %s", chrom);

    n = chromVals->chromSize / istep;
    PROTECT(res = NEW_NUMERIC(n));

    for (i = 0, j = 0, valptr=chromVals->valBuf; i < n; ++i) {
      double sum = 0;
      int count = 0;
      
      for (k = 0; k < istep; ++k, ++j, ++valptr) {
        if (bitReadOne(chromVals->covBuf , j)) {
          ++count;
          sum += *valptr;
        }
      }
      if (count > 0)
        REAL(res)[i] = sum;
      else
        REAL(res)[i] = defval;
    }
    
    UNPROTECT(1);
    bigWigValsOnChromFree(&chromVals);
  }
  
  UNPROTECT(4);

  return res;
}
