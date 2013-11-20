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
#include "bw_query.h"

#include <string.h>
#include <R.h>
#include <Rdefines.h>

/* Interface to abstract handling of fragmented bigWig objects
 * 
 * "obj" is either an R object of class 'bigWig' with an 'handle_ptr'
 * attribute (classic bigWig object) or it is a character vector 
 * containing the prefix and suffix to the path of each
 * bigWig fragment (path = <prefix><chrom><suffix>).
 */
bigWig_t * bigWig_for_chrom(SEXP obj, const char * chrom) {
  SEXP ptr;
  bigWig_t * bigWig = NULL;

  if (IS_CHARACTER(obj)) {
    struct errCatch * err;
    SEXP prefix, suffix;
    int slen;
    char * path;
    
    if (Rf_length(obj) != 2)
      error("bigWig fragment must be a set of one prefix/suffix pair");
    
    prefix = STRING_ELT(obj, 0);
    suffix = STRING_ELT(obj, 1);
    
    /* create path string */
    slen = strlen(chrom) + strlen(CHAR(prefix)) + strlen(CHAR(suffix)) + 1;
    path = R_alloc(slen, sizeof(char));
    sprintf(path, "%s%s%s", CHAR(prefix), chrom, CHAR(suffix));
    
    /* open file */
    err = errCatchNew();
    if (errCatchStart(err))
      bigWig = bigWigFileOpen(path);
    errCatchEnd(err);
    
    if (err->gotError) {
      /* TODO: Add some resource clean-up because Kent code doesn't
               clean up after itself on errors.
      */
      REprintf("error: %s\n", err->message->string);
      errCatchFree(&err);
      UNPROTECT(1);
      return NULL;
    }
    errCatchFree(&err);
    
    UNPROTECT(1);
    return bigWig;
  }
  
  PROTECT(ptr = GET_ATTR(obj, install("handle_ptr")));
  if (ptr == R_NilValue)
    error("invalid bigWig object");

  bigWig = R_ExternalPtrAddr(ptr);
  if (bigWig == NULL) {
    error("bigWig object has been unloaded");
  }
  
  UNPROTECT(1);
  return bigWig;
}

void bigWig_for_chrom_release(SEXP obj, bigWig_t * bw) {
  if (IS_CHARACTER(obj)) {
    /* locally loaded bigWig object, needs to be released */
    bbiFileClose(&bw);
  } /* else nothing to do */
}


/*
 *  Query and auxiliary functions  
 */

const char * char_elt(SEXP vec, int idx) {
  if (isFactor(vec)) {
    int f_idx = INTEGER(vec)[idx] - 1;
    return CHAR(STRING_ELT(GET_LEVELS(vec), f_idx));
  } else
    return CHAR(STRING_ELT(vec, idx));
}

/* in-place reverse of real vector */
void vec_reverse(SEXP rvec) {
  int len = length(rvec);
  int i, j;
  double * ptr = REAL(rvec);
  
  for (i = 0, j = len - 1; i < j; ++i, --j) {
    double tmp = ptr[i];
    ptr[i] = ptr[j];
    ptr[j] = tmp;
  }
}

void fill_row(SEXP matrix, SEXP row, int row_idx) {
  if (!isMatrix(matrix)) {
    REAL(matrix)[row_idx] = REAL(row)[0];
  } else {
    int nx = nrows(matrix);
    int ny = ncols(matrix);
    int i, k;
    double * m_ptr = REAL(matrix);
    double * r_ptr = REAL(row);

    /* R matrices are stored column-wise */
  
    for (i = 0, k = row_idx; i < ny; ++i, k += nx)
      m_ptr[k] = r_ptr[i];
  }
}

SEXP bigWig_probe_query(SEXP obj_plus, SEXP obj_minus, SEXP bed, SEXP op, SEXP step, SEXP use_strand, SEXP with_attributes, SEXP as_matrix, SEXP gap_value, SEXP abs_value) {
  bwStepOp bwOp;
  bigWig_t * bw = NULL;
  SEXP result;
  int i, N;
  const char * prev_chrom = NULL;
  int prev_strand_plus = 1;
  int istep;
  double d_gap;
  int do_abs = INTEGER(abs_value)[0] == TRUE;
  int is_matrix = INTEGER(as_matrix)[0] == TRUE;
  int no_step = 0;
  
  PROTECT(step = AS_INTEGER(step));
  istep = INTEGER(step)[0];
  no_step = (istep == NA_INTEGER);
  
  if (!no_step && istep <= 0)
    error("step size must be >= 1");
  
  PROTECT(gap_value = AS_NUMERIC(gap_value));
  d_gap = REAL(gap_value)[0];
  
  // for each bed ...
  SEXP chroms, starts, ends, strands = R_NilValue;
  int has_strand = 0;
  
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
  
  // initialize selected operation
  bw_select_op(&bwOp, CHAR(STRING_ELT(op, 0)), 1);
  
  // decide if output is list or matrix (assume that R side has validated the option)
  N = length(chroms);
  if (is_matrix) {
    int start = INTEGER(starts)[0];
    int end = INTEGER(ends)[0];
    
    if (no_step)
      PROTECT(result = NEW_NUMERIC(N));
    else {
      int size = (end - start)/istep;
      PROTECT(result = allocMatrix(REALSXP, N, size));
    }
  } else {
    PROTECT(result = NEW_LIST(N));
  }
  
  // process bed file (foreach.bed ...)
  for (i = 0; i < N; ++i) {
    SEXP res;
    const char * chrom = char_elt(chroms, i);
    int is_plus = 1;
    int start, end;
    
    if (has_strand) {
      const char * strand = char_elt(strands, i);
      
      if (!strcmp(strand, "-"))
        is_plus = 0;
    }
    
    if (prev_chrom == NULL || strcmp(chrom, prev_chrom) || prev_strand_plus != is_plus) {
      if (bw != NULL) {
        if (is_plus || obj_minus == R_NilValue)
          bigWig_for_chrom_release(obj_plus, bw);
        else
          bigWig_for_chrom_release(obj_minus, bw);
      }
      
      // get bigWig object
      if (is_plus || obj_minus == R_NilValue)
        bw = bigWig_for_chrom(obj_plus, chrom);
      else
        bw = bigWig_for_chrom(obj_minus, chrom);
      
      // update
      prev_chrom = chrom;
      prev_strand_plus = is_plus;
    } // else keep the same object
    
    start = INTEGER(starts)[i];
    end = INTEGER(ends)[i];
    if (no_step)
      istep = end - start;
    
    PROTECT(res = bw_step_query(bw, &bwOp, chrom, start, end, istep, d_gap, do_abs));
  
    // reverse if on negative strand
    if (is_plus == 0)
      vec_reverse(res);
    
    if (is_matrix) {
      fill_row(result, res, i);
      // release memory ??
    } else {
      SET_VECTOR_ELT(result, i, res);
      // if use attributes && !as_matrix => add attributes to vector
      // ...
    }
    
    UNPROTECT(1);
  }
  
  // if (use attributes && as_matrix => add step size to matrix)
  // ...
  
  
  
  // clean up
  if (bw != NULL) {
    if (prev_strand_plus || obj_minus == R_NilValue)
      bigWig_for_chrom_release(obj_plus, bw);
    else
      bigWig_for_chrom_release(obj_minus, bw);
  }

  UNPROTECT(5);
  
  return result;
}
