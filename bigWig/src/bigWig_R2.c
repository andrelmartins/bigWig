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
      
      /* copy error string into R object (will be freed by GC later) */
      SEXP errStr = mkChar(err->message->string);
      
      errCatchFree(&err);
      
      error("error: %s", CHAR(errStr)); /* use string to report error back to R */
      return NULL;
    }
    errCatchFree(&err);
    
    /* check if chromsome exists in file */
    if (bw_has_chrom(bigWig, chrom) == 0) {
      bbiFileClose(&bigWig);
      
      error("file '%s' has no information on chrom '%s'", path, chrom);
      return NULL;
    }
    
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
  
  if (bw_has_chrom(bigWig, chrom) == 0) {
    error("bigWig has no information on chromosome: '%s'", chrom);
  }
  
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

static SEXP R_bw_step_query(bigWig_t * bigwig, bwStepOp * op, const char * chrom, int start, int end, int step, double gap_value, int do_abs, double thresh) {
  int size = bw_step_query_size(start, end, step);
  SEXP result = NEW_NUMERIC(size);
  bw_step_query(bigwig, op, chrom, start, end, step, gap_value, do_abs, thresh, REAL(result));
  
  return result;
}

typedef SEXP (*step_query_func)(bigWig_t * bigwig, bwStepOp * op, const char * chrom, int start, int end, int step, double gap_value, int do_abs, int is_plus, void * uptr);

SEXP bigWig_region_query(SEXP obj_plus, SEXP obj_minus, SEXP bed, bwStepOp bwOp, SEXP step, SEXP use_strand, SEXP with_attributes, SEXP as_matrix, SEXP gap_value, SEXP abs_value, SEXP follow_strand, step_query_func step_query, void * uptr) {
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
  int use_attributes = INTEGER(with_attributes)[0] == TRUE;
  int do_rev_strand = INTEGER(follow_strand)[0] == TRUE;
  
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

    if (step_query != NULL)
      PROTECT(res = (*step_query)(bw, &bwOp, chrom, start, end, istep, d_gap, do_abs, is_plus, uptr));
    else
      PROTECT(res = R_bw_step_query(bw, &bwOp, chrom, start, end, istep, d_gap, do_abs, 0.0));
  
    // reverse if on negative strand
    if (is_plus == 0 && do_rev_strand)
      vec_reverse(res);
    
    if (is_matrix) {
      fill_row(result, res, i);
    } else {
      if (use_attributes) {
        SEXP att_chrom = NEW_STRING(1);
        SEXP att_start = NEW_INTEGER(1);
        SEXP att_end = NEW_INTEGER(1);
        SEXP att_step = NEW_INTEGER(1);
        
        SET_STRING_ELT(att_chrom, 0, mkChar(chrom));
        INTEGER(att_start)[0] = start;
        INTEGER(att_end)[0] = start + Rf_length(res) * istep;
        INTEGER(att_step)[0] = istep;
        
        setAttrib(res, install("chrom"), att_chrom);
        setAttrib(res, install("start"), att_start);
        setAttrib(res, install("end"), att_end);
        setAttrib(res, install("step"), att_step);

      }      
      SET_VECTOR_ELT(result, i, res);
    }
    
    UNPROTECT(1);
  }
  
  // if (use attributes && as_matrix => add step size to matrix)
  if (use_attributes && is_matrix) {
    SEXP step_size = NEW_INTEGER(1);
    INTEGER(step_size)[0] = istep;
    
    setAttrib(result, install("step"), step_size);
  }
  
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

SEXP bigWig_probe_query(SEXP obj_plus, SEXP obj_minus, SEXP bed, SEXP op, SEXP step, SEXP use_strand, SEXP with_attributes, SEXP as_matrix, SEXP gap_value, SEXP abs_value, SEXP follow_strand) {
  bwStepOp bwOp;
  
  // initialize selected operation
  bw_select_op(&bwOp, CHAR(STRING_ELT(op, 0)), 1);
  
  return bigWig_region_query(obj_plus, obj_minus, bed, bwOp, step, use_strand, with_attributes, as_matrix, gap_value, abs_value, follow_strand, NULL, NULL);
}

SEXP bigWig_bp_chrom_query(SEXP obj, SEXP op, SEXP chrom, SEXP step, SEXP with_attributes, SEXP gap_value, SEXP abs_value, SEXP bwMap) {
  bigWig_t * bw = NULL;
  SEXP result;
  int istep;
  double d_gap;
  int do_abs = INTEGER(abs_value)[0] == TRUE;
  int no_step = 0;
  int use_attributes = INTEGER(with_attributes)[0] == TRUE;
  bwStepOp bwOp;
  const char * c_chrom;
  
  //
  // TODO: handle bwMap
  //

  // initialize selected operation
  bw_select_op(&bwOp, CHAR(STRING_ELT(op, 0)), 0);

  c_chrom = char_elt(chrom, 0);
  
  PROTECT(step = AS_INTEGER(step));
  istep = INTEGER(step)[0];
  no_step = (istep == NA_INTEGER);
  
  if (no_step || istep <= 0)
    error("step size must be >= 1");
  
  PROTECT(gap_value = AS_NUMERIC(gap_value));
  d_gap = REAL(gap_value)[0];


  bw = bigWig_for_chrom(obj, c_chrom);

  PROTECT(result = NEW_NUMERIC(bw_step_query_size(0, bw_chrom_size(bw, c_chrom), istep)));
  bw_chrom_step_query(bw, &bwOp, c_chrom, istep, d_gap, do_abs, REAL(result));
  
  // attributes
  if (use_attributes) {
    SEXP att_chrom = NEW_STRING(1);
    SEXP att_start = NEW_INTEGER(1);
    SEXP att_end = NEW_INTEGER(1);
    SEXP att_step = NEW_INTEGER(1);
        
    SET_STRING_ELT(att_chrom, 0, mkChar(c_chrom));
    INTEGER(att_start)[0] = 0;
    INTEGER(att_end)[0] = bw_chrom_size(bw, c_chrom);
    INTEGER(att_step)[0] = istep;
        
    setAttrib(result, install("chrom"), att_chrom);
    setAttrib(result, install("start"), att_start);
    setAttrib(result, install("end"), att_end);
    setAttrib(result, install("step"), att_step);
  }

  bigWig_for_chrom_release(obj, bw);
  
  UNPROTECT(3);
  
  return result;
}


SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue;
  SEXP names = getAttrib(list, R_NamesSymbol);
  int i;

  for (i = 0; i < length(list); i++)
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
	    elmt = VECTOR_ELT(list, i);
	    break;
	}
  
  return elmt;
}

struct bw_map_query_data {
  int readLen;
  int readLeftEdge;
  const char * op_name;
  double thresh;
  
  // only used in bw_with_map_step_query_func
  SEXP bwMap;
};

SEXP bw_map_step_query_func(bigWig_t * bigwig, bwStepOp * op, const char * chrom, int start, int end, int step, double gap_value, int do_abs, int is_plus, void * uptr) {
  struct bw_map_query_data * data = (struct bw_map_query_data *) uptr;
  
  double thresh = data->thresh * step; // here gap_value := threshold.fraction
  
  if ((!is_plus && data->readLeftEdge == 1) || (is_plus && data->readLeftEdge == 0)) {
    // need to shift coordinates / handle edges
    int size = (end - start)/step;
    int actual_end = start + size * step;
    
    // position i in 'read coordinates' corresponds to position i - readLen + 1 in 'mappability coordinates'
    start = start - data->readLen + 1;
    actual_end = actual_end - data->readLen + 1;
    
    if (actual_end <= 0) {
      // everything is outside scope, fill with NAs
      SEXP result;
      double * ptr;
      int i;

      PROTECT(result = NEW_NUMERIC(size));

      ptr = REAL(result);
      for (i = 0; i < size; ++i)
        ptr[i] = 1;  // if outside chrom, then unmappable
      
      UNPROTECT(1);
      
      return result;
    } else if (start < 0) { // but actual_end > 0
      SEXP result;
      int cur_start;
      int cur_end;
      int i;
      
      PROTECT(result = NEW_NUMERIC(size));
            
      for (cur_start = start, cur_end = cur_start + step, i = 0; cur_start < actual_end; cur_start += step, cur_end += step, ++i) {
        if (cur_start >= 0) {
          int j;
          int len;

          // everything else is a regular block
          SEXP tmp = R_bw_step_query(bigwig, op, chrom, cur_start, end, step, gap_value, do_abs, thresh);
          PROTECT(tmp);
          
          len = Rf_length(tmp);
          if (i + Rf_length(tmp) > size) // TODO: we're dropping the last element ... check if this is the right thing to do ...
            len = size - i;
          
          for (j = 0; j < len; ++j)
            REAL(result)[i + j] = REAL(tmp)[j];
          
          UNPROTECT(1);
          break;
        } else if (cur_end <= 0) {
          REAL(result)[i] = 1; // if outside chrom, then unmappable
        } else {
          SEXP tmp;
          // incomplete intervals: (valid ops: sum, avg, thresh)
          // if sum or thresh, same thing is valid
          // if avg, use sum then compute average
          if (!strcmp(data->op_name, "avg")) {
            bwStepOp opsum;
            
            bw_select_op(&opsum, "sum", 0);
            
            tmp = R_bw_step_query(bigwig, &opsum, chrom, 0, cur_end, cur_end, gap_value, do_abs, thresh);
            
            REAL(result)[i] = REAL(tmp)[0] / step;
          } else if (!strcmp(data->op_name, "thresh")) {
            bwStepOp opsum;
            
            bw_select_op(&opsum, "sum", 0);
            
            tmp = R_bw_step_query(bigwig, &opsum, chrom, 0, cur_end, cur_end, gap_value, do_abs, thresh);
            
            REAL(result)[i] = (REAL(tmp)[0] >= thresh ? 1 : 0);
          } else {
            tmp = R_bw_step_query(bigwig, op, chrom, 0, cur_end, cur_end, gap_value, do_abs, thresh);
            
            REAL(result)[i] = REAL(tmp)[0];
          }
        }
      }
      
      UNPROTECT(1);
      
      return result;
    } else {
      // regular query
      return R_bw_step_query(bigwig, op, chrom, start, actual_end, step, gap_value, do_abs, thresh);
    }
  } else
    return R_bw_step_query(bigwig, op, chrom, start, end, step, gap_value, do_abs, thresh);
}

void fill_bw_map_query_data(SEXP obj, const char * op_name, struct bw_map_query_data * data) {
  SEXP thresh = getListElement(obj, "threshold.fraction");
  SEXP readLeftEdge = getListElement(obj, "read.left.edge");
  SEXP readLen = getListElement(obj, "read.len");

  PROTECT(readLen = AS_INTEGER(readLen));
  PROTECT(readLeftEdge = AS_LOGICAL(readLeftEdge));
  PROTECT(thresh = AS_NUMERIC(thresh));
  
  data->readLen = INTEGER(readLen)[0];
  data->readLeftEdge = INTEGER(readLeftEdge)[0];
  data->op_name = op_name;
  data->thresh = REAL(thresh)[0];

  UNPROTECT(3);
}

SEXP bwMap_bp_query(SEXP obj, SEXP bed, SEXP op, SEXP step, SEXP with_attributes, SEXP as_matrix) {
  bwStepOp bwOp;
  const char * op_name = CHAR(STRING_ELT(op, 0));
  SEXP gap_value;
  SEXP abs_value;
  SEXP use_strand;
  SEXP follow_strand;
  SEXP bw_obj = getListElement(obj, "bw");
  
  SEXP result;
  struct bw_map_query_data data;
  
  // initialize selected operation
  bw_select_op(&bwOp, op_name, 0);
  
  PROTECT(abs_value = NEW_LOGICAL(1));
  INTEGER(abs_value)[0] = FALSE;
  PROTECT(use_strand = NEW_LOGICAL(1));
  INTEGER(use_strand)[0] = TRUE;
  PROTECT(gap_value = NEW_NUMERIC(1));
  REAL(gap_value)[0] = 0; // by default, if no info then regions are mappable
  PROTECT(follow_strand = NEW_LOGICAL(1));
  INTEGER(follow_strand)[0] = FALSE;
  
  fill_bw_map_query_data(obj, op_name, &data);
  
  //
  result = bigWig_region_query(bw_obj, bw_obj, bed, bwOp, step, use_strand, with_attributes, as_matrix, gap_value, abs_value, follow_strand, bw_map_step_query_func, &data);
  
  UNPROTECT(4);
  
  return result;
}

SEXP bw_with_map_step_query_func(bigWig_t * bigwig, bwStepOp * op, const char * chrom, int start, int end, int step, double gap_value, int do_abs, int is_plus, void * uptr) {
  SEXP res, res_map;
  double * rptr;
  double * rmapptr;
  int i;
  struct bw_map_query_data * data = (struct bw_map_query_data*) uptr;
  bigWig_t * bw_map;
  bwStepOp bwMapOp;
  
  // run regular query
  PROTECT(res = R_bw_step_query(bigwig, op, chrom, start, end, step, gap_value, do_abs, 0.0));
  
  // run bwMap query
  bw_map = bigWig_for_chrom(data->bwMap, chrom);
  bw_select_op(&bwMapOp, data->op_name, 0);
  
  PROTECT(res_map = bw_map_step_query_func(bigwig, &bwMapOp, chrom, start, end, step, 0, 1, is_plus, uptr));
  bigWig_for_chrom_release(data->bwMap, bw_map);
  
  // combine results, i.e., fill all unmappable positions (res_map = 1) with NA
  rptr = REAL(res);
  rmapptr = REAL(res_map);
  for (i = 0; i < Rf_length(res); ++i, ++rptr, ++rmapptr) {
    if (*rmapptr == 1)
      *rptr = NA_REAL;
  }
  
  //
  UNPROTECT(2);
  
  return res;
}

SEXP bigWig_bp_query(SEXP obj_plus, SEXP obj_minus, SEXP bed, SEXP op, SEXP step, SEXP use_strand, SEXP with_attributes, SEXP as_matrix, SEXP gap_value, SEXP abs_value, SEXP follow_strand, SEXP bwMap) {
  bwStepOp bwOp;

  // initialize selected operation
  bw_select_op(&bwOp, CHAR(STRING_ELT(op, 0)), 0);
  
  if (bwMap != R_NilValue) {
    // fill in bwMap information
    struct bw_map_query_data data;
    
    fill_bw_map_query_data(bwMap, "thresh", &data);
    
    data.bwMap = bwMap;
    
    return bigWig_region_query(obj_plus, obj_minus, bed, bwOp, step, use_strand, with_attributes, as_matrix, gap_value, abs_value, follow_strand, bw_with_map_step_query_func, &data);
  }
  
  return bigWig_region_query(obj_plus, obj_minus, bed, bwOp, step, use_strand, with_attributes, as_matrix, gap_value, abs_value, follow_strand, NULL, NULL);
}

