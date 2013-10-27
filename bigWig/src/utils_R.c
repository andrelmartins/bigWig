/*
  C implementation of utility functions
*/

#include <R.h>
#include <Rdefines.h>

SEXP foreach_bed(SEXP bed, SEXP function, SEXP envir) {
  SEXP chroms, starts, ends, strands = R_NilValue;
  SEXP arg_idx, arg_chrom, arg_start, arg_end, arg_strand;
  SEXP fcall;
  int has_strand = 0;
  int i, N;

  PROTECT(bed);

  if(!isEnvironment(envir)) 
    error("'envir' should be an environment");

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

  /* build function call */
  PROTECT(arg_idx = NEW_INTEGER(1));
  PROTECT(arg_chrom = allocVector(STRSXP, 1));
  PROTECT(arg_start = NEW_INTEGER(1));
  PROTECT(arg_end = NEW_INTEGER(1));
  PROTECT(arg_strand = allocVector(STRSXP, 1));

  if (has_strand == 0)
    SET_STRING_ELT(arg_strand, 0, NA_STRING);

  PROTECT(fcall = lang6(function, arg_idx, arg_chrom, arg_start, arg_end, arg_strand));

  /* run loop */
  N = length(chroms);
  for (i = 0; i < N; ++i) {
    INTEGER(arg_idx)[0] = i + 1;
    INTEGER(arg_start)[0] = INTEGER(starts)[i];
    INTEGER(arg_end)[0] = INTEGER(ends)[i];
    
    /* chrom */
    if (isFactor(chroms)) {
      int idx = INTEGER(chroms)[i] - 1;
      SET_STRING_ELT(arg_chrom, 0, STRING_ELT(GET_LEVELS(chroms), idx));
    } else
      SET_STRING_ELT(arg_chrom, 0, STRING_ELT(chroms, i));
    
    /* strand */
    if (has_strand == 1) {
      if (isFactor(strands)) {
        int idx = INTEGER(strands)[i] - 1;
        SET_STRING_ELT(arg_strand, 0, STRING_ELT(GET_LEVELS(strands), idx));
      } else
        SET_STRING_ELT(arg_strand, 0, STRING_ELT(strands, i));
    }
    
    /* eval and ignore result */
    eval(fcall, envir);
    
    /* check for interrupts */
    R_CheckUserInterrupt();
  }
  
  /* clean up */
  UNPROTECT(9);
  
  return R_NilValue;
}
