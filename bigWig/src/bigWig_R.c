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

#include <R.h>
#include <Rdefines.h>

SEXP bigWig_load(SEXP filename) {
  return R_NilValue;
}
