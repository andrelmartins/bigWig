#include "bw_base.h"
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

int bw_has_chrom(bigWig_t * bw, const char * chromName) {
  struct bbiChromInfo * chrom, * chromList = bbiChromList(bw);
    
  for (chrom = chromList; chrom != NULL; chrom = chrom->next)
    if (!strcmp(chromName, chrom->name)) {
      bbiChromInfoFreeList(&chromList);
      return 1;
    }

  bbiChromInfoFreeList(&chromList);

  return 0;
}

long bw_chrom_size(bigWig_t * bw, const char * chromName) {
  struct bbiChromInfo * chrom, * chromList = bbiChromList(bw);
 
  for (chrom = chromList; chrom != NULL; chrom = chrom->next)
    if (!strcmp(chromName, chrom->name)) {
      long result = chrom->size;
      bbiChromInfoFreeList(&chromList);
      return result;
    }

  return -1;
}
