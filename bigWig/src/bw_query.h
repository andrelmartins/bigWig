#ifndef BW_QUERY_H
#define BW_QUERY_H

#include "bw_base.h"

typedef struct {
  double defaultValue;
  int do_abs;
  
  double total;
  double count;
  double thresh;
} bwStepOpData;

typedef void (* bw_op_clear)(bwStepOpData * data);
typedef void (* bw_op_add)(bwStepOpData * data, double isize, double ivalue);
typedef double (* bw_op_result)(bwStepOpData * data, int step);

typedef struct {
  bw_op_clear clear;
  bw_op_add add;
  bw_op_result result;
} bwStepOp;

void bw_select_op(bwStepOp * op, const char * bw_op_type, int probe_mode);
int bw_step_query_size(int start, int end, int step);
int bw_step_query(bigWig_t * bigwig, bwStepOp * op, const char * chrom, int start, int end, int step, double gap_value, int do_abs, double thresh, double * buffer);
int bw_chrom_step_query(bigWig_t * bigwig, bwStepOp * op, const char * chrom, int step, double gap_value, int do_abs, double * buffer);

#endif
