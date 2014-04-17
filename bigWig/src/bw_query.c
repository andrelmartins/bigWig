#include "bw_query.h"
#include <string.h>
 #include <math.h>
#include "common.h"
#include "bigWig.h"
#include "localmem.h"

void bw_op_clearz(bwStepOpData * data) {
  data->total = 0;
  data->count = 0;
}
void bw_op_clear_min(bwStepOpData * data) {
  data->total = INFINITY; //R_PosInf;
  data->count = 0;
}
void bw_op_clear_max(bwStepOpData * data) {
  data->total = -INFINITY; //R_NegInf;
  data->count = 0;
}

void bw_op_add_sum_probe(bwStepOpData * data, double isize, double ivalue) {
  if (data->do_abs)
    ivalue = fabs(ivalue);
  data->total += ivalue;
  data->count = 1.0;
}
void bw_op_add_sum_bp(bwStepOpData * data, double isize, double ivalue) {
  if (data->do_abs)
    ivalue = fabs(ivalue);
  data->total += ivalue * isize;
  data->count = 1.0;
}
void bw_op_add_avg_probe(bwStepOpData * data, double isize, double ivalue) {
  if (data->do_abs)
    ivalue = fabs(ivalue);
  data->total += ivalue;
  data->count += 1.0;
}
void bw_op_add_avg_bp(bwStepOpData * data, double isize, double ivalue) {
  if (data->do_abs)
    ivalue = fabs(ivalue);
  data->total += ivalue * isize;
  data->count = 1.0;
}
void bw_op_add_wavg_probe(bwStepOpData * data, double isize, double ivalue) {
  if (data->do_abs)
    ivalue = fabs(ivalue);
  data->total += ivalue * isize;
  data->count += isize;
}
void bw_op_add_min(bwStepOpData * data, double isize, double ivalue) {
  if (data->do_abs)
    ivalue = fabs(ivalue);
  if (data->total > ivalue)
    data->total = ivalue;
  data->count = 1.0;
}
void bw_op_add_max(bwStepOpData * data, double isize, double ivalue) {
  if (data->do_abs)
    ivalue = fabs(ivalue);
  if (data->total < ivalue)
    data->total = ivalue;
  data->count = 1.0;
}
void bw_op_add_thresh(bwStepOpData * data, double isize, double ivalue) {
  data->total += ivalue * isize;
  data->count = 1.0;
}

double bw_op_result_sum_min_max(bwStepOpData * data, int step) {
  if (data->count == 0.0)
    return data->defaultValue;
  return data->total;
}
double bw_op_result_avg_wavg(bwStepOpData * data, int step) {
  if (data->count == 0.0)
    return data->defaultValue;
  return data->total / data->count;
}
double bw_op_result_avg_bp(bwStepOpData * data, int step) {
  if (data->count == 0.0)
    return data->defaultValue;
  return data->total / step;
}
double bw_op_result_thresh(bwStepOpData * data, int step) {
  if (data->count == 0.0)
    return data->defaultValue;
  
  if (data->total >= data->thresh)
    return 1.0;
  else
    return 0.0;
}

void bw_select_op(bwStepOp * op, const char * bw_op_type, int probe_mode) {
  if (!strcmp("sum", bw_op_type)) {
    op->clear = bw_op_clearz;
    if (probe_mode)
      op->add = bw_op_add_sum_probe;
    else
      op->add = bw_op_add_sum_bp;
    op->result = bw_op_result_sum_min_max;
  } else if (!strcmp("avg", bw_op_type)) {
    op->clear = bw_op_clearz;
    if (probe_mode) {
      op->add = bw_op_add_avg_probe;
      op->result = bw_op_result_avg_wavg;
    } else {
      op->add = bw_op_add_avg_bp;
      op->result = bw_op_result_avg_bp;
    }
  } else if (!strcmp("wavg", bw_op_type)) {
    op->clear = bw_op_clearz;
    if (probe_mode)
      op->add = bw_op_add_wavg_probe;
    else {
      /* error */
    }
    op->result = bw_op_result_avg_wavg;
  } else if (!strcmp("min", bw_op_type)) {
    op->clear = bw_op_clear_min;
    op->add = bw_op_add_min;
    op->result = bw_op_result_sum_min_max;
  } else if (!strcmp("max", bw_op_type)) {
    op->clear = bw_op_clear_max;
    op->add = bw_op_add_max;
    op->result = bw_op_result_sum_min_max;
  } else if (!strcmp("thresh", bw_op_type)) { /* this op uses 'defaultValue' (aka gap_value) as threshold */
    op->clear = bw_op_clearz;
    op->add = bw_op_add_thresh;
    op->result = bw_op_result_thresh;
  } else {
      /* throw error! */
  }
}

double interval_size(double step_start, double step_end, double istart, double iend) {
  double start = (step_start > istart ? step_start : istart);
  double end = (step_end > iend ? iend : step_end);
  return end - start;
}

int bw_step_query_size(int start, int end, int step) {
  return (end - start)/step;
}

/* TODO:
  - mappability: pass vector of step sum of mappable values and use that as the step value
                 (or fill steps with NA if "step value" == 0)

  - profile / optimize : 
    - maybe rewrite this file in C++ and use a function template / functors for
      the different clear/add/result "traits"...
    - maybe move clipping to only the code that needs it
    - same with mapping ...
    - optimized versions for chrom_step_query ... (all "intervals" are of size 1)
 */

int bw_step_query(bigWig_t * bigwig, bwStepOp * op, const char * chrom, int start, int end, int step, double gap_value, int do_abs, double thresh, double * buffer) {
  bwStepOpData data;
  struct lm * localMem = lmInit(0); /* use default value */
  struct bbiInterval * intervals;
  int nIntervals;
  int size = bw_step_query_size(start, end, step);
  int idx;     

  data.defaultValue = gap_value;
  data.do_abs = do_abs;
  data.thresh = thresh;

  /* collect bigWig intervals in range */
  intervals = bigWigIntervalQuery(bigwig, (char*) chrom, start, end, localMem); /* maybe always use the expanded interval ?? to simplify the code, or have an option in the op selection to picks one or the other ... */
  nIntervals = slCount(intervals);
  
  /* initialize result vector */
  for (idx = 0; idx < size; ++idx)
    buffer[idx] = gap_value;
  
  /* fill in values */
  if (nIntervals > 0) {
    struct bbiInterval * interval;
    double left, right;
    double i_size;
   
    /* initialize */
    (*(op->clear))(&data);
    left = start;
    right = start + step;
    idx = 0;
    for (interval = intervals; interval != NULL && idx < size; interval = interval->next) {
      /* interval starts beyond current step */
      if (((double)interval->start) >= right) {
        /* save current value */
        buffer[idx] = (*(op->result))(&data, step);
        
        /* start next */
        (*(op->clear))(&data);          

        while (idx < size && ((double)interval->start) >= right) {
          ++idx;
          left += step;
          right += step;
        }
      }

      /* interval starts at or before the current step */
      if (((double)interval->start) < right) {
        i_size = interval_size(left, right, ((double)interval->start), ((double)interval->end));
        (*(op->add))(&data, i_size, interval->val);
      }
      
      /* interval ends beyond or at the current step */
      if (((double)interval->end) >= right && idx < size) {
        do {
          /* save current value */
          buffer[idx] = (*(op->result))(&data, step);
        
          /* start next */
          (*(op->clear))(&data);
          
          ++idx;
          left += step;
          right += step;

          i_size = interval_size(left, right, ((double)interval->start), ((double)interval->end));
          (*(op->add))(&data, i_size, interval->val);
        } while (((double)interval->end) >= right && idx < size);
        if (left >= ((double)interval->end)) {
          /* current interval now ends before left */
          (*(op->clear))(&data);
        }
      }
    }
    
    if (idx < size)
      buffer[idx] = (*(op->result))(&data, step);
  }
  
  lmCleanup(&localMem);
  
  return nIntervals;
}

/* not valid for probe mode! */
int bw_chrom_step_query(bigWig_t * bigwig, bwStepOp * op, const char * chrom, int step, double gap_value, int do_abs, double * buffer) {
  int n;
  int i, j, k;
  double * valptr;
  
  struct bigWigValsOnChrom *chromVals = bigWigValsOnChromNew();
  bwStepOpData data;
  
  data.defaultValue = gap_value;
  data.do_abs = do_abs;
  
  if (!bigWigValsOnChromFetchData(chromVals, (char*) chrom, bigwig))
    return -1;
  
  n = bw_step_query_size(0, chromVals->chromSize, step);
  
  for (i = 0, j = 0, valptr=chromVals->valBuf; i < n; ++i) {
    (*(op->clear))(&data);
    
    for (k = 0; k < step; ++k, ++j, ++valptr) {
      if (bitReadOne(chromVals->covBuf , j))
        (*(op->add))(&data, 1, *valptr);
    }

    buffer[i] = (*(op->result))(&data, step);
  }
  
  bigWigValsOnChromFree(&chromVals);

  return 0;
}


