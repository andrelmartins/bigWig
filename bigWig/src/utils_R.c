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

/*
 * 
 *  region stacking code 
 * 
 */

struct region {
  int index;

  int subtree_max_end;

  struct region * left;
  struct region * right;
};


static void insert(struct region ** root, struct region * rec, int key_new, int * starts, int this_end) {
  int key_this;
  struct region * ptr = *root;

  if (ptr == NULL) {
    *root = rec;
    return;
  }

  key_this = starts[ptr->index];
  if (this_end > ptr->subtree_max_end)
    ptr->subtree_max_end = this_end;

  if (key_this > key_new)
    insert(&(ptr->left), rec, key_new, starts, this_end);
  else
    insert(&(ptr->right), rec, key_new, starts, this_end);
}

struct end_frame {
  int max_end;
  struct end_frame * next;
};

static int place_record(struct end_frame ** frames, int start, int end) {
  int level = 0;
  struct end_frame * ptr;

  if (*frames == NULL) {
    level = 1;
    *frames = Calloc(1, struct end_frame);
    (*frames)->max_end = end;
    (*frames)->next = NULL;
    return level;
  }

  ptr = *frames;
  while(1) {
    ++level;

    if (ptr->max_end < start) {
      ptr->max_end = end;
      return level;
    }

    if (ptr->next == NULL) {
      ptr->next = Calloc(1, struct end_frame);
      ptr->next->max_end = end;
      ptr->next->next = NULL;
      return level + 1;
    }

    ptr = ptr->next;
  }
}

static void free_frames(struct end_frame * frames) {
  struct end_frame * next;

  while (frames != NULL) {
    next = frames->next;
    Free(frames);
    frames = next;
  }
}

static void recursive_stack_levels(struct end_frame ** frames, struct region * node, int * starts, int * ends, int * output) {

  if (node == NULL)
    return;

  /* visit left (smaller) */
  if (node->left != NULL)
    recursive_stack_levels(frames, node->left, starts, ends, output);

  /* visit node */
  output[node->index] = place_record(frames, starts[node->index], ends[node->index]);

  /* visit right (greater) */
  if (node->right != NULL)
    recursive_stack_levels(frames, node->right, starts, ends, output);
}

static void stack_levels(struct region * root, int * starts, int * ends, int * output) {
  struct end_frame * frames = NULL;

  recursive_stack_levels(&frames, root, starts, ends, output);
  free_frames(frames);
}

SEXP stack_bed(SEXP starts, SEXP ends) {
  SEXP result;
  int i, N;
  struct region * root = NULL;
  
  PROTECT(starts = AS_INTEGER(starts));
  PROTECT(ends = AS_INTEGER(ends));
  
  N = Rf_length(starts);
  PROTECT(result = NEW_INTEGER(N));
  
  /* add regions */
  for (i = 0; i < N; ++i) {
    struct region * rec = (struct region*) R_alloc(1, sizeof(struct region));

    rec->index = i;
    rec->left = NULL;
    rec->right = NULL;
    rec->subtree_max_end = INTEGER(ends)[i];

    insert(&root, rec, INTEGER(starts)[i], INTEGER(starts), INTEGER(ends)[i]);
  }
  
  /* run interval scheduling algorithm */
  stack_levels(root, INTEGER(starts), INTEGER(ends), INTEGER(result));
  
  UNPROTECT(3);
  
  return result;
}
