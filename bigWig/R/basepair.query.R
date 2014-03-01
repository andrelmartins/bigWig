#
# Query functions for "base pair" mode
#

region.bpQuery.bigWig <- function(bw, chrom, start, end, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL) {
  if (!is.null(bwMap) && !valid.strand(strand))
    stop("strand is required when using mappability information")

  valid.bw(bw)
  valid.bp.op(op)
  valid.query.range(start, end)
  valid.chrom(bw, chrom)
  
  if (!is.na(strand)) {
    stopifnot(valid.strand(strand))
    bed = data.frame(chrom, start, end, 0, 0, strand)
    .Call(bigWig_bp_query, bw, bw, bed, op, NA, TRUE, FALSE, TRUE, gap.value, abs.value, FALSE, bwMap)
  } else {
    bed = data.frame(chrom, start, end)
    .Call(bigWig_bp_query, bw, NULL, bed, op, NA, FALSE, FALSE, TRUE, gap.value, abs.value, FALSE, bwMap)
  }
}

bed.region.bpQuery.bigWig <- function(bw, bed, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL) {
  if (!is.null(bwMap) && !valid.strand(strand))
    stop("strand is required when using mappability information")

  valid.bw(bw)
  valid.bp.op(op)
  bed.valid.query.range(bed)
    
  if (!is.na(strand)) {
    stopifnot(valid.strand(strand))
    bed = cbind(bed, data.frame(0, 0, strand))
    .Call(bigWig_bp_query, bw, bw, bed, op, NA, TRUE, FALSE, TRUE, gap.value, abs.value, bwMap)
  } else {
    .Call(bigWig_bp_query, bw, NULL, bed, op, NA, FALSE, FALSE, TRUE, gap.value, abs.value, FALSE, bwMap)
  }
}

bed6.region.bpQuery.bigWig <- function(bw.plus, bw.minus, bed6, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL) {
  valid.bw(bw.plus)
  valid.bw(bw.minus)
  stopifnot(dim(bed6)[2] >= 6)
  stopifnot(all(valid.strand(as.character(bed6[,6]))))
  valid.bp.op(op)
  bed.valid.query.range(bed6)
  
  .Call(bigWig_bp_query, bw.plus, bw.minus, bed6, op, NA, TRUE, FALSE, TRUE, gap.value, abs.value, FALSE, bwMap)
}

# note: start, end are optional here (use NULL for both to get the entire choromosome)
step.bpQuery.bigWig <- function(bw, chrom, start, end, step, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL, with.attributes = TRUE) {
  if (!is.null(bwMap) && !valid.strand(strand))
    stop("strand is required when using mappability information")
  
  valid.bw(bw)
  valid.bp.op(op)
  valid.chrom(bw, chrom)

  if ((is.null(start) && !is.null(end)) || (!is.null(start) && is.null(end)))
    stop("either set both start and end to null (chromosome-wide query) or neither")
  
  if (is.null(start) && is.null(end)) {
    stopifnot(is.na(strand) || valid.strand(strand))
    result = .Call(bigWig_bp_chrom_query, bw, op, chrom, step, with.attributes, gap.value, abs.value, FALSE, bwMap)
    
    # for now handle bwMap here
    if (!is.null(bwMap)) {
      res_map = step.bpQuery.bwMap(bwMap, chrom, NULL, NULL, step, strand, op = "thresh", with.attributes = FALSE)
      
      result[res_map == 1] = NA
      res_map = NULL
    }

    return(result)
  }
    
  valid.query.range(start, end, step = step)
      
  if (!is.na(strand)) {
    stopifnot(valid.strand(strand))
    bed = data.frame(chrom, start, end, 0, 0, strand)
    .Call(bigWig_bp_query, bw, NULL, bed, op, step, TRUE, with.attributes, FALSE, gap.value, abs.value, FALSE, bwMap)[[1]]
  } else {
    bed = data.frame(chrom, start, end)
    .Call(bigWig_bp_query, bw, NULL, bed, op, step, FALSE, with.attributes, FALSE, gap.value, abs.value, FALSE, bwMap)[[1]]
  }
}

bed.step.bpQuery.bigWig <- function(bw, bed, step, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL, with.attributes = TRUE, as.matrix = FALSE) {
  if (!is.null(bwMap) && !valid.strand(strand))
    stop("strand is required when using mappability information")
    
  valid.bw(bw)
  valid.bp.op(op)
  bed.valid.query.range(bed, step = step)
  
  if (!is.na(strand)) {
    stopifnot(valid.strand(strand))
    bed = cbind(bed, data.frame(0, 0, strand))
    .Call(bigWig_bp_query, bw, bw, bed, op, step, TRUE, with.attributes, as.matrix, gap.value, abs.value, FALSE, bwMap)
  } else {
    .Call(bigWig_bp_query, bw, NULL, bed, op, step, FALSE, with.attributes, as.matrix, gap.value, abs.value, FALSE, bwMap)
  }
}

bed6.step.bpQuery.bigWig <- function(bw.plus, bw.minus, bed6, step, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL, with.attributes = TRUE, as.matrix = FALSE, follow.strand = FALSE) {
  valid.bw(bw.plus)
  valid.bw(bw.minus)
  stopifnot(dim(bed6)[2] >= 6)
  stopifnot(all(valid.strand(as.character(bed6[,6]))))
  valid.bp.op(op)
  bed.valid.query.range(bed6, step = step)

  if (as.matrix) {
    sizes = bed6[,3] - bed6[,2]
    stopifnot(all(sizes == sizes[1]))
  }
  
  .Call(bigWig_bp_query, bw.plus, bw.minus, bed6, op, step, TRUE, with.attributes, as.matrix, gap.value, abs.value, follow.strand, bwMap)
}
