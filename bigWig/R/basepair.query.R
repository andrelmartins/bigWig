#
# Query functions for "base pair" mode
#
valid.chrom <- function(bw, chrom) {
  stopifnot(is.character(chrom))

  if (class(bw) == "bigWig") {
    if (!(chrom %in% bw$chroms))
     stop("bigWig does not contain information on chromosome: ", chrom)
    return(TRUE)
  } else if (is.character(bw) && length(bw) == 2) {
    path = paste(bw[1], chrom, bw[2], sep='')
    if (file.exists(path))
      return(TRUE)
    stop("can't find bigWig file for chromosome: ", chrom)
  } else
    stop("invalid bigWig object")
}

valid.bp.op <- function(op) {
  if (all(op != c("sum", "avg", "min", "max")))
    stop("invalid base pair operation: ", op)
}

valid.strand <- function(strand) {
  !is.na(strand) & (strand == "+" | strand == "-")
}

region.bpQuery.bigWig <- function(bw, chrom, start, end, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL) {
  if (!is.null(bwMap) && !valid.strand(strand))
    stop("strand is required when using mappability information")

  valid.bp.op(op)
  valid.query.range(start, end)
  valid.chrom(bw, chrom)
  
  if (!is.na(strand)) {
    bed = data.frame(chrom, start, end, 0, 0, strand)
    .Call(bigWig_bp_query, bw, bw, bed, op, NA, TRUE, FALSE, TRUE, gap.value, abs.value, bwMap)
  } else {
    bed = data.frame(chrom, start, end)
    .Call(bigWig_bp_query, bw, NULL, bed, op, NA, FALSE, FALSE, TRUE, gap.value, abs.value, bwMap)
  }
}

bed.region.bpQuery.bigWig <- function(bw, bed, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL) {
  if (!is.null(bwMap) && !valid.strand(strand))
    stop("strand is required when using mappability information")
  
  valid.bp.op(op)
  bed.valid.query.range(bed)
    
  if (!is.na(strand)) {
    bed = cbind(bed, data.frame(0, 0, strand))
    .Call(bigWig_bp_query, bw, bw, bed, op, NA, TRUE, FALSE, TRUE, gap.value, abs.value, bwMap)
  } else {
    .Call(bigWig_bp_query, bw, NULL, bed, op, NA, FALSE, FALSE, TRUE, gap.value, abs.value, bwMap)
  }
}

bed6.region.bpQuery.bigWig <- function(bw.plus, bw.minus, bed6, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL) {
  stopifnot(dim(bed6) >= 6)
  stopifnot(all(valid.strand(as.character(bed6[,6]))))
  valid.bp.op(op)
  bed.valid.query.range(bed6)
  
  .Call(bigWig_bp_query, bw.plus, bw.minus, bed6, op, NA, TRUE, FALSE, TRUE, gap.value, abs.value, bwMap)
}

# note: start, end are optional here (use NULL for both to get the entire choromosome)
step.bpQuery.bigWig <- function(bw, chrom, start, end, step, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL, with.attributes = TRUE) {
  if (!is.null(bwMap) && !valid.strand(strand))
    stop("strand is required when using mappability information")
  
  valid.bp.op(op)
  valid.chrom(bw, chrom)

  if ((is.null(start) && !is.null(end)) || (!is.null(start) && is.null(end)))
    stop("either set both start and end to null (chromosome-wide query) or neither")
  
  if (is.null(start) && is.null(end)) {
    result = .Call(bigWig_bp_chrom_query, bw, op, chrom, step, with.attributes, gap.value, abs.value, bwMap)
    
    if (!is.na(strand) && strand == '-') {
      ats = attributes(result)
      result = rev(result)
      attributes(result) <- ats
    }

    return(result)
  }
    
  valid.query.range(start, end, step = step)
      
  if (!is.na(strand)) {
    bed = data.frame(chrom, start, end, 0, 0, strand)
    .Call(bigWig_bp_query, bw, NULL, bed, op, step, TRUE, with.attributes, FALSE, gap.value, abs.value, bwMap)[[1]]
  } else {
    bed = data.frame(chrom, start, end)
    .Call(bigWig_bp_query, bw, NULL, bed, op, step, FALSE, with.attributes, FALSE, gap.value, abs.value, bwMap)[[1]]
  }
}

bed.step.bpQuery.bigWig <- function(bw, bed, step, strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL, with.attributes = FALSE, as.matrix = FALSE) {
  if (!is.null(bwMap) && !valid.strand(strand))
    stop("strand is required when using mappability information")
  
  valid.bp.op(op)
  bed.valid.query.range(bed, step = step)
  
  if (!is.na(strand)) {
    bed = cbind(bed, data.frame(0, 0, strand))
    .Call(bigWig_bp_query, bw, bw, bed, op, step, TRUE, with.attributes, as.matrix, gap.value, abs.value, bwMap)
  } else {
    .Call(bigWig_bp_query, bw, NULL, bed, op, step, FALSE, with.attributes, as.matrix, gap.value, abs.value, bwMap)
  }
}

bed6.step.bpQuery.bigWig <- function(bw.plus, bw.minus, bed6, step, op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL, with.attributes = FALSE, as.matrix = FALSE) {
  stopifnot(dim(bed6)[2] >= 6)
  stopifnot(all(valid.strand(as.character(bed6[,6]))))
  valid.bp.op(op)
  bed.valid.query.range(bed6, step = step)

  if (as.matrix) {
    sizes = bed6[,3] - bed6[,2]
    stopifnot(all(sizes == sizes[1]))
  }
  
  .Call(bigWig_bp_query, bw.plus, bw.minus, bed6, op, step, TRUE, with.attributes, as.matrix, gap.value, abs.value, bwMap)
}
