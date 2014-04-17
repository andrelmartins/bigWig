#
# functions for mappability information
#

load.bwMap <- function(filename, read.len, read.left.edge, threshold.fraction = 0, udcDir=NULL) {
  read.len = as.integer(read.len)
  stopifnot(read.len > 0)
  stopifnot(is.logical(read.left.edge))
  stopifnot(threshold.fraction >= 0 && threshold.fraction <= 1)
  
  bw = load.bigWig(filename, udcDir = udcDir)
  
  if (bw$max != 1 || (bw$min != 1 && bw$min != 0)) {
    unload.bigWig(bw)
    stop("mappability bigWig file should only have zero or one values (one for unmappable positions)")
  }
  
  res = list(bw = bw, read.len = read.len, read.left.edge = read.left.edge, threshold.fraction = threshold.fraction)
  
  class(res) <- "bwMap"
  return(res)
}

unload.bwMap <- function(bwMap) {
  stopifnot(class(bwMap) == "bwMap")
  unload.bigWig(bwMap$bw)
}

print.bwMap <- function(x, ...) {
  cat.bool <- function(name, value) {
    cat(" ", name, ": ", sep='')
    if (value)
      cat("yes\n")
    else
      cat("no\n")
  }

  cat("bwMap\n")
  cat(" chromCount:", length(x$bw$chroms), "\n")
  for (i in 1:length(x$bw$chroms))
    cat("    ", x$bw$chroms[i], x$bw$chromSizes[i], "\n")
  cat(" basesCovered:", prettyNum(x$bw$basesCovered, big.mark=','), "\n")
  cat(" readLength:", x$read.len, "\n")
  cat.bool("readLeftEdge", x$read.left.edge)
  cat(" stepFractionThreshold:", x$threshold.fraction, "\n")
}

#
# Query functions for "base pair" mode
#

valid.bwMap.op <- function(op) {
  if (all(op != c("sum", "avg", "thresh")))
    stop("invalid base pair operation: ", op)
}

region.bpQuery.bwMap <- function(bwMap, chrom, start, end, strand, op = "thresh") {
  if (!valid.strand(strand))
    stop("strand is required when using mappability information")

  valid.bwMap.op(op)
  valid.query.range(start, end)
  
  if (!any(bwMap$bw$chroms == chrom)) {
    stop("bigWig does not contain information on chromosome: ", chrom)
  }

  bed = data.frame(chrom, start, end, 0, 0, strand)
  .Call(bwMap_bp_query, bwMap, bed, op, NA, FALSE, TRUE)
}

bed6.region.bpQuery.bwMap <- function(bwMap, bed6, op = "thresh") {
  stopifnot(dim(bed6) >= 6)
  stopifnot(all(valid.strand(as.character(bed6[,6]))))
  valid.bwMap.op(op)
  bed.valid.query.range(bed6)
  
  .Call(bwMap_bp_query, bwMap, bed6, op, NA, FALSE, TRUE)
}

# note: start, end are optional here (use NULL for both to get the entire choromosome)
step.bpQuery.bwMap <- function(bwMap, chrom, start, end, step, strand, op = "thresh", with.attributes = TRUE) {
  if (!valid.strand(strand))
    stop("strand is required when using mappability information")
  
  valid.bwMap.op(op)
  
  if (!any(bwMap$bw$chroms == chrom)) {
    stop("bigWig does not contain information on chromosome: ", chrom)
  }
  
  if (is.null(start) && is.null(end)) {
    chromIdx = which(bwMap$bw$chroms == chrom)
    bed = data.frame(chrom, 0, bwMap$bw$chromSizes[chromIdx], 0, 0, strand)
    
    return(.Call(bwMap_bp_query, bwMap, bed, op, step, with.attributes, FALSE)[[1]])
  }
  if (is.null(start) || is.null(end))
    stop("either set both start and end to null (chromosome-wide query) or neither")

  valid.query.range(start, end, step = step)
    
  bed = data.frame(chrom, start, end, 0, 0, strand)
  .Call(bwMap_bp_query, bwMap, bed, op, step, with.attributes, FALSE)[[1]]
}

bed6.step.bpQuery.bwMap <- function(bwMap, bed6, step, op = "thresh", with.attributes = FALSE, as.matrix = FALSE) {
  stopifnot(dim(bed6)[2] >= 6)
  stopifnot(all(valid.strand(as.character(bed6[,6]))))
  valid.bwMap.op(op)
  bed.valid.query.range(bed6, step = step)

  if (as.matrix) {
    sizes = bed6[,3] - bed6[,2]
    stopifnot(all(sizes == sizes[1]))
  }
  
  .Call(bwMap_bp_query, bwMap, bed6, op, step, with.attributes, as.matrix)
}
