#
# Query functions for "probe" mode
#
original.expr <- function(arg, n) {
  newname <- substitute(arg)
 
   # Travel up the frame stack until we hit the top.
   for(i in 1:n) {
     oldname <- do.call("substitute", list(as.name(newname), parent.frame(i)))
     newname <- oldname
   }
   deparse(newname)
}

valid.bw <- function(bw) {
  if (!(class(bw) == "bigWig" || (is.character(bw) && length(bw) == 2 && (!any(is.na(bw) | is.null(bw)))))) {
    stop("invalid bigWig: ", original.expr(bw, 2))
  }
}

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

valid.probe.op <- function(op) {
  if (all(op != c("sum", "avg", "wavg", "min", "max")))
    stop("invalid probe operation: ", op)
}

valid.query.range <- function(start, end, index = NA, step = NA) {
  if (end <= start) {
    if (!is.na(index))
      stop("bed:", index, ": end must be > start: ", start, ", ", end)
    else
      stop("end must be > start: ", start, ", ", end)
  }
  if (start < 0 || end < 1) {
    if (!is.na(index))
      stop("bed:", index, ": start and end must be positive: ", start, ", ", end)
    else
      stop("start and end must be positive: ", start, ", ", end)
  }
  
  if (!is.na(step) && (end - start) %/% step == 0) {
    if (!is.na(index))
      stop("bed:", index, ": query region is less than a step wide")
    else
      stop("query region is less than a step wide")
  }
  
  if (!is.na(step) && (end - start) %% step > 0) {
    if (!is.na(index))
      warning("bed:", index, ": query region is not an integer multiple of step")
    else
      warning("query region is not an integer multiple of step")
  }
}

valid.strand <- function(strand) {
  is.character(strand) & (strand == "+" | strand == "-")
}

bed.valid.query.range <- function(bed, step = NA) {
  foreach.bed(bed, function(i, chrom, start, end, strand) {
    valid.query.range(start, end, index = i, step = step)
  })
}

region.probeQuery.bigWig <- function(bw, chrom, start, end, op = "wavg", abs.value = FALSE, gap.value = NA) {
  valid.bw(bw)
  valid.probe.op(op)
  valid.query.range(start, end)
  valid.chrom(bw, chrom)
  
  bed = data.frame(chrom, start, end)
  .Call(bigWig_probe_query, bw, NULL, bed, op, NA, FALSE, FALSE, TRUE, gap.value, abs.value, FALSE)
}

bed.region.probeQuery.bigWig <- function(bw, bed, op = "wavg", abs.value = FALSE, gap.value = NA) {
  valid.bw(bw)
  bed.valid.query.range(bed)
  valid.probe.op(op)
  .Call(bigWig_probe_query, bw, NULL, bed[, 1:3], op, NA, FALSE, FALSE, TRUE, gap.value, abs.value, FALSE)
}

bed6.region.probeQuery.bigWig <- function(bw.plus, bw.minus, bed6, op = "wavg", abs.value = FALSE, gap.value = NA) {
  valid.bw(bw.plus)
  valid.bw(bw.minus)
  stopifnot(dim(bed6)[2] >= 6)
  stopifnot(all(valid.strand(as.character(bed6[,6]))))
  bed.valid.query.range(bed6)
  valid.probe.op(op)
  .Call(bigWig_probe_query, bw.plus, bw.minus, bed6, op, NA, TRUE, FALSE, TRUE, gap.value, abs.value, FALSE)
}

# note: start, end are optional here (use NULL for both to get the entire choromosome)
step.probeQuery.bigWig <- function(bw, chrom, start, end, step, op = "wavg", abs.value = FALSE, gap.value = NA, with.attributes = TRUE) {
  valid.bw(bw)
  valid.probe.op(op)
  valid.chrom(bw, chrom)

  if ((is.null(start) && !is.null(end)) || (!is.null(start) && is.null(end)))
    stop("either set both start and end to null (chromosome-wide query) or neither")

  if (is.null(start) && is.null(end)) {
    # if we got a path, load the bigWig
    if (is.character(bw))
      bw = load.bigWig(paste(bw[1], chrom, bw[2], sep=''))
  
    chromIdx = which(bw$chroms == chrom)
    bed = data.frame(chrom, 0, bw$chromSizes[chromIdx])
    
    return(.Call(bigWig_probe_query, bw, NULL, bed, op, step, FALSE, with.attributes, FALSE, gap.value, abs.value, FALSE)[[1]])
  }
  
  valid.query.range(start, end, step = step)
  
  bed = data.frame(chrom, start, end)
  .Call(bigWig_probe_query, bw, NULL, bed, op, step, FALSE, with.attributes, FALSE, gap.value, abs.value, FALSE)[[1]]
}

bed.step.probeQuery.bigWig <- function(bw, bed, step, op = "wavg", abs.value = FALSE, gap.value = NA, with.attributes = FALSE, as.matrix = FALSE) {
  valid.bw(bw)
  bed.valid.query.range(bed, step = step)
  valid.probe.op(op)
  if (as.matrix) {
    sizes = bed[,3] - bed[,2]
    stopifnot(all(sizes == sizes[1]))
  }
  
  .Call(bigWig_probe_query, bw, NULL, bed[, 1:3], op, step, FALSE, with.attributes, as.matrix, gap.value, abs.value, FALSE)
}

bed6.step.probeQuery.bigWig <- function(bw.plus, bw.minus, bed6, step, op = "wavg", abs.value = FALSE, gap.value = NA, with.attributes = FALSE, as.matrix = FALSE, follow.strand = FALSE) {
  valid.bw(bw.plus)
  valid.bw(bw.minus)
  bed.valid.query.range(bed6, step = step)
  stopifnot(dim(bed6)[2] >= 6)
  stopifnot(all(valid.strand(as.character(bed6[,6]))))
  
  valid.probe.op(op)
  if (as.matrix) {
    sizes = bed6[,3] - bed6[,2]
    stopifnot(all(sizes == sizes[1]))
  }
  
  .Call(bigWig_probe_query, bw.plus, bw.minus, bed6, op, step, TRUE, with.attributes, as.matrix, gap.value, abs.value, follow.strand)
}
