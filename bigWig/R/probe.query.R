#
# Query functions for "probe" mode
#

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

bed.step.probeQuery.bigWig <- function(bw, bed, step, op = "wavg", abs.value = FALSE, gap.value = NA, with.attributes = TRUE, as.matrix = FALSE) {
  valid.bw(bw)
  bed.valid.query.range(bed, step = step)
  valid.probe.op(op)
  if (as.matrix) {
    sizes = bed[,3] - bed[,2]
    stopifnot(all(sizes == sizes[1]))
  }
  
  .Call(bigWig_probe_query, bw, NULL, bed[, 1:3], op, step, FALSE, with.attributes, as.matrix, gap.value, abs.value, FALSE)
}

bed6.step.probeQuery.bigWig <- function(bw.plus, bw.minus, bed6, step, op = "wavg", abs.value = FALSE, gap.value = NA, with.attributes = TRUE, as.matrix = FALSE, follow.strand = FALSE) {
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
