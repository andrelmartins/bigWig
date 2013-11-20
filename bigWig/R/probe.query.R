#
# Query functions for "probe" mode
#

valid.probe.op <- function(op) {
  if (all(op != c("sum", "avg", "wavg", "min", "max")))
    stop("invalid probe operation: ", op)
}

region.probeQuery.bigWig <- function(bw, chrom, start, end, op = "sum", abs.value = FALSE, gap.value = NA) {
  valid.probe.op(op)
  if (!any(bw$chroms == chrom)) {
    warning("bigWig does not contain information on chromosome: ", chrom)
    return(gap.value)
  }
  bed = data.frame(chrom, start, end)
  as.vector(.Call(bigWig_probe_query, bw, bed, op, NA, FALSE, FALSE, TRUE, gap.value, abs.value))
}

bed.region.probeQuery.bigWig <- function(bw, bed, op = "sum", abs.value = FALSE, gap.value = NA) {
  valid.probe.op(op)
  as.vector(.Call(bigWig_probe_query, bw, bed, op, NA, FALSE, FALSE, TRUE, gap.value, abs.value))
}

# note: start, end are optional here (use NULL for both to get the entire choromosome)
step.probeQuery.bigWig <- function(bw, chrom, start, end, step, op = "sum", abs.value = FALSE, gap.value = NA, with.attributes = TRUE) {
  valid.probe.op(op)
  if (!any(bigWig$chroms == chrom)) {
    warning("bigWig does not contain information on chromosome: ", chrom)
    return(rep(gap.value, (end - start) %/% step))
  }
  if (is.null(start) && is.null(end)) {
  }
  if (is.null(start) || is.null(end))
    stop("either set both start and end to null (chromosome wide query) or neither")
  
  bed = data.frame(chrom, start, end)
  .Call(bigWig_probe_query, bw, bed, op, step, FALSE, with.attributes, FALSE, gap.value, abs.value)
}

bed.step.probeQuery.bigWig <- function(bw, bed, step, op = "sum", abs.value = FALSE, gap.value = NA, with.attributes = FALSE, as.matrix = FALSE) {
  valid.probe.op(op)
  if (as.matrix) {
    sizes = bed[,3] - bed[,2]
    stopifnot(all(sizes) == sizes[1])
  }
  
  .Call(bigWig_probe_query, bw, bed, op, step, FALSE, with.attributes, as.matrix, gap.value, abs.value)
}

bed6.step.probeQuery.bigWig <- function(bw.plus, bw.minus, bed6, step, op = "sum", abs.value = FALSE, gap.value = NA, with.attributes = FALSE, as.matrix = FALSE) {
  valid.probe.op(op)
  if (as.matrix) {
    sizes = bed[,3] - bed[,2]
    stopifnot(all(sizes) == sizes[1])
  }
  
  .Call(bigWig_probe_query, bw, bed, op, step, TRUE, with.attributes, as.matrix, gap.value, abs.value)
}
