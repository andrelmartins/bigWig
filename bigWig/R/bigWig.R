#
# R interface
#

load.bigWig <- function(filename, udcDir=NULL) {
  res = .Call(bigWig_load, filename, udcDir)
  class(res) <- "bigWig"
  return(res)
}

unload.bigWig <- function(bigWig) {
  invisible(.Call(bigWig_unload, bigWig))
}

query.bigWig <- function(bigWig, chrom, start, end, clip = TRUE) {
  if (!any(bigWig$chroms == chrom))
    warning("bigWig does not contain information on chromosome: ", chrom)
  res <- .Call(bigWig_query, bigWig, chrom, start, end, clip)
  if (!is.null(res))
    colnames(res) <- c("start", "end", "value")
  return(res)
}

print.bigWig <- function(x, ...) {
  cat.bool <- function(name, value) {
    cat(" ", name, ": ", sep='')
    if (value)
      cat("yes\n")
    else
      cat("no\n")
  }
  cat("bigWig\n")
  cat(" version:", x$version, "\n")
  cat.bool("isCompressed", x$isCompressed)
  cat.bool("isSwapped", x$isSwapped)
  cat(" primaryDataSize:", prettyNum(x$primaryDataSize, big.mark=','), "\n")
  cat(" primaryIndexSize:", prettyNum(x$primaryIndexSize, big.mark=','), "\n")
  cat(" zoomLevels:", x$zoomLevels, "\n")
  cat(" chromcount:", length(x$chroms), "\n")
  for (i in 1:length(x$chroms))
    cat("    ", x$chroms[i], x$chromSizes[i], "\n")
  cat(" basesCovered:", prettyNum(x$basesCovered, big.mark=','), "\n")
  cat(" mean: ", x$mean, "\n")
  cat(" min: ", x$min, "\n")
  cat(" max: ", x$max, "\n")
  cat(" std: ", x$std, "\n")
}

bedQuery.bigWig <- function(bed, bwPlus, bwMinus = NULL, gapValue = NULL, weighted = F, aggregator = "sum") {
  return(.Call(bigWig_bed_query, bed, bwPlus, bwMinus, gapValue, weighted, aggregator))
}

queryByStep.bigWig <- function(bigWig, chrom, start, end, step, do.sum=F, default.null = T, defaultValue = 0) {
  stopifnot(step >= 1)
  if (!any(bigWig$chroms == chrom))
    warning("bigWig does not contain information on chromosome: ", chrom)

  result = .Call(bigWig_query_by_step, bigWig, chrom, start, end, step, do.sum, defaultValue)
  if (is.null(result) && default.null == FALSE)
    result = rep(defaultValue, (end - start) %/% step)
  
  return(result)
}

chromStepSum.bigWig <- function(bigwig, chrom, step, defaultValue) {
  stopifnot(step >= 1)
  stopifnot(any(bigwig$chroms == chrom))
  return(.Call(bigWig_chrom_step_sum, bigwig, chrom, step, defaultValue))
}
