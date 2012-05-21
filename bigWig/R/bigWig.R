#
# R interface
#

load.bigWig <- function(filename, udcDir=NULL) {
  res = .Call("bigWig_load", filename, udcDir)
  class(res) <- "bigWig"
  return(res)
}

unload.bigWig <- function(bigWig) {
  invisible(.Call("bigWig_unload", bigWig))
}

query.bigWig <- function(bigWig, chrom, start, end) {
  res <- .Call("bigWig_query", bigWig, chrom, start, end)
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
