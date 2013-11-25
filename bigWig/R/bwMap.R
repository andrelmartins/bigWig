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
    error("mappability bigWig file should only have zero or one values (one for unmappable positions)")
  }
  
  res = list(bw = bw, read.len = read.len, read.left.edge = read.left.edge, threshold.fraction = threshold.fraction)
  
  class(res) <- "bwMap"
  return(res)
}

unload.bwMap <- function(bwMap) {
  stopifnot(class(bwMap) != "bwMap")
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
