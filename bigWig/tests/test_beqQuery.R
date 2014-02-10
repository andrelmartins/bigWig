library(bigWig)

fpath = system.file("extdata", "test3.bigWig", package="bigWig")
tmp = load.bigWig(fpath)

bedQuery_.bigWig <- function(bed, bwPlus, bwMinus = NULL, gapValue = NULL, weighted = F, aggregator = "sum") {
  op = aggregator
  if (op == "mean") {
    if (weighted) {
      stop("old behaviour for weighted averages is no longer supported!")
    } else
      op = "avg"
  }
  
  gv = NA
  if (!is.null(gapValue))
    gv = gapValue

  if (dim(bed)[2] >= 6 && !is.null(bwMinus)) {
    bed6.region.probeQuery.bigWig(bwPlus, bwMinus, bed, op, abs.value = FALSE, gap.value = gv)
  } else {
    bed.region.probeQuery.bigWig(bwPlus, bed, op, abs.value = FALSE, gap.value = gv)
  }
}

for (op in c("sum", "min", "max", "mean")) {
  for (gv in c(NA, 0, 1, 3)) {
    vOld = bedQuery.bigWig(data.frame("chrU", c(1, 15, 100), c(10, 45, 110)), tmp, aggregator=op, weighted=FALSE, gapValue = gv)
    vNew = bedQuery_.bigWig(data.frame("chrU", c(1, 15, 100), c(10, 45, 110)), tmp, aggregator=op, weighted=FALSE, gapValue = gv)

    stopifnot(all.equal(vOld, vNew))
  }
}
