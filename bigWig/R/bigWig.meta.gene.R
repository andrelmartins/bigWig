subsample.meta.gene <- function(bed, bwPlus, bwMinus, step, half.inner.steps, margin) {

  sizes = bed[,3] - bed[,2]
  min.len = 2*half.inner.steps * step
  good.idxs = which(sizes >= min.len)
  if (min(sizes) < min.len)
    warning("min(sizes) < min.len: ", min(sizes), " ", min.len, ", # = ", length(good.idxs))
  bed = bed[good.idxs,]
  
  N = dim(bed)[1]
  nPermut = 1000
  sampleFrac = 0.1

  windowSize = (2 * margin) %/% step + 2*half.inner.steps
  result = matrix(nrow = nPermut, ncol = windowSize)
  M = as.integer(round(N * sampleFrac, 0))

  # collect values
  values = matrix(nrow = N, ncol = windowSize)
  foreach.bed(bed, function(i, chrom, start, end, strand) {
    bw = bwPlus
    if (strand == '-')
      bw = bwMinus

    # collect left values
    rstart = start - margin
    rend = start + half.inner.steps * step

    left.values = queryByStep.bigWig(bw, chrom, rstart, rend, step, TRUE, default.null=FALSE)

    # collect right values
    rstart = end - half.inner.steps * step
    rend = end + margin
    
    right.values = queryByStep.bigWig(bw, chrom, rstart, rend, step, TRUE, default.null=FALSE)

    # merge
    row.values = c(left.values, right.values)

    if (strand == '-')
      row.values = rev(row.values)

    values[i,] <<- abs(row.values)
  })

  # compute permutations
  for (i in 1:nPermut) {
    idx <- sample(N, size = M, replace = T)
    result[i, ] = colSums(values[idx, ]) / M / step
  }
  ci9 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.875))
  ci1 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.125))
  ci5 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.5))

  res = list(result, ci9, ci1, ci5)
  attr(res, "step") <- step
  attr(res, "inner.steps") <- 2*half.inner.steps
  attr(res, "set.size") <- dim(bed)[1]
  return(res)
}

plot.meta.gene <- function(result, xlab="Distance from TSS (bp)", ylab="Median Signal Density", main = NULL, ylim = NULL, col = 'red', only.median = FALSE, ...) {

  N = length(result[[4]])
  step = attr(result, "step")
  inner.steps = attr(result, "inner.steps")
  set.size = attr(result, "set.size")

  margin.steps = (N - inner.steps) %/% 2

  x = ((1:N) - margin.steps) * step

  if (is.null(ylim)) {
    if (only.median)
      ylim = c(min(result[[4]]), max(result[[4]]))
    else
      ylim = c(min(result[[2]], result[[3]], result[[4]]), 
        max(result[[2]], result[[3]], result[[4]]))
  }
  
  plot(x, result[[4]], type = "l", col = col, lwd = 3, xlab = xlab, 
       ylab = ylab, ylim = ylim, main = main, ...)

  if (!only.median) {
    polygon(c(x, rev(x)), c(result[[2]], rev(result[[3]])), col = "lightgrey", 
            border = NA)
    lines(x, result[[4]], col = col, lwd = 3)
  }
  
  abline(v=0, lty=3, col='black', lwd=3)
  abline(v=inner.steps%/%2 * step, lty=3, col='black', lwd=2)
  abline(v=inner.steps * step, lty=3, col='black', lwd=3)

  mtext(paste("N =", set.size), adj=0, font=2)
}
