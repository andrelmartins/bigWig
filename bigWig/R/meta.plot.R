#
# Charles Danko Meta-plot
#
# Adapted by Andre Martins

collect.counts <- function(bigWig, chrom, start, end, step, do.sum = FALSE) {
  res <- .Call(bigWig_query_by_step, bigWig, chrom, start, end, step, do.sum, 0)
  if (is.null(res))
    return(0)
  return(res)
}


collect.many <- function(bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS = FALSE, do.sum = FALSE) {
  windowSize = (2*halfWindow) %/% step
  midPoint = (bed[,2] + bed[,3]) / 2

  if (at.TSS && dim(bed)[2] >= 6) {
    midPoint = bed[,2]
    neg.idxs = bed[,6] == '-'
    midPoint[neg.idxs] = bed[neg.idxs, 3] - 1
  }
  
  start = (midPoint - halfWindow)
  end = start + windowSize*step

  N = dim(bed)[1]
  result = matrix(nrow=N, ncol=windowSize)
  
  strands = NULL
  if (dim(bed)[2] >= 6)
    strands = as.character(bed[, 6])
  else
    strands = rep("+", N)

  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    strand = strands[i]

    if (strand == "+") {
      bigWig = bigWig.plus
      
      row = collect.counts(bigWig, chrom, start[i], end[i], step, do.sum)

      result[i, ] = abs(row)
    } else {
      bigWig = bigWig.minus
      if (is.null(bigWig.minus))
        bigWig = bigWig.plus
      

      row = collect.counts(bigWig, chrom, start[i], end[i], step, do.sum)

      result[i, ] =  abs(rev(row))
    }
  }

  result
}

meta.accum <- function(bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS = FALSE, do.sum = FALSE) {
  colSums(collect.many(bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS, do.sum))/(dim(bed)[1])
}


meta.subsample <- function(bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS=FALSE, do.sum = FALSE) {
  N = dim(bed)[1]
  
  nPermut = 1000
  sampleFrac = 0.1 # fraction of data per sample

  windowSize = (2*halfWindow) %/% step
  
  result = matrix(nrow=nPermut, ncol=windowSize)
  M = as.integer(round(N * sampleFrac, 0))

  values = collect.many(bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS=at.TSS, do.sum = do.sum)
  
  for (i in 1:nPermut) {
    idx <- sample(N, size=M, replace=T)
    
    result[i, ] = colSums(values[idx, ])/M
  }#renormalize to?? reads in window/windowSize*1000/librarySize

  #
  ## Get CI.
  ci9 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.875))
  ci1 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.125))
  ci5 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.5))

  ## You can then plot ci5; 
  ## ci9 and ci1 are the 75% confidence interval (equivalent to the 'boxes' in a boxplot).

  return(list(result, ci9,ci1,ci5))
}

meta.subsample.fragmented <- function(bed, bwFolder, bwSuffix, halfWindow, step, at.TSS=FALSE, do.sum = FALSE) {
  N = dim(bed)[1]
  
  nPermut = 1000
  sampleFrac = 0.1 # fraction of data per sample

  windowSize = (2*halfWindow) %/% step
  
  result = matrix(nrow=nPermut, ncol=windowSize)
  M = as.integer(round(N * sampleFrac, 0))

  #
  values.chrom = lapply(levels(bed[,1]), function(chrom) {
    bed.chrom = bed[bed[,1] == chrom,]
    filename = paste(bwFolder, "/", chrom, ".", bwSuffix, sep='')

    bwChrom = load.bigWig(filename)

    res = collect.many(bed.chrom, bwChrom, bwChrom, halfWindow, step = step, do.sum = do.sum)

    unload.bigWig(bwChrom)

    return(res)
  })

  # merge all into single big sample
  values = do.call("rbind", values.chrom)
  
  for (i in 1:nPermut) {
    idx <- sample(N, size=M, replace=T)
    
    result[i, ] = colSums(values[idx, ])/M
  }#normalize?? reads in window/windowSize*1000/librarySize


  #
  ## Get CI.
  ci9 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.875))
  ci1 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.125))
  ci5 = sapply(1:windowSize, function(idx) quantile(result[, idx], 0.5))

  ## You can then plot ci5; 
  ## ci9 and ci1 are the 75% confidence interval (equivalent to the 'boxes' in a boxplot).

  return(list(result, ci9,ci1,ci5))
}

meta.plot <- function(result, step, xlab="Distance to center (bp)", ylab="Signal", main=NULL, ylim=NULL, ...) {
  N = length(result[[4]])
  x = ((1:N) - N/2)* step

  # ylim
  if (is.null(ylim))
    ylim = c(min(result[[2]], result[[3]], result[[4]]), max(result[[2]], result[[3]], result[[4]]))

  # establish plot area
  plot(x, result[[4]], type="l", col="red", lwd=3, xlab=xlab, ylab=ylab, ylim=ylim, main=main, ...)

  # draw shade area
  polygon(c(x, rev(x)), c(result[[2]], rev(result[[3]])), col="lightgrey", border=NA)

  # redraw main plot line on top
  lines(x, result[[4]], col="red", lwd=3)
}


meta.plot.GROseq <- function(result.plus, result.minus, step, xlab="Distance to center (bp)", ylab="Average Reads", main=NULL, ylim = NULL, ...) {
  N = length(result.plus[[4]])
  x = ((1:N) - N/2)* step

  # ylim
  if (is.null(ylim))
    ylim = c(min(-result.minus[[2]], -result.minus[[3]], -result.minus[[4]]), max(result.plus[[2]], result.plus[[3]], result.plus[[4]]))

  # establish plot area
  plot(x, result.plus[[4]], type="l", col="red", lwd=3, xlab=xlab, ylab=ylab, ylim=ylim, main=main, ...)

  # draw shade area
  polygon(c(x, rev(x)), c(result.plus[[2]], rev(result.plus[[3]])), col=colors()[419], border=NA)
  polygon(c(x, rev(x)), c(-result.minus[[3]], rev(-result.minus[[2]])), col=colors()[407], border=NA)

  # redraw main plot line on top
  lines(x, result.plus[[4]], col="red", lwd=3)
  lines(x, -result.minus[[4]], col="blue", lwd=3)
  lines(x, rep(0, N), lwd=3, lty=2)
}

meta.plot.GROseq.TSS <- function(bed, bigWig.plus, bigWig.minus, halfWindow, step, ...) {
  sizes = bed[,3] - bed[,2]
  cat("*", sum(sizes > halfWindow), "TSSs selected\n")
  cat("* forward signal ...\n")
  fwd = meta.subsample(bed[sizes > halfWindow,], bigWig.plus, bigWig.minus, halfWindow, step, at.TSS=T)
  cat("* reverse signal ...\n")
  rev = meta.subsample(bed[sizes > halfWindow,], bigWig.minus, bigWig.plus, halfWindow, step, at.TSS=T)
  cat("* ploting ...\n")
  meta.plot.GROseq(fwd, rev, step, ...)
}

pair.ylim <- function(res1, res2) {
  c(min(sapply(c(res1[2:4], res2[2:4]), min)),
    max(sapply(c(res1[2:4], res2[2:4]), max)))
}

quad.ylim <- function(res1.p, res1.m, res2.p, res2.m) {
  c(-max(sapply(c(res1.m[2:4], res2.m[2:4]), max)),
    max(sapply(c(res1.p[2:4], res2.p[2:4]), max)))
}

meta.normalize <- function(result, scaleFactor) {
  lapply(result, function(res) res * scaleFactor)
}

calc.scale.factor <- function(bed, bigWig.p.1, bigWig.m.1, bigWig.p.2, bigWig.m.2) {
  getCounts <- function(bed, bigWig.p, bigWig.m) {
    N = dim(bed)[1]

    counts = vector(mode="numeric", length=N)
    for (i in 1:N) {
      chrom = as.character(bed[i,1])
      start = as.integer(bed[i, 2])
      end = as.integer(bed[i, 3])
      strand = as.character(bed[i, 6])

      bigWig = bigWig.p
      if (strand == "-")
        bigWig = bigWig.m

           bigWig.select(bigWig, chrom)

      data = bigWig.query.region(bigWig, start, end)
      if (!is.null(data))
        counts[i] = abs(sum(data[,3]))
    }

    counts
  }

  counts.1 = getCounts(bed, bigWig.p.1, bigWig.m.1)
  counts.2 = getCounts(bed, bigWig.p.2, bigWig.m.2)

  idxs = counts.1 > 0 & counts.2 > 0
  
  return(counts.1[idxs] / counts.2[idxs])
}
