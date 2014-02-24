#
# Functions to process data matrices into lists of top, middle, bottom vectors
#

quantile.profile <- function(mat, quantiles = c(0.875, 0.5, 0.125)) {
  stopifnot(length(quantiles) == 3)
  stopifnot(all(quantiles < 1 & quantiles > 0))
  
  qTop = quantiles[1]
  qMid = quantiles[2]
  qBottom = quantiles[3]
  
  N = dim(mat)[2]
  
  cTop = sapply(1:N, function(idx) quantile(mat[, idx], qTop))
  cMid = sapply(1:N, function(idx) quantile(mat[, idx], qMid))
  cBottom = sapply(1:N, function(idx) quantile(mat[, idx], qBottom))
  
  return(list(top = cTop, middle = cMid, bottom = cBottom))
}

subsampled.quantile.profile <- function(mat, quantiles = c(0.875, 0.5, 0.125), fraction = 0.10, n.samples = 10000) {
  stopifnot(length(quantiles) == 3)
  stopifnot(all(quantiles < 1 & quantiles > 0))
  
  # create sub-samples
  N = dim(mat)[1]
  M = dim(mat)[2]
  
  result = matrix(nrow=n.samples, ncol=M)
  K = as.integer(round(N * fraction, 0))

  for (i in 1:n.samples) {
    idx <- sample(N, size=K, replace=FALSE)
    
    result[i, ] = colMeans(mat[idx, ], na.rm = TRUE)
  }

  return(quantile.profile(result, quantiles))
}

confint.profile <- function(mat, alpha = 0.05) {
  N = dim(mat)[2]
  cMid = colMeans(mat, na.rm = TRUE)
  sderrs = apply(mat, 2, function(col) {
    sdev = sd(col, na.rm = TRUE)
    M = sum(!is.na(col))
    sdev / sqrt(M)
  })
  
  delta = pnorm(1 - alpha/2) * sderrs
  
  return(list(top = cMid + delta, middle = cMid, bottom = cMid - delta))
}

bootstrapped.confint.profile <- function(mat, alpha = 0.05, n.samples = 300) {
  stat <- function(tbl, idxs) {
    colMeans(tbl[idxs,], na.rm=TRUE)
  }
  
  res = boot(mat, stat, R = n.samples)
  aux = sapply(1:dim(mat)[2], function(idx) {
    tmp = boot.ci(res, type = "norm", index = idx, conf = 1 - alpha)
    c(tmp$normal[,2], tmp$t0, tmp$normal[,3])
  })
  
  return(list(top = aux[3,], middle = aux[2,], bottom = aux[1,]))
}

#
# Convinience function
#

# Profile object structure:
# . name
# . x0 (in units of number of steps)
# . step (size in bp)
# . top, middle, bottom vectors
#

profile.bigWig <- function(bed, bw.plus, bw.minus = NULL, step = 1, name = "Signal", matrix.op = NULL, profile.op = subsampled.quantile.profile, ...) {
  #
  # 1. collect data
  N = dim(bed)[2]
  mat = NULL
  if (N >= 6) {
    stopifnot(bw.minus != NULL)
    mat = bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed, step, abs.value = TRUE, as.matrix = TRUE, follow.strand = TRUE)
  } else {
    if (!is.null(bw.minus))
      warning("bw.minus != NULL but BED contains no strand information")
    mat = bed.step.bpQuery.bigWig(bw.plus, bed, step, as.matrix = TRUE)
  }
  
  #
  # 2. apply matrix transformation
  if (!is.null(matrix.op))
    mat = matrix.op(mat, ...)
  
  #
  # 3. apply profile operation
  result = profile.op(mat, ...)
  
  #
  # 4. create result
  X0 = 0 # can't really tell what X0 was from the input arguments
  
  c(list(name = name, X0 = X0, step = step), result)
}
