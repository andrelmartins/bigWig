#
# bed manipulation functions
#

foreach.bed <- function(bed, func, envir = parent.frame()) {
  if (is.null(bed) || dim(bed)[1] == 0) {
    warning("foreach.bed called with empty bed file.")
    invisible(NULL)
  } else
    invisible(.Call(foreach_bed, bed, func, envir))
}

#
# Bed transformation functions
#

# Arguments: BED (source bed file), upstreamWindow, downstreamWindow
# Output: transformed BED file: interval sizes = upstreamWindow + downstreamWindow + 1

transform.bed <- function(anchor, bed, upstreamWindow, downstreamWindow) {
  start = anchor - upstreamWindow
  end = anchor + downstreamWindow + 1
  
  N = dim(bed)[2]
  if (N >= 6) { # use strands
    is_minus = bed[,6] == '-'
    start[is_minus] = anchor[is_minus] - downstreamWindow
    end[is_minus] = anchor[is_minus] + upstreamWindow + 1
  }
  start = as.integer(start)
  end = as.integer(end)
  
  res = NULL
  if (N > 3)
    res = data.frame(bed[,1], start, end, bed[,4:N])
  else
    res = data.frame(bed[,1], start, end)
  
  colnames(res) <- colnames(bed)
  rownames(res) <- rownames(bed)
  
  return(res)
}

center.bed <- function(bed, upstreamWindow, downstreamWindow) {
  anchor = bed[,2] + ((bed[,3] - bed[,2]) %/% 2)
  
  transform.bed(anchor, bed, upstreamWindow, downstreamWindow)
}

fiveprime.bed <- function(bed, upstreamWindow, downstreamWindow) {
  anchor = bed[,2]
  if (dim(bed)[2] >= 6) {
    is.minus = bed[,6] == '-'
    anchor[is.minus] = bed[is.minus, 3]
  }

  transform.bed(anchor, bed, upstreamWindow, downstreamWindow)
}

threeprime.bed <- function(bed, upstreamWindow, downstreamWindow) {
  anchor = bed[,3] - 1
  if (dim(bed)[2] >= 6) {
    is.minus = bed[,6] == '-'
    anchor[is.minus] = bed[is.minus, 2]
  }
  
  transform.bed(anchor, bed, upstreamWindow, downstreamWindow)
}

downstream.bed <- function(bed, downstreamWindow) {
  N = dim(bed)[2]

  start = bed[,2]
  end = start + downstreamWindow
  
  if (N >= 6) { # use strands
    is_minus = bed[,6] == '-'
    start[is_minus] = bed[is_minus, 3] - downstreamWindow
    end[is_minus] = bed[is_minus, 3]
    
    if (any(start < 0)) {
      warning(sum(start < 0), " rows truncated at zero.")
      start[start < 0] = 0
    }
  }
  
  start = as.integer(start)
  end = as.integer(end)
  
  res = NULL
  if (N > 3)
    res = data.frame(bed[,1], start, end, bed[,4:N])
  else
    res = data.frame(bed[,1], start, end)
  
  colnames(res) <- colnames(bed)
  rownames(res) <- rownames(bed)
  
  return(res)
}

upstream.bed <- function(bed, upstreamWindow) {
  N = dim(bed)[2]

  start = bed[,2] - upstreamWindow
  end = bed[,2]
  
  if (N >= 6) { # use strands
    is_minus = bed[,6] == '-'
    start[is_minus] = bed[is_minus, 3]
    end[is_minus] = bed[is_minus, 3] + upstreamWindow
  }

  if (any(start < 0)) {
    warning(sum(start < 0), " rows truncated at zero.")
    start[start < 0] = 0
  }
  
  start = as.integer(start)
  end = as.integer(end)
  
  res = NULL
  if (N > 3)
    res = data.frame(bed[,1], start, end, bed[,4:N])
  else
    res = data.frame(bed[,1], start, end)
  
  colnames(res) <- colnames(bed)
  rownames(res) <- rownames(bed)
  
  return(res)
}
