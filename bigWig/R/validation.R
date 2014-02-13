original.expr <- function(arg, n) {
  newname <- substitute(arg)
  
  # Travel up the frame stack until we hit the top.
  for(i in 1:n) {
    oldname <- do.call("substitute", list(as.name(newname), parent.frame(i)))
    newname <- oldname
  }
  deparse(newname)
}

valid.bw <- function(bw) {
  if (!(class(bw) == "bigWig" || (is.character(bw) && length(bw) == 2 && (!any(is.na(bw) | is.null(bw)))))) {
    stop("invalid bigWig: ", original.expr(bw, 2))
  }
}

valid.chrom <- function(bw, chrom) {
  stopifnot(is.character(chrom))
  
  if (class(bw) == "bigWig") {
    if (!(chrom %in% bw$chroms))
    stop("bigWig does not contain information on chromosome: ", chrom)
    return(TRUE)
  } else if (is.character(bw) && length(bw) == 2) {
    path = paste(bw[1], chrom, bw[2], sep='')
    if (file.exists(path))
    return(TRUE)
    stop("can't find bigWig file for chromosome: ", chrom)
  } else
  stop("invalid bigWig object")
}

valid.query.range <- function(start, end, index = NA, step = NA) {
  if (end <= start) {
    if (!is.na(index))
    stop("bed:", index, ": end must be > start: ", start, ", ", end)
    else
    stop("end must be > start: ", start, ", ", end)
  }
  if (start < 0 || end < 1) {
    if (!is.na(index))
    stop("bed:", index, ": start and end must be positive: ", start, ", ", end)
    else
    stop("start and end must be positive: ", start, ", ", end)
  }
  
  if (!is.na(step) && (end - start) %/% step == 0) {
    if (!is.na(index))
    stop("bed:", index, ": query region is less than a step wide")
    else
    stop("query region is less than a step wide")
  }
  
  if (!is.na(step) && (end - start) %% step > 0) {
    if (!is.na(index))
    warning("bed:", index, ": query region is not an integer multiple of step")
    else
    warning("query region is not an integer multiple of step")
  }
}

valid.strand <- function(strand) {
  is.character(strand) & (strand == "+" | strand == "-")
}

valid.probe.op <- function(op) {
  if (all(op != c("sum", "avg", "wavg", "min", "max")))
  stop("invalid probe operation: ", op)
}

valid.bp.op <- function(op) {
  if (all(op != c("sum", "avg", "min", "max")))
  stop("invalid base pair operation: ", op)
}
