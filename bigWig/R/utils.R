
foreach.bed <- function(bed, func, envir = parent.frame()) {
  if (is.null(bed) || dim(bed)[1] == 0) {
    warning("foreach.bed called with empty bed file.")
    invisible(NULL)
  } else
    invisible(.Call(foreach_bed, bed, func, envir))
}
