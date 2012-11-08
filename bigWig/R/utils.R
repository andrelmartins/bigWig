
foreach.bed <- function(bed, func, envir = parent.frame()) {
  invisible(.Call(foreach_bed, bed, func, envir))
}
