
foreach.bed <- function(bed, func, envir = parent.frame()) {
  N = dim(bed)[1]

  has.strand = dim(bed)[2] >= 6
  
  for (i in 1:N) {
    chrom = as.character(bed[i, 1])
    start = as.integer(bed[i, 2])
    end = as.integer(bed[i, 3])

    strand = NA
    if (has.strand)
      strand = as.character(bed[i, 6])

    do.call(func, list(i, chrom, start, end, strand), envir=envir)
  }
}
