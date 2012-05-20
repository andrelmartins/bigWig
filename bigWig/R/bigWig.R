#
# R interface
#

load.bigWig <- function(filename) {
  res = .Call("bigWig_load", filename)
  class(res) <- "bigWig"
  return(res)
}

unload.bigWig <- function(bigWig) {
  invisible(.Call("bigWig_unload", bigWig))
}
