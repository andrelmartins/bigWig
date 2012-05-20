#
# R interface
#

load.bigWig <- function(filename, udcDir=NULL) {
  res = .Call("bigWig_load", filename, udcDir)
  class(res) <- "bigWig"
  return(res)
}

unload.bigWig <- function(bigWig) {
  invisible(.Call("bigWig_unload", bigWig))
}
