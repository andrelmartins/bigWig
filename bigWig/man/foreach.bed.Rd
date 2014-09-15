\name{foreach.bed}
\alias{foreach.bed}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Foreach loop for BED files (data.frames)}
\description{A quick way to apply functions to every entry in a BED file
  (data.frame). Similar to running \code{apply}, but faster since it's
  taylored to the structure of a BED file. Avoids a lot of boiler plate code.
}
\usage{
foreach.bed(bed, func, envir = parent.frame())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed}{Data frame corresponding to a BED file.}
  \item{func}{Function to apply to each entry in \code{bed}. Function
    must have four arguments: 'index', 'chrom', 'start', 'end' and
    'strand'. Index will be a one-based integer corresponding to the
    current BED line. Chrom is a character string with the chromosome
    name. Start and end are the coordinates for the current entry
    (remember that BED files are zero-based left-open intervals). Strand
    is a character string with the entry's strand (usually '+' or '-')
    or NA if the \code{bed} has less than 6 columns.}
  \item{envir}{Environment where the function is evaluated. Default
    value is \code{parent.frame()} which corresponds to the environment
    where the 'foreach.bed' was called, giving access (through <<-) to
    the local variables.
}
}
\details{Applies the supplied function to each entry in the \code{bed}
  file (actually the data.frame corresponding to the BED file).
}
\value{No return value.}

\author{André Luís Martins}
\note{Remember to use <<- when doing assignments of variables outside
  the loop.}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
#
# Silly example to compute sizes of Bed entries
#
# Obviously "sizes = bed[,3] - bed[,2]" is much faster than this.
#
sizes.bed <- function(bed) {
  N = dim(bed)[1]
  sizes = vector(mode="integer", length=N)

  foreach.bed(bed, function(i, chrom, start, end, strand) {
    sizes[i] <<- end - start
  })

  return(sizes)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigWig }
\keyword{ foreach }
