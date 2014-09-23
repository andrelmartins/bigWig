\name{probeQuery.bigWig}
\alias{probeQuery.bigWig}
\alias{region.probeQuery.bigWig}
\alias{bed.region.probeQuery.bigWig}
\alias{bed6.region.probeQuery.bigWig}
\alias{step.probeQuery.bigWig}
\alias{bed.step.probeQuery.bigWig}
\alias{bed6.step.probeQuery.bigWig}

\title{Query bigWig Files Containing Probe Data}
\description{Collect values (or transformations of values) from intervals (probes) that overlap the query region(s).}
\usage{

region.probeQuery.bigWig(bw, chrom, start, end, 
                         op = "wavg", abs.value = FALSE, gap.value = NA)

bed.region.probeQuery.bigWig(bw, bed, 
                             op = "wavg", abs.value = FALSE, gap.value = NA)

bed6.region.probeQuery.bigWig(bw.plus, bw.minus, bed6, 
                              op = "wavg", abs.value = FALSE, gap.value = NA)

step.probeQuery.bigWig(bw, chrom, start, end, step, 
                      op = "wavg", abs.value = FALSE, gap.value = NA, 
                      with.attributes = TRUE)

bed.step.probeQuery.bigWig(bw, bed, step, 
                           op = "wavg", abs.value = FALSE, gap.value = NA, 
                          with.attributes = TRUE, as.matrix = FALSE)

bed6.step.probeQuery.bigWig(bw.plus, bw.minus, bed6, step, 
                            op = "wavg", abs.value = FALSE, gap.value = NA, 
                            with.attributes = TRUE, as.matrix = FALSE, 
                            follow.strand = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bw}{either an R object of class 'bigWig' or a character vector containing the prefix and suffix to the path of each bigWig fragment (path = <prefix><chrom><suffix>).}
  \item{bw.plus}{same as 'bw', but for use with plus strand queries.}
  \item{bw.minus}{same as 'bw', but for use with minus strand queries.}
  \item{chrom}{the chromosome to which the query interval belongs to.}
  \item{start}{the start of the query interval (zero-based, left open coordinates like in BED files). In \code{step.probeQuery.bigWig} it can be NULL, indicating a full chromosome query.}
  \item{end}{the end of the query interval (zero-baed, left open coordinates like in BED files). In \code{step.probeQuery.bigWig} it can be NULL, indicating a full chromosome query.}
  \item{bed}{data.frame with 'chrom', 'start' and 'end' coordinates as found in a BED3 file. Further columns are ignored.}
  \item{bed6}{data.frame with 'chrom', 'start', 'end' and 'strand' (column 6) information as found in a BED6 file. Other columns are ignored.}
  \item{step}{step size in base pairs.}
  \item{op}{name of the operation to apply to each query interval: "sum", "avg", "wavg", "min", "max". "wavg" computes the weighted average of overlapping probes, where each probe is weighted by the fraction of the (probe covered) query interval that it covers.}
  \item{abs.value}{logical value; if TRUE, the absolute values of each probe is used as input the interval operation.}
  \item{gap.value}{value assigned to query intervals that have no intercepting probes.}
  \item{as.matrix}{logical value; if TRUE, in 'BED' type queries, the result is a matrix, otherwise it's a list of vectors, one per query 'BED' entry.}
  \item{with.attributes}{logival value; if TRUE, results are annotated with their source coordinates and/or step size.}
  \item{follow.strand}{logical value; if TRUE, results vectors/rows corresponding to minus strand queries are reversed.}
}

\value{Returns either a vector (non-bed functions) or either a list or a matrix (bed/bed6) functions.}

\author{André Luís Martins}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See \code{\link{load.bigWig}} to load a bigWig.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigWig }
\keyword{ query }% __ONLY ONE__ keyword per line
