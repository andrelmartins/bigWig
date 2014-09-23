\name{bpQuery.bigWig}
\alias{bpQuery.bigWig}
\alias{region.bpQuery.bigWig}
\alias{bed.region.bpQuery.bigWig}
\alias{bed6.region.bpQuery.bigWig}
\alias{step.bpQuery.bigWig}
\alias{bed.step.bpQuery.bigWig}
\alias{bed6.step.bpQuery.bigWig}

\title{Query bigWig Files Containing Base Pair Data}
\description{Collect values (or transformations of values) from intervals that overlap the query region(s).}

\usage{
region.bpQuery.bigWig(bw, chrom, start, end,
                      strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, 
                      bwMap = NULL)

bed.region.bpQuery.bigWig(bw, bed, 
                          strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, 
                          bwMap = NULL)

bed6.region.bpQuery.bigWig(bw.plus, bw.minus, bed6, 
                           op = "sum", abs.value = FALSE, gap.value = 0, bwMap = NULL)

step.bpQuery.bigWig(bw, chrom, start, end, step, 
                    strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, 
                    bwMap = NULL, with.attributes = TRUE)

bed.step.bpQuery.bigWig(bw, bed, step, 
                        strand = NA, op = "sum", abs.value = FALSE, gap.value = 0, 
                        bwMap = NULL, with.attributes = TRUE, as.matrix = FALSE)

bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step, 
                         op = "sum", abs.value = FALSE, gap.value = 0, 
                         bwMap = NULL, with.attributes = TRUE, as.matrix = FALSE, 
                         follow.strand = FALSE)

}

\arguments{
  \item{bw}{either an R object of class 'bigWig' or a character vector containing the prefix and suffix to the path of each bigWig fragment (path = <prefix><chrom><suffix>).}
  \item{bw.plus}{same as 'bw', but for use with plus strand queries.}
  \item{bw.minus}{same as 'bw', but for use with minus strand queries.}
  \item{chrom}{the chromosome to which the query interval belongs to.}
  \item{start}{the start of the query interval (zero-based, left open coordinates like in BED files). In \code{step.bpQuery.bigWig} it can be NULL, indicating a full chromosome query.}
  \item{end}{the end of the query interval (zero-baed, left open coordinates like in BED files). In \code{step.bpQuery.bigWig} it can be NULL, indicating a full chromosome query.}
  \item{strand}{+ or - character indicating the strand of the supplied coordinates.}
  \item{bed}{data.frame with 'chrom', 'start' and 'end' coordinates as found in a BED3 file. Further columns are ignored.}
  \item{bed6}{data.frame with 'chrom', 'start', 'end' and 'strand' (column 6) information as found in a BED6 file. Other columns are ignored.}
  \item{step}{step size in base pairs.}
  \item{op}{name of the operation to apply to each query interval: "sum", "avg", "wavg", "min", "max". bigWig intervals are treated as if each position within the interval has the same value as the interval.}
  \item{abs.value}{logical value; if TRUE, the absolute values of each probe is used as input the interval operation.}
  \item{gap.value}{value assigned to query intervals that have no recorded values.}
  \item{as.matrix}{logical value; if TRUE, in 'BED' type queries, the result is a matrix, otherwise it's a list of vectors, one per query 'BED' entry.}
  \item{with.attributes}{logival value; if TRUE, results are annotated with their source coordinates and/or step size.}
  \item{bwMap}{R object of class 'bwMap' containing the information required to call steps or regions as unmappable. **not implemented yet**}
  \item{follow.strand}{logical value; if TRUE, results vectors/rows corresponding to minus strand queries are reversed.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
\code{step.bpQuery.bigWig} can be used to, efficiently, collect chromosome-wide values by passing NULL as the \code{start} and \code{end} arguments.

}
\value{Returns either a vector (non-bed functions) or either a list or a matrix (bed/bed6) functions.}

\author{André Luís Martins}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See \code{\link{load.bigWig}} to load a bigWig.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigWig }
\keyword{ query }% __ONLY ONE__ keyword per line
