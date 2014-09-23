\name{bwMap}
\alias{load.bwMap}
\alias{unload.bwMap}
\alias{print.bwMap}
\alias{region.bpQuery.bwMap}
\alias{bed6.region.bpQuery.bwMap}
\alias{step.bpQuery.bwMap}
\alias{bed6.step.bpQuery.bwMap}

\title{Query mappability bigWig files}
\description{Collect mappability information for query region(s).}
\usage{

load.bwMap(filename, read.len, read.left.edge, threshold.fraction = 0, udcDir=NULL)

unload.bwMap(bwMap)

region.bpQuery.bwMap(bwMap, chrom, start, end, strand, op = "thresh")

bed6.region.bpQuery.bwMap(bwMap, bed6, op = "thresh")

step.bpQuery.bwMap(bwMap, chrom, start, end, step, strand, 
                   op = "thresh", with.attributes = TRUE)

bed6.step.bpQuery.bwMap(bwMap, bed6, step, 
                        op = "thresh", with.attributes = FALSE, as.matrix = FALSE)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{a character string giving the name of the file to
    load. It can be a valid URL.}
  \item{read.len}{integer number representing the length (in base pairs) of reads for which this mappability file was constructed.}
  \item{read.left.edge}{logical value indicating if a read is represented by it's left-most edge (5' position) or it's right-most edge (3' position).}
  \item{threshold.fraction}{numeric value indicating the maximum fraction of unmappable bases in a query region for that region to still be considered mappable.}
  \item{udcDir}{a character string giving the name of the folder to
    use as a local cache when accessing remote files. Set to NULL to
    use the default location (/tmp/udcCache).
  }
  \item{bwMap}{a bigWig Mappability object}
  \item{chrom}{the chromosome to which the query interval belongs to.}
  \item{start}{the start of the query interval (zero-based, left open coordinates like in BED files). In \code{step.bpQuery.bwMap} it can be NULL, indicating a full chromosome query.}
  \item{end}{the end of the query interval (zero-baed, left open coordinates like in BED files). In \code{step.bpQuery.bwMap} it can be NULL, indicating a full chromosome query.}
  \item{strand}{+ or - character indicating the strand of the supplied coordinates.}
  \item{op}{name of the operation to apply to each query interval: "sum", "avg", "thresh". bigWig intervals are treated as if each position within the interval has the same value as the interval.}
  \item{step}{step size in base pairs.}
  \item{bed6}{data.frame with 'chrom', 'start', 'end' and 'strand' (column 6) information as found in a BED6 file. Other columns are ignored.}
  \item{with.attributes}{logival value; if TRUE, results are annotated with their source coordinates and/or step size.}
  \item{as.matrix}{logical value; if TRUE, in 'BED' type queries, the result is a matrix, otherwise it's a list of vectors, one per query 'BED' entry.}
}

\details{

By convention, a mappability bigWig indicates mappable positions as 0 and unmappable as 1. Raw data in underlying bigWig file marks a position as unmappable if the read that starts at that position is unmappable in the plus strand (obviously dependent on read length).

\code{step.bpQuery.bwMap} can be used to collect chromosome-wide values by passing NULL as the \code{start} and \code{end} arguments.

}
\value{Returns either a vector (non-bed functions) or either a list or a matrix (bed/bed6) functions.}

\author{André Luís Martins}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See \code{\link{load.bigWig}} to load a bigWig.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bwMap }
\keyword{ bpQuery }% __ONLY ONE__ keyword per line
