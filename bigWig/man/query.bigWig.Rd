\name{query.bigWig}
\alias{query.bigWig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Query bigWig Files}
\description{Collect the intersection of the query interval with the
  interval-value pairs contained in the bigWig file.}
\usage{
query.bigWig(bigWig, chrom, start, end, clip = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bigWig}{a bigWig object to be queried.}
  \item{chrom}{the chromosome to which the query interval belongs to.}
  \item{start}{the start of the query interval (zero-based, left open
    coordinates like in BED files).}
  \item{end}{the end of the query interval (zero-baed, left open
    coordinates like in BED ifles).}
  \item{clip}{logical value; if TRUE bigWig regions are clipped to the query interval.}
}

\value{Returns a N by 3 matrix with columns (start, end, value) where N
  is the number of intervals contained in the bigWig file that intersect
  with the query interval. Result coordinates are clipped to the query 
  interval (unless \code{clip == FALSE}).}

\author{André Luís Martins}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See \code{\link{load.bigWig}} to load a bigWig.}


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigWig }
\keyword{ query }% __ONLY ONE__ keyword per line
