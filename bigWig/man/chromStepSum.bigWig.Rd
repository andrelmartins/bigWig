\name{chromStepSum.bigWig}
\alias{chromStepSum.bigWig}
\alias{queryByStep.bigWig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Stepped Query bigWig Files}
\description{Collect stepped summary values across a chromosome or region.}
\usage{
chromStepSum.bigWig(bigwig, chrom, step, defaultValue)

queryByStep.bigWig(bigWig, chrom, start, end, step, do.sum = F) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bigwig}{a bigWig object to be queried.}
  \item{chrom}{the chromosome to query.}
  \item{step}{the size of the query blocks.}
  \item{defaultValue}{the default value to return for blocks with no
    data.}
  \item{start}{the start of the query interval (zero-based, left open
    coordinates like in BED files).}
  \item{end}{the end of the query interval (zero-baed, left open
    coordinates like in BED ifles).}
  \item{do.sum}{indicates if the value per query step should be the sum
    or the average (default).}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{chromStepSum.bigWig} Returns a vector of size \code{chromSize
    \%/\% step} filled with the sums of the per base values in each
  block step.
  \item{queryByStep.bigWig} Returns a vector of the size \code{(end -
    start) \%/\% step} filled with the sums of the per bigWig interval values in each block step.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{André Luís Martins}
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See \code{\link{load.bigWig}} to load a bigWig. See
  \code{\link{query.bigWig}} to query a single interval. See
  \code{\link{bedQuery.bigWig}} to query a BED file.}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigWig }
\keyword{ query }
\keyword{ chrom }
