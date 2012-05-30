\name{bedQuery.bigWig}
\alias{bedQuery.bigWig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Query bigWig Files for BED intervals}
\description{Collect min/max/sum/mean values for each interval in the
  BED file.}
\usage{
bedQuery.bigWig(bed, bwPlus, bwMinus = NULL, gapValue = NULL, weighted = F, aggregator = "sum")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bed}{a \code{data.frame} containing the information from a BED
    file. If it has 6 columns and bwMinus is supplied, then the queries
    will use the appropriate bigWig file (bwPlus for '+' strand, bwMinus
    for '-' strand).}
  \item{bwPlus}{a bigWig object to be queried.}
  \item{bwMinus}{a bigWig object to be queried when BED intervals are in
    the '-' strand.}
  \item{gapValue}{value to report when no bigWig intervals fall within a
    given BED interval (otherwise \code{NA}) and to use for gaps when
    \code{aggregator = "mean"} and \code{weighted = TRUE}.}
  \item{weighted}{indicates whether bigWig interval values should be
    weighted by the fraction of the query interval that they account for
    when \code{agregator = "mean"}.}
  \item{aggregator}{choice of aggregation function to use on each BED
    interval. Possible values are "min", "max", "sum" and "mean".}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Returns a numeric vector with length equal to the number of
  intervals in the BED file.}
\references{
%% ~put references to the literature/web site here ~
}
\author{André Luís Martins}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See \code{\link{load.bigWig}} to load a bigWig. See
  \code{\link{query.bigWig}} to query a single interval.}
\examples{


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigWig }
\keyword{ query }
\keyword{ bed }

