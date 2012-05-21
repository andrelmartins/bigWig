\name{unload.bigWig}
\alias{unload.bigWig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Unload bigWig Files}
\description{Release memory associated with the bigWig object.}
\usage{
unload.bigWig(bigWig)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{bigWig}{a bigWig object}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{André Luís Martins
%%  ~~who you are~~
}
\note{}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See \code{\link{load.bigWig}} to load a bigWig.}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (bigWig) 
{
    invisible(.Call("bigWig_unload", bigWig))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigwig }
\keyword{ unload }
\keyword{ free }
