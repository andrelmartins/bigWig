\name{load.bigWig}
\alias{load.bigWig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Load bigWig Files}
\description{Load UCSC bigWig files. Only the file header is loaded into
  memory, enough to perform subsequent queries. As with bigWig command
  line tools, the file can be on a remote server.}
\usage{
load.bigWig(filename, udcDir = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{a character string giving the name of the file to
    load. It can be a valid URL.}
  \item{udcDir}{a character string giving the name of the folder to
    use as a local cache when accessing remote files. Set to NULL to
    use the default location (/tmp/udcCache).
    }
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
\author{André Luís Martins}
\note{Return value contains a pointer to a C allocated structure stored
  as an attribute (handle_pte). Do not change this value directly!©
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See \code{\link{unload.bigWig}} to safely unload a bigWig}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (filename) 
{
    res = .Call("bigWig_load", filename)
    class(res) <- "bigWig"
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigWig }
\keyword{ load }
\keyword{ read }
