\name{load.bigWig}
\alias{load.bigWig}

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

\value{
Returns an object of the bigWig class, which contains a C pointer (see note bellow) and the following accessible parts:
\item{basesCovered}{Number of bp covered by bigWig data.}
\item{chromSizes}{Integer vector with number of bases per chromosome (same order as chroms vector.}
\item{chroms}{Character vector with chromosome names for which there is data in the bigWig file.}
\item{isCompressed}{Logical value indicating state of internal bigWig file compression.}
\item{isSwapped}{Logical value indicating if binary representation is byte swapped (big vs little endian).}
\item{max}{Numeric value containing the maximum value stored in the bigWig file.}
\item{mean}{Numeric value containing the mean value across all covered bases in the bigWig file.}
\item{min}{Numeric value containing the minimum value stored in the bigWig file.}
\item{primaryDataSize}{Integer number of bytes used by data.}
\item{primaryIndexSize}{Integer number of bytes used by index.}
\item{std}{Standard deviation of bigWig values over covered bases.}
\item{version}{Integer version number for bigWig file format.}
\item{zoomLevels}{Integer number of zoom level supported by bigWig file.}
}
\references{
Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed data sets. \emph{Bioinformatics}. 2010;26(17):2204-7
}
\author{André Luís Martins}
\note{Return value contains a pointer to a C allocated structure stored
  as an attribute (handle_pte). Do not change this value directly!
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See \code{\link{unload.bigWig}} to safely unload a bigWig}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigWig }
\keyword{ load }
\keyword{ read }
