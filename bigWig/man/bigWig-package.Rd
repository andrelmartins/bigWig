\name{bigWig-package}
\alias{bigWig-package}
\alias{bigWig}
\docType{package}
\title{
bigWig Query Library
~~ bigWig ~~
}
\description{R interface to query bigWig files. C implementation relies
  on UCSC Genome Browser code.
}
\details{
\tabular{ll}{
Copyright: \tab The code in src/jkweb is Copyright (c) 2000-2002 Jim Kent (c) 2003-2014 Regents of the University of California.  All other code is Copyright (c) 2012-2014 Cornell University.\cr
Package: \tab bigWig\cr
Type: \tab Package\cr
Version: \tab 0.2-9\cr
Date: \tab 2014-11-20\cr
License: \tab BSD_3_clause + file LICENSE\cr
LazyLoad: \tab yes\cr
}
bigWig files should be loaded with the \code{load.bigWig} function and
then queries can be conducted with \code{query.bigWig} and related
functions. There is support for loading both local and remote files
(http,ftp). Some utility functions are also provided to build meta plots.
}
\author{
Andre L. Martins, copyright holder Cornell University. src/jkweb by Jim Kent, copyright by Jim Kent and Regents of the University of California.

Maintainer: André Luís Martins <alm253@cornell.edu>
}
\references{
Kent WJ, Zweig AS, Barber G, Hinrichs AS, Karolchik D. BigWig and BigBed: enabling browsing of large distributed data sets. \emph{Bioinformatics}. 2010;26(17):2204-7
}
\keyword{ package }

