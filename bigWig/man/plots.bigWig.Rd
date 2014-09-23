\name{plots.bigWig}
\alias{plot.metaprofile}

\title{Plot meta-profile(s)}
\description{Draw meta-profile plot with appropriate error regions.}
\usage{
\method{plot}{metaprofile}(x, minus.profile = NULL, X0 = plus.profile$X0, 
       draw.error = TRUE, col = c("red", "blue", "lightgrey", "lightgrey"), 
       ylim = NULL, xlim = NULL, xlab = "Distance (bp)", ylab = plus.profile$name, ...)
}

\arguments{
\item{x}{Meta-profile instance for sense strand}
\item{minus.profile}{Optional meta-profile instance for reverse strand}
\item{X0}{Numeric offset in base pairs (bp) to shift (subtract) "zero" position.}
\item{draw.error}{Logical value indicating if profile error polygon should be drawn.}
\item{col}{Vector of colors to use for, respectively, sense strand profile line, reverse strand profile line, sense strand error polygon, reverse strand error polygon.}
\item{ylim}{The (y1, y2) limits of the plot.}
\item{xlim}{The (x1, x2) limits of the plot.}
\item{xlab}{Label for x-axis.}
\item{ylab}{Label for y-axis.}
\item{...}{other graphical parameters (see ‘par’).}
}

\author{André Luís Martins}

\seealso{See \code{\link{metaprofile.bigWig}} to generate a meta-profile from a bigWig instance.}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ bigWig }
\keyword{ metaprofile }
\keyword{ plot }
