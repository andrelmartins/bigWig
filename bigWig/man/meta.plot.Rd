\name{meta.plot}
\alias{meta.plot}
\alias{meta.plot.GROseq}
\alias{meta.subsample}
\alias{meta.normalize}
\alias{meta.accum}
\title{Meta Plot Functions}
\description{Functions to collect data and draw meta plots using bigWig
  files.}

\usage{
meta.accum(bed, bigWig.plus, bigWig.minus, halfWindow, step)

meta.subsample(bed, bigWig.plus, bigWig.minus, halfWindow, step)

meta.plot(result, step, xlab="Distance to center (bp)", ylab="Signal", main=NULL, ylim=NULL, ...)

meta.plot.GROseq(result.plus, result.minus, step, xlab="Distance to center (bp)", ylab="Average Reads", main=NULL, ylim = NULL, ...)

meta.normalize(result, scaleFactor)
}
\arguments{
  \item{bed}{a BED file containing the coordinates for the elements to
    be used in generating meta plots. The center of each BED interval is
    used as the center of the data window.}
  \item{bigWig.plus}{the bigWig object to be used for forward strand BED
    intervals (or all intervals if the BED file has no strand
    information).}
  \item{bigWig.minus}{the bigWig object to be used for reverse stand BED
    intervals (only used if the BED file has strand information).}
  \item{halfWindow}{half the size of the query window.}
  \item{step}{positive integer value indicating the number of base pairs
    in each plot/query point.}
  \item{result}{the return value from \code{meta.subsample}.}
  \item{result.plus}{the return value from \code{meta.subsample} to be
    used as the sense orientation plot (orientation with respect to BED
    interval elements).}
  \item{result.minus}{the return value from \code{meta.subsample} to be
    used as the anti-sense orientation plot (orientation with respect to BED
    interval elements).}
  \item{scaleFactor}{value to scale the results from a previous call to
    \code{meta.subsample}.}
  \item{main}{meta plot title}
  \item{ylim}{specify to force a range on the y axis, otherwise scaled
    to fit data}
  \item{xlab}{label of the x axis}
  \item{ylab}{label of the y axis}
  \item{...}{additional parameters passed to the main plot function.}
}

\details{
  \code{meta.accum} returns the sum of all query windows.

  \code{meta.subsample} performs one thousand 10\% subsamplings of the
  BED intervals, computing the sum of all query windows for each. Then
  computes the 0.875, 0.125 and 0.5 quantile values at each point.

  \code{meta.plot} plots the results of \code{meta.subsample}, shading
  the area between the 0.875 and 0.125 quantile values.

  \code{meta.plot.GROseq} works similarly to meta.plot, but combines two
  results in the same plot, one for the sense and one for the anti-sense
  strands.

  \code{meta.normalize} can be used to scale the results of
  \code{meta.subsample}.
}

\value{\code{meta.accum} returns a vector with the sum at each window
  step of the mean values in each query window.

  \code{meta.subsample} returns a list with four elements: the value
  corresponding to \code{meta.accum} for each subsample; the values for
  the 0.875 quantile; the values for the 0.125 quantile and the values
  for the 0.5 quantile.
}

\author{André Luís Martins}
