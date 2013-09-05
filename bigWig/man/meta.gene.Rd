\name{meta.gene}
\alias{subsample.meta.gene}
\alias{plot.meta.gene}

\title{Meta Plot Functions for Genes}
\description{Functions to collect data and draw meta plots using bigWig
  files over gene bodies.}

\usage{
subsample.meta.gene(bed, bwPlus, bwMinus, step, half.inner.steps, margin)

plot.meta.gene(result, xlab="Distance from TSS (bp)", ylab="Median Signal Density", main = NULL, ylim = NULL, col = 'red', only.median = FALSE, ...)
}

\arguments{
  \item{bed}{a BED file containing the coordinates for the genes to
    be used in generating meta plots.}
  \item{bwPlus}{the bigWig object to be used for forward strand BED
    intervals (or all intervals if the BED file has no strand
    information).}
  \item{bwMinus}{the bigWig object to be used for reverse stand BED
    intervals (only used if the BED file has strand information). If set
    to NULL in conjunction with a \code{bed} argument with strand
    information, \code{collect.many} and functions that use it, will use the
    values from bigWig.plus and reverse the orientation.}
  \item{step}{positive integer value indicating the number of base pairs
    in each plot/query point.}
  \item{half.inner.steps}{how many steps to collect inside the gene
    bodies from each end.}
  \item{margin}{number of bp to collect outside gene bodies.}
   \item{main}{meta plot title}
  \item{ylim}{specify to force a range on the y axis, otherwise scaled
    to fit data}
  \item{col}{color used for the main meta-plot line (median of subsamples).}
  \item{only.median}{only draw median subsample line.}
  \item{...}{additional parameters passed to the main plot function.}
}

\details{
  \code{subsample.meta.gene} collects data for meta-plots by performing
  one thousand 10\% subsamplings of the BED intervals (after culling to
  ensure that the inner steps don't overlap), computing the average
  density per step. Then computes the 0.875, 0.125 and 0.5 quantile
  values over subsamples at each point.

  \code{plot.meta.gene} draws meta-plot with coordinates anchored on
  TSS, shading the area between the 0.875 and 0.125 quantile values (if
  only.median is FALSE). Dashed lines show the two edges and inner
  mid-point.
}

\value{\code{subsample.meta.gene} returns a list with four elements: the value
  corresponding to \code{meta.accum} for each subsample; the values for
  the 0.875 quantile; the values for the 0.125 quantile and the values
  for the 0.5 quantile.
}

\author{André Luís Martins}
