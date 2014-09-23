\name{profiles}
\alias{quantiles.metaprofile}
\alias{subsampled.quantiles.metaprofile}
\alias{confinterval.metaprofile}
\alias{bootstrapped.confinterval.metaprofile}
\alias{metaprofile.bigWig}
\title{Data Profile Utility Functions}
\description{Functions to create summary profiles out of data matrices.}

\usage{

quantiles.metaprofile(mat, quantiles = c(0.875, 0.5, 0.125))

subsampled.quantiles.metaprofile(mat, quantiles = c(0.875, 0.5, 0.125), fraction = 0.10, 
                                 n.samples = 1000)

confinterval.metaprofile(mat, alpha = 0.05)

bootstrapped.confinterval.metaprofile(mat, alpha = 0.05, n.samples = 300)

metaprofile.bigWig(bed, bw.plus, bw.minus = NULL, step = 1, name = "Signal", 
                   matrix.op = NULL, profile.op = subsampled.quantiles.metaprofile, ...)

}
\arguments{
  \item{mat}{the input data matrix; each row corresponds to a query region, columns to steps.}
  \item{quantiles}{vector of size three with top, middle and bottom quantile breaks to use in creating the summary profile.}
  \item{fraction}{fraction of the data (query regions) to include in each subsample.}
  \item{n.samples}{number of data samples to generate.}
  \item{alpha}{alpha value for confidence intervals (confidence level = 1 - alpha ).}
  \item{bed}{the input BED data.frame defining the set of query regions.}
  \item{bw.plus}{either an R object of class 'bigWig' or a character vector containing the prefix and suffix to the path of each bigWig fragment (path = <prefix><chrom><suffix>).}
  \item{bw.minus}{same as 'bw.plus', but for use with minus strand queries.}
  \item{step}{step size in base pairs.}
  \item{name}{character vector describing the data.}
  \item{matrix.op}{matrix scalling function to apply to the data.}
  \item{profile.op}{summary profile funcion.}
  \item{...}{extra arguments to be passed to matrix.op and/or profile.op.}
}

\value{*.profile functions return a list with 'top', 'middle' and 'bottom' profiles and the 'step' value if present as an attribute in the input matrix. 'step' item is required for the \code{plot.profile.bigWig} function. profile.bigWig returns a list with 'name', 'X0' (= 0), 'step' and the 'top', 'middle' and 'bottom' profiles.
}

\author{André Luís Martins}
