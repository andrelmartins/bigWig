\name{plots.bigWig}
\alias{plot.track.bigWig}
\title{Plot Functions}
\description{Functions to plot data collected using bigWig functions.}

\usage{
plot.track.bigWig(counts.plus, counts.minus = NULL, newpage = TRUE, ylim
= NULl, fill = c("#E41A1C", "#377EB8"))
}
\arguments{
  \item{counts.plus}{Values collected from a 'step' function
    (\code{step.bpQuery.bigWig} or \code{step.probeQuery.bigWig}) (for
    the plus strand).}
  \item{counts.minus}{Values collected from a 'step' function
    (\code{step.bpQuery.bigWig} or \code{step.probeQuery.bigWig}) (for
    the minus strand).}
  \item{newpage}{Logical value, if TRUE a new page is created before
    drawing.}
  \item{ylim}{numeric vector of length two, which, if not NULL, defines
    the plotting range.}
  \item{fill}{colors to use for each strand.}
}

\details{
  Code relies on 'grid' package to draw graphics. Furthermore, it relies
  on attributes supplied by the source data collection function
  (\code{step.bpQuery.bigWig} or \code{step.probeQuery.bigWig}) to
  establish coordinates.
}

\author{André Luís Martins}
