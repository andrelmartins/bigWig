\name{bigWig-deprecated}
\alias{bigWig-deprecated}
\title{Deprecated functions in package \sQuote{bigWig}}

\description{
  These functions are provided for compatibility with older versions
  of \sQuote{bigWig} only, and will be defunct at the next release.
}

\details{

  The following functions are deprecated and will be made defunct; use
  the replacement indicated below:
  \itemize{

    \item{bedQuery.bigWig: \code{\link{bed.region.probeQuery.bigWig}} or \code{\link{bed6.region.probeQuery.bigWig}}
    
    Note: options \code{aggregator="mean", weighted=TRUE} in \code{bedQuery.bigWig} will not yield the same result as \code{op = "wavg"} in the new functions!
    }

  \item{queryByStep.bigWig: \code{\link{step.probeQuery.bigWig}}}
  
  \item{chromStepSum.bigWig: \code{\link{step.bpQuery.bigWig}}}
  }
  
  The following functions are deprecated and will be made defunct. For replacements
  see "profiles" man page:
  \itemize{
    \item{meta.accum(bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS = FALSE, do.sum = FALSE)}
    \item{meta.subsample(bed, bigWig.plus, bigWig.minus, halfWindow, step, at.TSS = FALSE, do.sum = FALSE)}
    \item{meta.subsample.fragmented(bed, bwFolder, bwSuffix, halfWindow, step, at.TSS = FALSE, do.sum = FALSE)}
    \item{meta.plot(result, step, xlab="Distance to center (bp)", ylab="Signal", main=NULL, ylim=NULL, ...)}
    \item{meta.plot.GROseq(result.plus, result.minus, step, xlab="Distance to center (bp)", ylab="Average Reads", main=NULL, ylim = NULL, ...)}
    \item{meta.plot.GROseq.TSS(bed, bigWig.plus, bigWig.minus, halfWindow, step, ...)}
    \item{meta.normalize(result, scaleFactor)}
  }
}
