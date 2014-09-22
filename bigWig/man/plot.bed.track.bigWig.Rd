\name{plot.bed.track.bigWig}
\alias{plot.bed.track.bigWig}
\title{Plot Functions}
\description{Function to plot a BED track.}

\usage{
plot.bed.track.bigWig(bed, chrom, start, end, newpage = TRUE, margin = unit(4, "pt"), height = NULL, draw.direction = TRUE, score.to.color = NULL, mark.color = "white")
}
\arguments{
  \item{bed}{data.frame containing BED regions.}
  \item{chrom}{character value with chromosome name for region to be plotted.}
  \item{start}{integer value with (zero based) start coordinate for region to be plotted.}
  \item{end}{integer value with (one based) end coordinate for region to be plotted.}
  \item{newpage}{Logical value, if TRUE a new page is created before drawing.}
  \item{margin}{unit value of margin between rows of BED regions.}
  \item{height}{unit value of height of a BED region.}
  \item{draw.direction}{logical value, if TRUE, arrows are drawn over BED regions to indicate direction (if room is avalable and BED has strand information).}
  \item{score.to.color}{function that transforms scores to R colors, if present (in 5th column).}
  \item{mark.color}{R color used to draw directional and edge clip marks.}
}

\details{
  Code relies on 'grid' package to draw graphics. Default region color is 'black'.
}

\author{André Luís Martins}
