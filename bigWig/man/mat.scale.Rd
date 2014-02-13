\name{mat.scale}
\alias{rpkm.scale}
\alias{densityToOne.scale}
\alias{maxToOne.scale}
\alias{zeroToOne.scale}
\title{Data Matrix Scalling Utility Functions}
\description{Functions to scale data matrices.}

\usage{

rpkm.scale(mat, step, libSize)

densityToOne.scale(mat, na.on.zero = TRUE)

maxToOne.scale(mat)

zeroToOne.scale(mat)

}
\arguments{
  \item{mat}{the input data matrix; each row corresponds to a query region, columns to steps}
  \item{step}{step size in base pairs}
  \item{libSize}{total library mapped read count}
  \item{na.on.zero}{logical indicating if steps with zero counts should be marked as NA}
}

\details{
  "RPKM stands for Reads Per Kilobase of transcript per Million mapped reads. FPKM stands for Fragments Per Kilobase of transcript per Million
  mapped reads. " from \url{http://cufflinks.cbcb.umd.edu/faq.html}.
}

\value{Data matrix scaled according to the requested method.
}

\author{André Luís Martins}
