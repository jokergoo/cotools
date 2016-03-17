\name{genomicCorr.pintersect}
\alias{genomicCorr.pintersect}
\title{
Intersection between two sets of genomic regions
}
\description{
Intersection between two sets of genomic regions
}
\usage{
genomicCorr.pintersect(query, reference, method = NULL, ...)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{method}{function in which the input is the percentage of each \code{query} that is overlaped in \code{reference} and output is a scalar}
  \item{...}{pass to \code{\link{percentOverlaps}}}

}
\details{
For each region in \code{query}, it calculates the percent that is covered by \code{reference}.

The returned value is percent which is how much \code{query} is covered by \code{reference} (by default).

Be careful with the \code{strand} in your GRanges object!!
}
\examples{
# There is no example
NULL
}
