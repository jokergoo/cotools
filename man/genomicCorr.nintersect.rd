\name{genomicCorr.nintersect}
\alias{genomicCorr.nintersect}
\title{
Intersections between two sets of genomic regions
}
\description{
Intersections between two sets of genomic regions
}
\usage{
genomicCorr.nintersect(query, reference, ...)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{...}{pass to \code{\link[GenomicRanges]{countOverlaps}}}

}
\details{
It calculates number of regions in \code{query} that overlap with \code{reference}.

Please note this value is not equal to the number of intersections betweenn two sets of regions,
because one region in \code{query} may overlap with more than one
regions in \code{reference}

Be careful with the \code{strand} in your GRanges object!!
}
\examples{
# There is no example
NULL

}
