\name{genomicCorr.sintersect}
\alias{genomicCorr.sintersect}
\title{
Intersection between two sets of genomic regions
}
\description{
Intersection between two sets of genomic regions
}
\usage{
genomicCorr.sintersect(query, reference, restrict = NULL)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{restrict}{subset of sites that should be only looked into, a \code{\link[GenomicRanges]{GRanges}} object}

}
\details{
It calculates the total length of overlapped regions in \code{query}.

If the interest is e.g. the number of CpG sites both in \code{query} and in \code{reference}
\code{restrict} can be set with a GRanges object which contains positions of CpG sites.

Be careful with the \code{strand} in your GRanges object!!
}
\examples{
# There is no example
NULL

}
