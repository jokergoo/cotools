\name{genomicCorr.reldist}
\alias{genomicCorr.reldist}
\title{
Relative distance between two sets of genomic regions
}
\description{
Relative distance between two sets of genomic regions
}
\usage{
genomicCorr.reldist(query, reference)
}
\arguments{

  \item{query}{genomic region 1, a \code{\link[GenomicRanges]{GRanges}} object}
  \item{reference}{genomic region 2, a \code{\link[GenomicRanges]{GRanges}} object}

}
\details{
For regions in \code{query} and \code{reference}, they are all degenerated as single points
which are the middle points of regions. For each middle point in \code{query}, it looks 
for the nearest point in \code{reference}. The statistics is defined as the ratio of the distance
to nearest point in \code{reference} to the distance of two nearest points in \code{reference} which 
cover the middle point in \code{query}. If \code{query} and \code{reference} are not correlated at all,
It is expected that the ratio follows a uniform distribution. So final statisitics are the KS-statistics
between the real distribution of rations to the uniform distribution.

  \preformatted{
    ----|*************|----- reference
    -----***|--------------- query
         ratio = 3/13  }
}
\references{
Favoriv A, et al. Exploring massive, genome scale datasets with the GenometriCorr package. PLoS Comput Biol. 2012 May; 8(5):e1002529
}
\examples{
# There is no example
NULL
}
