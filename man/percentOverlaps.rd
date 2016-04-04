\name{percentOverlaps}
\alias{percentOverlaps}
\title{
For every interval in \code{query}, percent that is covered by \code{subject}
}
\description{
For every interval in \code{query}, percent that is covered by \code{subject}
}
\usage{
percentOverlaps(query, subject, ...)
}
\arguments{

  \item{query}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{subject}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{...}{pass to \code{\link[GenomicRanges]{findOverlaps}}}

}
\value{
a numeric vector in which every element correspond to one interval in \code{query}.

Be careful with \code{strand} in your \code{\link[GenomicRanges]{GRanges}} object!!
}
\examples{
# There is no example
NULL

}
