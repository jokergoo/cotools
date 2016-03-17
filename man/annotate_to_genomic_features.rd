\name{annotate_to_genomic_features}
\alias{annotate_to_genomic_features}
\title{
simple annotation to genomic features
}
\description{
simple annotation to genomic features
}
\usage{
annotate_to_genomic_features(gr, genomic_features,
    name = NULL, type = c("percent", "number"), prefix = "overlap_to_", ...)
}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{genomic_features}{a single GRanges object or a list of GRanges objects}
  \item{name}{names for the genomic features if there is no name in \code{genomic_features}}
  \item{type}{How to calculate the values for the annotation. 'number' means numbers of genomic features that each region in \code{gr} overlap; 'percent' means the  percent of each region in \code{gr} that is overlapped by genomic features}
  \item{prefix}{prefix for names of the annotation columns}
  \item{...}{pass to \code{\link[GenomicRanges]{countOverlaps}} or \code{\link{percentOverlaps}}}

}
\details{
it adds new columns in \code{gr} which tell you how \code{gr} is overlaped by each of \code{genomic_features}

Note for the annotation strand information is ignored
}
\examples{
# There is no example
NULL
}
