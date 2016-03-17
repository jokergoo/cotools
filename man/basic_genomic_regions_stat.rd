\name{basic_genomic_regions_stat}
\alias{basic_genomic_regions_stat}
\title{
Basic statistics on genomic regions
}
\description{
Basic statistics on genomic regions
}
\usage{
basic_genomic_regions_stat(gr_list, annotation = NULL, annotation_color = NULL,
    main = NULL, species = "hg19", type = c("proportion", "number", "median_width"),
    by_chr = FALSE)
}
\arguments{

  \item{gr_list}{a list of \code{\link[GenomicRanges]{GRanges}}.}
  \item{annotation}{a vector which contains levels of samples, better have names which correspond to the names of \code{gr_list}}
  \item{annotation_color}{colors corresponding to levels of annotations}
  \item{main}{title of the plot}
  \item{species}{species, necessary if \code{type} equals to \code{proportion}.}
  \item{type}{type of statistics}
  \item{by_chr}{take all chromosomes as a whole or calculate statistics for every chromosome?}

}
\details{
For \code{type} settings:

\describe{
  \item{proportion}{proportion of total length of regions compared to the whole genome}
  \item{number}{number of regions}
  \item{median_width}{median width of regions}
}
}
\value{
A data frame which contains statistics for each chromosome in each sample.
}
\examples{
# There is no example
NULL
}
