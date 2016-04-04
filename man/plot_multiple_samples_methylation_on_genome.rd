\name{plot_multiple_samples_methylation_on_genome}
\alias{plot_multiple_samples_methylation_on_genome}
\title{
methylation for more than one samples
}
\description{
methylation for more than one samples
}
\usage{
plot_multiple_samples_methylation_on_genome(sample_id, annotation,
    chromosome = paste0("chr", 1:22), species = "hg19", nw = 1000, ...)
}
\arguments{

  \item{sample_id}{a vector of sample ids}
  \item{annotation}{annotation of samples (e.g. subtypes)}
  \item{chromosome}{chromosome}
  \item{species}{species}
  \item{nw}{number of windows}
  \item{...}{pass to \code{\link[gtrellis]{initialize_layout}}}

}
\details{
The whole genome is segented by \code{nw} windows
}
\examples{
# There is no example
NULL

}
