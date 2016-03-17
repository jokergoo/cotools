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

  \item{sample_id}{sample ids}
  \item{annotation}{annotation of samples}
  \item{annotation_color}{colors}
  \item{chromosome}{chromosome}
  \item{species}{species}
  \item{window_width}{window width}
  \item{style}{style for visualization}
  \item{...}{pass to \code{\link[gtrellis]{initialize_layout}}}

}
\examples{
# There is no example
NULL
}
