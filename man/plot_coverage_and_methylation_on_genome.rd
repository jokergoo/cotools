\name{plot_coverage_and_methylation_on_genome}
\alias{plot_coverage_and_methylation_on_genome}
\title{
coverage and methylation for one sample
}
\description{
coverage and methylation for one sample
}
\usage{
plot_coverage_and_methylation_on_genome(sid, chromosome = paste0("chr", 1:22),
    species = "hg19", nw = 10000, ...)
}
\arguments{

  \item{sid}{sample id}
  \item{chromosome}{chromosome}
  \item{species}{species}
  \item{nw}{number of windows}
  \item{...}{pass to \code{\link[gtrellis]{initialize_layout}}}

}
\details{
There will be a track for methylation and one track for coverage.
}
\examples{
# There is no example
NULL
}
