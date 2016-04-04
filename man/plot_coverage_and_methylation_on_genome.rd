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

  \item{sid}{a single sample id}
  \item{chromosome}{chromosome}
  \item{species}{species}
  \item{nw}{number of windows}
  \item{...}{pass to \code{\link[gtrellis]{initialize_layout}}}

}
\details{
The whole genome is segented by \code{nw} windows and mean methylation and mean CpG coverage
are visualized as two tracks.
}
\examples{
# There is no example
NULL

}
