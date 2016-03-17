\name{wgbs_qcplot}
\alias{wgbs_qcplot}
\title{
basic qc plot for methylation
}
\description{
basic qc plot for methylation
}
\usage{
wgbs_qcplot(sample_id, chromosome = paste0("chr", 1:22))
}
\arguments{

  \item{sample_id}{a single sample id}
  \item{chromosome}{chromosome}

}
\details{
it will produce five plots

10000 CpG sites are randomly sampled to make the plot
}
\examples{
# There is no example
NULL
}
