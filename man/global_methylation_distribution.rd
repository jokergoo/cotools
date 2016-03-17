\name{global_methylation_distribution}
\alias{global_methylation_distribution}
\title{
global methylation distribution
}
\description{
global methylation distribution
}
\usage{
global_methylation_distribution(sample_id, annotation,
    annotation_color = structure(seq_along(unique(annotation)), names = unique(annotation)),
    ha = NULL,
    chromosome = paste0("chr", c(1:22, "X")), by_chr = FALSE, max_cov = 100,
    background = NULL, p = 0.001)
}
\arguments{

  \item{sample_id}{sample id}
  \item{annotation}{subtype}
  \item{annotation_color}{color for subtypes}
  \item{ha}{heatmap annotations}
  \item{chromosome}{chromosome}
  \item{by_chr}{whether by chr}
  \item{max_cov}{maximum coverage}
  \item{background}{background}
  \item{p}{probability to randomly sample CpG sites}

}
\examples{
# There is no example
NULL
}
