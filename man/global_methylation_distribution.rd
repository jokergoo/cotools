\name{global_methylation_distribution}
\alias{global_methylation_distribution}
\title{
Global methylation distribution
}
\description{
Global methylation distribution
}
\usage{
global_methylation_distribution(sample_id, annotation,
    annotation_color = structure(seq_along(unique(annotation)), names = unique(annotation)),
    ha = NULL, chromosome = paste0("chr", 1:22), by_chr = FALSE, max_cov = 100,
    background = NULL, p = 0.001)
}
\arguments{

  \item{sample_id}{a vector of sample ids}
  \item{annotation}{subtype information}
  \item{annotation_color}{color for subtypes}
  \item{ha}{additional annotation can be specified as a \code{\link[ComplexHeatmap]{HeatmapAnnotation}} object}
  \item{chromosome}{chromosomes}
  \item{by_chr}{whether make the plot by chromosome}
  \item{max_cov}{maximum coverage (used to get rid of extremely high coverage which affect visualization of CpG coverage distribution)}
  \item{background}{background to look into}
  \item{p}{probability to randomly sample CpG sites}

}
\details{
It visualize distribution of methylation valus and CpG coverages through heatmaps.
}
\examples{
# There is no example
NULL

}
