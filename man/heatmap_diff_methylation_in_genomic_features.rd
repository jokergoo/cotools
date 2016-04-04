\name{heatmap_diff_methylation_in_genomic_features}
\alias{heatmap_diff_methylation_in_genomic_features}
\title{
Heatmap for mean methylation in genomic features
}
\description{
Heatmap for mean methylation in genomic features
}
\usage{
heatmap_diff_methylation_in_genomic_features(gr, annotation,
    annotation_color = structure(seq_along(unique(annotation)), names = unique(annotation)),
    txdb = NULL, gf_list = NULL, gf_type = "percent",
    min_mean_range = 0.2, cutoff = 0.05, adj_method = "BH", title = NULL,
    cluster_cols = c("subgroup", "all", "none"), ha = NULL, ...)
}
\arguments{

  \item{gr}{object returned by \code{\link{get_mean_methylation_in_genomic_features}}}
  \item{annotation}{subtype of samples}
  \item{annotation_color}{colors of subtypes}
  \item{txdb}{A \code{\link[GenomicFeatures]{TxDb}} object}
  \item{gf_list}{a list of genomic features which are used as row annotations}
  \item{gf_type}{how to overlap genomic features}
  \item{min_mean_range}{minimal range between mean value in subtypes}
  \item{cutoff}{if subtype information is provided, p-value for the oneway ANOVA test}
  \item{adj_method}{how to calculate adjusted p-values}
  \item{title}{title of the plot}
  \item{cluster_cols}{how to cluster columns}
  \item{ha}{column annotations added to the heatmap}
  \item{...}{pass to \code{\link[ComplexHeatmap]{Heatmap}}}

}
\examples{
# There is no example
NULL

}
