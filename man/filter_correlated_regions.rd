\name{filter_correlated_regions}
\alias{filter_correlated_regions}
\title{
filter correlation regions
}
\description{
filter correlation regions
}
\usage{
filter_correlated_regions(chromosome = paste0("chr", 1:22), template,
    cutoff = 0.05, adj_method = "BH", meth_diameter_cutoff = 0.25, meth_IQR_cutoff = 0.25,
    anova_cutoff = 0.05)
}
\arguments{

  \item{chromosome}{chromosomes}
  \item{template}{template to find cr files}
  \item{cutoff}{cutoff of adjusted p-values}
  \item{adj_method}{method for calculating adjusted p-values}
  \item{meth_diameter_cutoff}{cutoff for diameters}
  \item{meth_IQR_cutoff}{cutoff for IQR, if there is no subtype information, IQR is used to remove less variable methylation}
  \item{anova_cutoff}{cutoff for ANOVA test}

}
\examples{
# There is no example
NULL

}
