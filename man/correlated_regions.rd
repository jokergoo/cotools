\name{correlated_regions}
\alias{correlated_regions}
\title{
correlated regions
}
\description{
correlated regions
}
\usage{
correlated_regions(sample_id, expr, txdb, chr, extend = 50000,
    cov_filter = function(x) sum(x > 0) > length(x)/2,
    cor_method = "spearman", factor = NULL, window_size = 5, max_width = 10000,
    raw_meth = FALSE, cov_cutoff = 3, min_dp = 4, col = NULL)
}
\arguments{

  \item{sample_id}{sample id}
  \item{expr}{expression matrix}
  \item{txdb}{\code{transcritpDb} object}
  \item{chr}{chromosome}
  \item{extend}{extension of gene model, both upstream and downstream}
  \item{cov_filter}{function to filter on coverage}
  \item{cor_method}{method to calculate correlation}
  \item{factor}{subtype}
  \item{window_size}{how many CpGs in a window}
  \item{max_width}{maximum width of a window}
  \item{raw_meth}{whether use raw methylation value (unsmoothed)}
  \item{cov_cutoff}{cutoff for coverage}
  \item{min_dp}{minimal non-NA values for calculating correlations}
  \item{col}{color for subtypes}

}
\details{
based on \code{\link{correlated_regions_per_gene}}
}
\examples{
# There is no example
NULL
}
