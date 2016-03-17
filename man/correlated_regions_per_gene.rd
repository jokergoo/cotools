\name{correlated_regions_per_gene}
\alias{correlated_regions_per_gene}
\title{
correlated regions in extended gene model
}
\description{
correlated regions in extended gene model
}
\usage{
correlated_regions_per_gene(site, meth, cov, expr, chr, cov_cutoff = 3, min_dp = 4,
    cor_method = "spearman", window_size = 5, factor = NULL, max_width = 10000)
}
\arguments{

  \item{site}{CpG sites}
  \item{meth}{methylation matrix corresponding to \code{site}}
  \item{cov}{coverage}
  \item{expr}{expression for current gene}
  \item{chr}{chromosome}
  \item{cov_cutoff}{cutoff for coverage}
  \item{min_dp}{minimal number of non-NA values for calculating correlations}
  \item{cor_method}{method for calcualting correlations}
  \item{window_size}{how many CpG sites in a window}
  \item{factor}{subtype}
  \item{max_width}{maximum width of a window}

}
\examples{
# There is no example
NULL
}
