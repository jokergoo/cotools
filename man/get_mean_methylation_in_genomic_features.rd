\name{get_mean_methylation_in_genomic_features}
\alias{get_mean_methylation_in_genomic_features}
\title{
calculate average methylation value in a list of regions
}
\description{
calculate average methylation value in a list of regions
}
\usage{
get_mean_methylation_in_genomic_features(sample_id, gf_list, average = TRUE, p = 0.001,
    chromosome = paste0("chr", 1:22),
    filter_fun = function(s) length(s) >= 10 && any(diff(s) < 50))
}
\arguments{

  \item{sample_id}{a list of sample IDs}
  \item{gf_list}{a list of genomic features in \code{\link[GenomicRanges]{GRanges}} class}
  \item{average}{whether to calcualte average methylation in a interval? if not, the function will randomly sample some CpG sites in the interval}
  \item{p}{if average is FALSE, the probability to pick cpg sites}
  \item{chromosome}{chromosome name}
  \item{filter_fun}{filtering function on sites in each intersection}

}
\value{
a list of \code{GRanges} objects.
}
\examples{
# There is no example
NULL
}
