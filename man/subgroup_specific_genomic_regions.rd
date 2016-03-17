\name{subgroup_specific_genomic_regions}
\alias{subgroup_specific_genomic_regions}
\title{
find subgroup specific regions
}
\description{
find subgroup specific regions
}
\usage{
subgroup_specific_genomic_regions(gr, factors,
    present = 0.7, absent = 0.3, type = NULL)
}
\arguments{

  \item{gr}{a GRanges object coming from \code{\link{common_regions}}}
  \item{factors}{classes which correspond to samples in \code{gr}}
  \item{present}{how to define a common region is specifically present in one subgroup? The subgroup specificity is calculated based on the precent matrix stored in \code{gr}. For each subgroup which is defined in \code{factors}, if the value is a single numeric value, it means the mean value should be larger than this. It can also be a function for which the input is the vector of percent in corresponding subgroup and the output should be a single logical value.}
  \item{absent}{how to define a common region is specifically absent in one subgroup. Format is same as \code{present}}
  \item{type}{It uses a string containing 1 and 0 to represent types of specificity. E.g. '1100' means present in subgroup 1 and 2 while absent in subgroup 3 and 4. By default, it will output all combination of subgroup specificities.}

}
\value{
In fact, this function split the original \code{gr} into a list in which each element
contains regions corresponding to different subgroup specificity.

The returned value can be sent to \code{\link{plot_subgroup_specificity_heatmap}} to visualize the specificity.
}
\examples{
# There is no example
NULL
}
