\name{plot_subgroup_specificity_heatmap}
\alias{plot_subgroup_specificity_heatmap}
\title{
heatmap for visualizing subgroup specific genomic regions
}
\description{
heatmap for visualizing subgroup specific genomic regions
}
\usage{
plot_subgroup_specificity_heatmap(gr_list, genomic_features = NULL)
}
\arguments{

  \item{gr_list}{a list of GRanges which is coming from \code{\link{subgroup_specific_genomic_regions}}}
  \item{genomic_features}{Genomic features that are used for annotatate to regions in \code{gr_list}, it should be a list of GRanges objects}

}
\details{
columns are clustered inside each subgroup.
}
\examples{
# There is no example
NULL
}
