\name{add_subtype_specificity}
\alias{add_subtype_specificity}
\title{
ad subtype specificity columns
}
\description{
ad subtype specificity columns
}
\usage{
add_subtype_specificity(cr, cutoff = 0.05, suffix = "_ss")
}
\arguments{

  \item{cr}{cr}
  \item{cutoff}{cutoff for ANOVA test}

}
\details{
1 is defined as the methylation is higher than all other subtypes and the difference is significant.
-1 is defined as the methylation is lower than all other subtypes and the difference is significant.
All the others are defined as 0.
}
\examples{
# There is no example
NULL
}
