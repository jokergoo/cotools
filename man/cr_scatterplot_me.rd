\name{cr_scatterplot_me}
\alias{cr_scatterplot_me}
\title{
Scatter plot between methylation and expression
}
\description{
Scatter plot between methylation and expression
}
\usage{
cr_scatterplot_me(cr, expr, gi = NULL, text_column = NULL,
    xlab = "Methylation", ylab = "Expression")
}
\arguments{

  \item{cr}{correlated regions from \code{\link{correlated_regions}}}
  \item{expr}{expression matrix which was used to detect \code{cr}}
  \item{gi}{gene id}
  \item{text_column}{which column in \code{cr} should be put as annotation text in the plot}
  \item{xlab}{xlab in the plot}
  \item{ylab}{ylab in the plot}

}
\details{
Scatterplot for all CRs corresponding to the gene will be made.
}
\examples{
# There is no example
NULL

}
