\name{methylation_hooks}
\alias{methylation_hooks}
\title{
Hook functions to extract methylation
}
\description{
Hook functions to extract methylation
}
\usage{
methylation_hooks(..., RESET = FALSE, READ.ONLY = NULL, LOCAL = FALSE)
}
\arguments{

  \item{...}{Arguments for the parameters, see "details" section}
  \item{RESET}{reset to default values}
  \item{READ.ONLY}{whether only return read-only options}
  \item{LOCAL}{switch local mode}

}
\details{
there are following hooks:

\describe{
  \item{set}{set a chromosome. The function accepts a single chromosome name and  returns an object which is used as the first argument in other functions}
  \item{meth}{how to extract methylation value. The function should have three arguments, the object returned from \code{set()}, index of rows and index of columns}
  \item{raw}{how to extract row methylation value}
  \item{site}{how to extract site}
  \item{coverage}{how to extract coverage}
  \item{GRanges}{howt to extract sites as a GRanges object}
}
}
\examples{
# There is no example
NULL
}
