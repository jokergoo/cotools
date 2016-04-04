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
Methylation from whole genome bisulfite seuqencing is always huge and it does not
make sense to read them all into the memory. This hook sets how to read the methylation
data and how to return methylation value (e.g. CpG coverage, methylation rate...)

There are following hooks:

\describe{
  \item{set}{set a chromosome. The function accepts a single chromosome name and  returns an object which is used as the first argument in other functions}
  \item{meth}{how to extract methylation value. The function should have three arguments: the object returned from \code{set()}, index of rows and index of columns. Normally, the first argument (\code{obj}) can be ignored when calling this hook. Note the methylation matrix should be column names (used as sample id in other functions)}
  \item{raw}{how to extract raw methylation value, same setting as \code{meth}}
  \item{site}{the function should return a vector of CpG sites}
  \item{coverage}{how to extract CpG coverage, same setting as \code{meth}.}
  \item{GRanges}{howt to extract CpG sites as a \code{\link[GenomicRanges]{GRanges}} object.}
}

Note: positions of CpG sites in a chromosome should be sorted.
}
\examples{
# There is no example
NULL

}
