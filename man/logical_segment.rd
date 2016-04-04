\name{logical_segment}
\alias{logical_segment}
\title{
segment by continuous logical values
}
\description{
segment by continuous logical values
}
\usage{
logical_segment(l)
}
\arguments{

  \item{l}{logical vector}

}
\details{
the logical vector will be segmented according to their values.
It returns intervals for continuous \code{\link{TRUE}} values
}
\value{
a data frame in which the first column is the index of start sites
the second column is the index of end sites
}
\examples{
# There is no example
NULL

}
