\name{extract_sites}
\alias{extract_sites}
\title{
Extract subset of sites which are in a set of intervals
}
\description{
Extract subset of sites which are in a set of intervals
}
\usage{
extract_sites(start, end, site, index = TRUE, filter_fun = NULL)
}
\arguments{

  \item{start}{start position, a vector}
  \item{end}{end position, a vector. Note there should be no overlap between all \code{[start, end]} (You may use \code{\link[IRanges]{reduce}} to merge the overlapping intervals.)}
  \item{site}{positions of all sites, should be sorted.}
  \item{index}{whether return the index in the position vector or just the position itself?}
  \item{filter_fun}{filter sites in the interval. If there are more than one intervals, sites in each interval will be filtered by \code{filter_fun}. Argument \code{s} contains positions of sites in every interval. For example, you can filter intervals if number of sites in  the interval is small (< 10), by \code{filter_fun = function(s) length(s) >= 10}}

}
\details{
Providing a huge vector of genomic positions, we want to extract subset of sites which
locate in a specific region. Normally, we will use:

  \preformatted{
	site = sort(sample(10000000, 1000000))
	start = 123456
	end = 654321
	subsite = site[site >= start & site <= end]  }

Unfortunately, in above code, the whole vector \code{site} will be scaned for four times
(\code{>=}, \code{<=}, \code{&} and \code{[}).
If you want to look for sites in more than one regions (e.g. 1000 regions), in every
loop, the whole \code{site} vector will be re-scanned again and again which is very time-consuming.

Here we have \code{\link{extract_sites}} function which uses binary search to do subsetting.
Of course, \code{site} should be sorted non-decreasing in the first place.

  \preformatted{
	subsite = extract_sites(start, end, site, index = FALSE)  }

Not only for single interval, you can also extract sites in a list of genomic regins,
by setting \code{start} and \code{end} as a vector.

  \preformatted{
	start = c(123456, 234567, 345678)
	end = c(133456, 244567, 355678)
	subsite = extract_sites(start, end, site)  }

\code{filter_fun} can filter sites in each interval by looking at the number of sites in it.
Following code filters out intervals that have less than 10 sites.

  \preformatted{
	subsite = extract_sites(start, end, site, filter_fun = function(x) length(x) >= 10)  }

You can choose to return index only or positions.

  \preformatted{
	subsite = extract_sites(start, end, site, index = FALSE)
	head(subsite)
	subsite_index = extract_sites(start, end, site, index = TRUE)
	head(subsite_index)
	head(site[subsite_index])  }

Follows compare the normal subsetting and binary search.

We first extract all the index that we want to test so that it would not affect the comparison.
In following code, we explicitly transform \code{site} and \code{pos} to \code{double} mode because binary search
is implemented in C-level, and it expects parameter stored in \code{double} mode.

  \preformatted{
	pos = do.call("rbind", lapply(1:1000, function(i) sort(sample(max(site), 2))))
	head(pos)

	t1 = system.time(for(i in 1:1000) \{
		which(site >= pos[i, 1] & site <= pos[i, 2])
	\})
	t1

	site2 = as.double(site)
	pos2 = as.double(pos)
	dim(pos2) = dim(pos)
	t2 = system.time(for(i in 1:1000) \{
		extract_sites(pos2[i, 1], pos2[i, 2], site2)
	\})
	t2  }
}
\value{
a vector of positions or index. If there is no sites in the interval, it will return \code{NULL}.
}
\examples{
# There is no example
NULL

}
