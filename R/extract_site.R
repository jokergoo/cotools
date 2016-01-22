
# == title
# Extract subset of sites which are in a set of intervals
#
# == param
# -start      start position, a vector
# -end        end position, a vector. Note there should be no overlap between all ``[start, end]``
#             (You may use `IRanges::reduce` to merge the overlapping intervals.)
# -site       positions of all sites, should be sorted.
# -index      whether return the index in the position vector or just the position itself?
# -filter_fun filter sites in the interval. If there are more than one intervals, sites
#             in each interval will be filtered by ``filter_fun``. Argument ``s`` contains positions of sites in every interval.
#             For example, you can filter intervals if number of sites in 
#             the interval is small (< 10), by ``filter_fun = function(s) length(s) >= 10``
#
# == details
# Providing a huge vector of genomic positions, we want to extract subset of sites which
# locate in a specific region. Normally, we will use:
#
# 	site = sort(sample(10000000, 1000000))
# 	start = 123456
# 	end = 654321
# 	subsite = site[site >= start & site <= end]
#
# Unfortunately, in above code, the whole vector ``site`` will be scaned for four times
# (``>=``, ``<=``, ``&`` and ``[``).
# If you want to look for sites in more than one regions (e.g. 1000 regions), in every
# loop, the whole ``site`` vector will be re-scanned again and again which is very time-consuming.
#
# Here we have `extract_sites` function which uses binary search to do subsetting.
# Of course, ``site`` should be sorted non-decreasing in the first place.
#
# 	subsite = extract_sites(start, end, site, index = FALSE)
#
# Not only for single interval, you can also extract sites in a list of genomic regins,
# by setting ``start`` and ``end`` as a vector.
#
# 	start = c(123456, 234567, 345678)
# 	end = c(133456, 244567, 355678)
# 	subsite = extract_sites(start, end, site)
#
# ``filter_fun`` can filter sites in each interval by looking at the number of sites in it.
# Following code filters out intervals that have less than 10 sites.
#
# 	subsite = extract_sites(start, end, site, filter_fun = function(x) length(x) >= 10)
#
# You can choose to return index only or positions.
#
# 	subsite = extract_sites(start, end, site, index = FALSE)
# 	head(subsite)
# 	subsite_index = extract_sites(start, end, site, index = TRUE)
# 	head(subsite_index)
# 	head(site[subsite_index])
#
# Follows compare the normal subsetting and binary search.
#
# We first extract all the index that we want to test so that it would not affect the comparison.
# In following code, we explicitly transform ``site`` and ``pos`` to ``double`` mode because binary search
# is implemented in C-level, and it expects parameter stored in ``double`` mode.
#
# 	pos = do.call("rbind", lapply(1:1000, function(i) sort(sample(max(site), 2))))
# 	head(pos)
#
# 	t1 = system.time(for(i in 1:1000) {
# 		which(site >= pos[i, 1] & site <= pos[i, 2])
# 	})
# 	t1
#
# 	site2 = as.double(site)
# 	pos2 = as.double(pos)
# 	dim(pos2) = dim(pos)
# 	t2 = system.time(for(i in 1:1000) {
# 		extract_sites(pos2[i, 1], pos2[i, 2], site2)
# 	})
# 	t2
#
#
# == value
# a vector of positions or index. If there is no sites in the interval, it will return ``NULL``.
#
extract_sites = function(start, end, site, index = TRUE, filter_fun = NULL) {
	if(!is.double(site)) site = as.double(site)
	if(!is.double(start)) start = as.double(start)
	if(!is.double(end)) end = as.double(end)
	
	.extract = function(start, end) {
		n = length(site)
		
		# --------o------o---- site
		# --+--+--------------
		if(end < site[1]) return(NULL)
		
		# --------o------o---- site
		# -----------------+-+
		if(start > site[n]) return(NULL)
		
		# --------o------o---- site
		# -----+------------+-
		if(start <= site[1] & end >= site[n]) return(seq_along(site))
		
		if(start <= site[1]) {
			i1 = 1
		} else {
			i1 = .Internal(findInterval(site, start, FALSE, FALSE))
		}

		if(end >= site[n]) {
			i2 = n
		} else {
			i2 = .Internal(findInterval(site, end, FALSE, FALSE))
		}
		
		# --------o------o---- site
		# ----------+--+------
		if(i1 == i2) {
			if(site[i1] < start && site[i2+1] > end) return(NULL)
		}
		
		if(start > site[1] && site[i1] < start) i1 = i1 + 1
		return(i1:i2)
	}
	
	id = integer(0)
	if(is.null(filter_fun)) {
		for(i in seq_along(start)) {
			ind = .extract(start[i], end[i])
			id = c(id, ind)
		}
	} else {
		for(i in seq_along(start)) {
			ind = .extract(start[i], end[i])
			if(filter_fun(site[ind])) {
				id = c(id, ind)
			}
		}
	}
	
	if(!index) {
		return(as.integer(site[id]))
	} else {
		return(id)
	}
}
