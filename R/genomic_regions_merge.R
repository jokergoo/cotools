
reduce_regions = function(gr, max.gap = 500, ratio = 1) {
	if(length(gr) == 0) {
		return(GRanges())
	}
	all_seqnames = unique(as.vector(seqnames(gr)))
	df = NULL
	for(sn in all_seqnames) {
		gr_subset = sort(gr[seqnames(gr) == sn])
		d = dynamic_region_finder(start(gr_subset), end(gr_subset), max.gap = max.gap, ratio = ratio)
		df = rbind(df, cbind(chr = rep(sn, nrow(d)), d))
	}
	GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]))
}

dynamic_region_finder = function(start, end, max.gap = 500, ratio = 1, value_df = NULL) {

	 nr = length(start)

	 if(!is.null(value_df) && nrow(value_df) != nr) {
	 	stop("Number of rows of `value_df` should be equal to `start`.")
	 }
	 
	 i_d = 0
	 while(1) {
	 	i_d = i_d + 1
	 	
	 	new_region = merge_region(start, end, max.gap = max.gap, ratio = ratio, value_df = value_df)
	 	qqcat("@{i_d} iterations for dynamic merge, now regions: @{dim(new_region)[1]}...\n")

	 	# no changes
	 	if(nr == nrow(new_region)) {
	 		break
	 	}

	 	nr = nrow(new_region)
	 	start = new_region$start
	 	end = new_region$end
	 	if(!is.null(value_df)) {
	 		value_df = new_region[, 3:ncol(new_region), drop = FALSE]
	 	}
	 }

	 return(new_region)
}

# == title
# merge two close regions
# 
# == param
# -start    position of start site
# -end      position of end site
# -max.gap  maximum gap
# -ratio    distance = ratio*width
# -value_df other value coresponding to each region
#
# == details
# If regions are merged, corresponding value in `value_df` will also be merged
merge_region = function(start, end, max.gap = 500, ratio = 1, value_df = NULL) {
	
	n_region = length(start)

	if(!is.null(value_df) && nrow(value_df) != n_region) {
	 	stop("Number of rows of `value_df` should be equal to `start`.")
	 }

	 # no need to merge
	if(length(start) <= 1) {
		if(is.null(value_df)) {
			return(data.frame(start = start, end = end))
		} else {
			return(data.frame(start = start, end = end, value_df))
		}
	}

	# <<== put merge part as cpp code!!
	width = end - start + 1
	if(class(ratio) == "bp") {
		dist = ratio
	} else {
		dist = ratio*width
	}
	dist = ifelse(dist > max.gap, max.gap, dist)
	extended_start = start - dist
	extended_end = end + dist
	extended_start[extended_start < min(start)] = min(start)
	extended_end[extended_end > max(end)] = max(end)

	# `l` means whether current region should be merged with next one
	l = rep(FALSE, length = n_region - 1)
	for(i in seq_len(n_region)) {
		if(i %% 10000 == 0) {
			qqcat("[merging] @{i}/@{n_region} regions checked\n")
		}

		if(i == 1) {

			# backward looking, which `start` is in [extended_start, extended_end]
			for(j in seq(i+1, n_region)) {
				if(start[j]-1 <= extended_end[i]) {  # overlap or connected, start[j]: start of the second row, extended_end[i]: end of the first row
					l[(i+1):j - 1] = l[(i+1):j - 1] | TRUE
				} else {
					break
				}
			}

		} else if(i == n_region) {

			# forward looking, which `end` is in [extended_start, extended_end]
			for(j in seq(i-1, 1)) {
				if(end[j]+1 >= extended_start[i]) {  # overlap
					l[j:(i-1)] = l[j:(i-1)] | TRUE
				} else {
					break
				}
			}

		} else {

			# forward looking, which `end` is in [extended_start, extended_end]
			for(j in seq(i-1, 1)) {
				if(end[j]+1 >= extended_start[i]) {  # overlap or connected
					l[j:(i-1)] = l[j:(i-1)] | TRUE
				} else {
					break
				}
			}

			# backward looking, which `start` is in [extended_start, extended_end]
			for(j in seq(i+1, n_region)) {
				if(start[j]-1 <= extended_end[i]) {  # overlap or connected
					l[(i+1):j - 1] = l[(i+1):j - 1] | TRUE
				} else {
					break
				}
			}

		}
		
	}

	# merge regions based

	
	#width = end - start + 1
	#shorter_width = pmin(width[-n_region], width[-1])
	#distance = start[-1] - end[-n_region]

	# merge
	# merge_index == 1 means merge region 1 and 2
	# merge_index == 2 means merge region 2 and 3
	#merge_index = logical_segment(distance <= max.gap & distance <= ratio*shorter_width)
	merge_index = logical_segment(l)
	# do not need to merge
	if(nrow(merge_index) == 0) {
		if(is.null(value_df)) {
			return(data.frame(start = start, end = end))
		} else {
			rownames(value_df) = NULL
			return(data.frame(start = start, end = end, value_df))
		}
	} 
	
	mat = t(apply(merge_index, 1, function(x) {
			if(is.null(value_df)) {
				return(c(start[ x[1] ], end[ x[2] + 1 ]))
			} else {
				v = c(start[ x[1] ], end[ x[2] + 1 ])
				v2 = sapply(value_df, function(y) sum(y[ (x[1]):(x[2]+1) ]))
				return(c(v, v2))
			}
		}))

	ind = NULL
	for(i in seq_len(nrow(merge_index))) {
		ind = c(ind, (merge_index[i,1]):(merge_index[i,2]+1))
	}
	
	complement_ind = setdiff(seq_along(start), ind)   # this is not an efficient way!  <=======

	if(length(complement_ind)) {
		if(is.null(value_df)) {
			mat = rbind(mat, cbind(start[complement_ind], end[complement_ind]))
		} else {
			mat = rbind(mat, cbind(start[complement_ind], end[complement_ind], as.matrix(value_df[complement_ind, ])))
		}
	}
	if(is.null(value_df)) {
		colnames(mat) = c("start", "end")
	} else {
		colnames(mat) = c("start", "end", names(value_df))
	}
	rownames(mat) = NULL
	mat = mat[order(mat[, 1]), ,drop = FALSE]
	return(as.data.frame(mat))
}

# == title
# segment by continuous logical values
#
# == param
# -l logical vector
#
# == details
# the logical vector will be segmented according to their values.
# It returns intervals for continuous `TRUE` values
#
# == values
# a data frame in which the first column is the index of start sites
# the second column is the index of end sites
logical_segment = function(l) {
	 w = which(l)

	 # no cluster
	 if(! length(w)) {
	 	return(data.frame(start_index = numeric(0), end_index = numeric(0)))
	 }

	 d = diff(w) > 1

	 # only one cluster
	 if(! any(d)) {
	 	return(data.frame(start_index = w[1], end_index = w[length(w)]))
	 } 

	 # more than one clusters
	 # second step: dynamic merge
	 #   if distance between two clusters are smaller than short cluster, merge the two clusters
	 di = which(d)
	 end = w[di]
	 start = w[di+1]
	 start = c(w[1], start)
	 end = c(end, w[length(w)])

	return(data.frame(start_index = start, end_index = end))
}

