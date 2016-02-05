
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


bp = function(x) {
	x = as.integer(x)
	class(x) = c(class(x), "bp")
	x
}

kb = function(x) {
	bp(x*1000)
}

mb = function(x) {
	bp(x*1000000)
}

print.bp = function(x) {
	x = paste0(x, "bp")
	print(x, quote = FALSE)
}


# value will be sum up
reduce2 = function(gr, max_gap = 1000, gap = bp(1000)) {
	if(length(gr) %in% c(0, 1)) {
		return(gr)
	}

	if(inherits(gap, "bp")) {
		gr_reduced = reduce(gr, min.gapwidth = gap, with.revmap = TRUE)
		if(is.null(mcols(gr))) {
			mcols(gr_reduced) = NULL	
		} else {
			mat = as.matrix(mcols(gr))
			cn = colnames(mat)
			mat = do.call("rbind", lapply(mcols(gr_reduced)[, 1], function(ind) colSums(mat[ind, , drop = FALSE])))
			colnames(mat) = cn
			mcols(gr_reduced) = mat
		}
		return(gr_reduced)
	} else {
		n = length(gr)
		gr_extend = gr
		gr_width = width(gr)
		start(gr_extend) = start(gr) - round(gr_width*gap)
		end(gr_extend) = end(gr) + round(gr_width*gap)
		gr_reduced = reduce2(gr_extend, max_gap = max_gap, gap = bp(max_gap))
		n2 = length(gr_reduced)

		i = 1
		if(n2 == 1) {
			return(gr_reduced)
		} else if(n == n2) {
			return(gr_reduced)
		} else {
			qqcat("  reduce interation @{i <- i + 1}..\n")
			gr = reduce2(gr_reduced, max_gap = max_gap, gap = gap)
		}
	}
}