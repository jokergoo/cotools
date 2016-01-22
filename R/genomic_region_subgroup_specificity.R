
# == title
# Find common genomic regions across several samples
#
# == param
# -gr_list a list of `GenomicRanges::GRanges`, better have names
# -min_width minimum width of the common regions
# -min_coverage minimum cross-sample coverage for the common regions
#
# == details
# A common region is defined a region which at least i samples overlap with.
# The output can be sent to `subgroup_specific_genomic_regions` to find subgroup
# specific regions.
#
# == value
# a GRanges object contains coordinates of common regions. The columns in meta data
# are percent of the common region which is covered by regions in every sample.
common_regions = function(gr_list, min_width = 0, min_coverage = 1) {
	
	if(min_coverage < 1) {
		stop("`mincoverage` should be larger than 0.\n")
	}
	
	# merge all gr into one data frame
	gr_merge = gr_list[[1]]
	for(i in seq_along(gr_list)[-1]) {
		
		gr_merge = c(gr_merge, gr_list[[i]])

		message(qq("merging @{i}/@{length(gr_list)} samples"))
	}
	
	message("calculating cross-sample coverage")
	cov = coverage(gr_merge)
	gr_cov = as(cov, "GRanges")
	
	gr_common = reduce(gr_cov[mcols(gr_cov)$score >= min_coverage])

	if(length(gr_common) == 0) {
		message("No common regions was detected.\n")
		return(NULL)
	}
	
	gr_common = gr_common[ width(gr_common) >= min_width ]

	message(qq("there are @{length(gr_common)} common regions"))
	
	# calculate the percentage for each CR covered by gr_list
	message("overlapping `gr_list` to common regions.")
	gr_common = annotate_to_genomic_features(gr_common, gr_list, type = "percent")
	
	attr(gr_common, "generated_by") = "common_regions"

	return2(gr_common)
}

#### make some diagnosis plots

# == title
# find subgroup specific regions
#
# == param
# -gr a GRanges object coming from `common_regions`
# -factors classes which correspond to samples in ``gr``
# -present how to define a common region is specifically present in one subgroup?
#          The subgroup specificity is calculated based on the precent matrix stored in ``gr``.
#          For each subgroup which is defined in ``factors``,
#          if the value is a single numeric value, it means the mean value should be larger than this.
#          It can also be a function for which the input is the vector of percent in corresponding subgroup
#          and the output should be a single logical value.
# -absent how to define a common region is specifically absent in one subgroup. Format is same as ``present``
# -type It uses a string containing 1 and 0 to represent types of specificity. E.g. '1100' means
#       present in subgroup 1 and 2 while absent in subgroup 3 and 4. By default, it will output
#       all combination of subgroup specificities.
#
# == value
# In fact, this function split the original ``gr`` into a list in which each element
# contains regions corresponding to different subgroup specificity.
#
# The returned value can be sent to `plot_subgroup_specificity_heatmap` to visualize the specificity.
subgroup_specific_genomic_regions = function(gr, factors, 
	present = 0.7, absent = 0.3, type = NULL) {

	if(is.null(attr(gr, "generated_by"))) {
		stop("`gr` should come from `common_regions()`.\n")
	} else if(attr(gr, "generated_by") != "common_regions") {
		stop("`gr` should come from `common_regions()`.\n")
	}

	if(length(factors) != ncol(mcols(gr))) {
		stop("Length of `factors` should be same as number of samples in `gr`\n")
	}

	level = unique(factors)
	if(is.null(type)) {
		ss = expand.grid(rep(list(0:1), length(level)))
		type = apply(ss, 1, paste, collapse = "")
	} else {
		ss = as.matrix(as.data.frame(strsplit(type, "")))
		ss = t(ss)
		dm = dim(ss)
		ss = as.integer(ss)
		dim(ss) = dm
		
		if(ncol(ss) != length(level)) {
			stop("`nchar` of `type` should be same as the number of levels.\n")
		}
	}
	
	mat = mcols(gr)[, grepl("^overlap_to_", colnames(mcols(gr)))]
	mat = as.matrix(mat)
	
	gr_list = vector("list", length(type))
	names(gr_list) = type
	
	for(i in seq_along(type)) {
		message(qq("extracting pattern of '@{type[i]}'..."))
		l = sapply(seq_len(nrow(mat)), function(k) {
			x = mat[k, ]
			for(j in seq_along(level)) {
				x2 = x[factors == level[j]]
				if(ss[i, j]) {  # present
					if(is.function(present)) {
						l2 = present(x2)
					} else {
						l2 = mean(x2) >= present
					}
				} else {  # absent
					if(is.function(absent)) {
						l2 = absent(x2)
					} else {
						l2 = mean(x2) <= absent
					}
				}
				
				if(!l2) {
					return(FALSE)
				}
			}
			return(TRUE)
		})
		
		gr_list[[i]] = gr[l]
	}
	
	attr(gr_list, "factors") = factors
	attr(gr_list, "generated_by") = "subgroup_specific_genomic_regions"

	return2(gr_list)
}

# == title
# heatmap for visualizing subgroup specific genomic regions 
#
# == param
# -gr_list a list of GRanges which is coming from `subgroup_specific_genomic_regions`
# -genomic_features Genomic features that are used for annotatate to regions in ``gr_list``, it should be a list of GRanges objects
#
# == details
# columns are clustered inside each subgroup.
plot_subgroup_specificity_heatmap = function(gr_list, genomic_features = NULL) {
	
	if(is.null(attr(gr_list, "generated_by"))) {
		stop("`gr` should come from `subgroup_specific_genomic_regions()`.\n")
	} else if(attr(gr_list, "generated_by") != "subgroup_specific_genomic_regions") {
		stop("`gr` should come from `subgroup_specific_genomic_regions()`.\n")
	}

	factors = attr(gr_list, "factors")
	
	# column name
	cn = gsub("^overlap_to_", "", grep("^overlap_to_", colnames(mcols(gr_list[[1]])), value = TRUE))

	annotation = data.frame(subgroup = factors)
	rownames(annotation) = cn
	
	level = unique(factors)

	mat = NULL
	gr_combine = GRanges()
	type = NULL

	gr_list_name = names(gr_list)
	for(i in seq_along(gr_list)) {
		if(length(gr_list[[i]]) < 2) next
		# cluster rows for each sub-matrix
		sub_matrix = mcols(gr_list[[i]])[, grepl("^overlap_to_", colnames(mcols(gr_list[[i]])))]
		sub_matrix = as.matrix(sub_matrix)
		colnames(sub_matrix) = cn
		message(qq("cluster sub-matrix @{gr_list_name[i]} on rows (@{length(gr_list[[i]])} rows)."))
		rclust = hclust(dist(sub_matrix))
		mat = rbind(mat, sub_matrix[rclust$order, ])
		gr_combine = c(gr_combine, gr_list[[i]][rclust$order])
		type = c(type, rep(gr_list_name[i], length(gr_list[[i]])))
	}
	mcols(gr_combine) = NULL

	# globally cluster inside each subgroup
	mat2 = NULL
	for(le in level) {
		message(qq("cluster @{le} on columns."))
		m = mat[, factors == le]
		cclust = hclust(dist(t(m)))
		mat2 = cbind(mat2, m[, cclust$order])
	}
	
	ht_list = Heatmap(mat2, name = "pct to CR", col = colorRamp2(c(0, 1), c("white", "blue")), 
		top_annotation = HeatmapAnnotation(df = annotation), top_annotation_height = unit(5, "mm"),
		cluster_rows = FALSE, cluster_columns = FALSE, split = type)
	
	w = width(gr_combine)
	
	ht_list = ht_list + 
		Heatmap(type, name = "type", width = unit(4, "mm")) + 
		Heatmap(w, name = "width", col = colorRamp2(c(min(w), quantile(w, 0.9)), c("#EEEEEE", "#000000")), width = unit(8, "mm"),
			cluster_columns = FALSE)

	if(!is.null(genomic_features)) {
		if(inherits(genomic_features, "GRanges")) {
			genomic_features = list("genomic_feature" = genomic_features)
		}
		gr_combine = annotate_to_genomic_features(gr_combine, genomic_features, prefix = "")
		
		ht_list = ht_list + Heatmap(as.matrix(mcols(gr_combine)), name = "anno", col = colorRamp2(c(0, 1), c("white", "orange")),
			width = unit(2, "mm")*length(genomic_features))
	}

	message(qq("generating heatmap, (@{nrow(mat2)}, @{ncol(mat2)})"))
	ht_list
}
