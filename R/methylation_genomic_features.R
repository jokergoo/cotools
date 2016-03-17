### function for methylation analysis


# == title
# heatmap for mean methylation in genomic features and annotation to other genomic features
#
# == param
# -gr object returned by `get_mean_methylation_in_genomic_features`
# -annotation subtype of samples
# -annotation_color colors of subtypes
# -txdb A `GenomicFeatures::TcDb` object
# -gf_list a list of genomic features which are used as row annotations
# -gf_type how to overlap genomic features
# -min_mean_range minimal range between mean value in subtypes
# -cutoff if subtype information is provided, p-value for the oneway ANOVA test
# -adj_method how to calculate adjusted p-values
# -title title of the plot
# -cluster_cols how to cluster columns
# -ha heatmap annotations
# -... pass to `ComplexHeatmap::Heatmap`
# 
heatmap_diff_methylation_in_genomic_features = function(gr, annotation, 
	annotation_color = structure(seq_along(unique(annotation)), names = unique(annotation)), 
	txdb = NULL, gf_list = NULL, gf_type = "percent", 
	min_mean_range = 0.2, cutoff = 0.05, adj_method = "BH", title = NULL, 
	cluster_cols = c("subgroup", "all", "none"), ha = NULL, ...) {
	
	# if(length(annotation) != ncol(mcols(gr))) {
	# 	stop("Length of `annotation` should be equal to ncol of `mcols(gr)`\n")
	# }
	
	cluster_cols = match.arg(cluster_cols)[1]

	nn = length(gr)
	
	mat = as.matrix(mcols(gr))
	mat = mat[, colnames(mat) != "ncpg", drop = FALSE]
	mcols(gr) = mat
	nr0 = nrow(mat)
	qqcat("@{nr0} rows in `gr`\n")
	
	if(min_mean_range > 0) {

		x = tapply(seq_len(ncol(mat)), annotation, function(ind) {
				rowMeans(mat[, ind, drop = FALSE])
			})
		dim(x) = NULL
		x = as.matrix(data.frame(x))
		l = apply(x, 1, diameter) > min_mean_range

		qqcat("@{sum(l)}/@{length(l)} rows in `gr` after filtering by `min_mean_range`\n")
		gr = gr[l]
		mat = as.matrix(mcols(gr))
	}
	
	if(cutoff < 1) {
		p = numeric(nrow(mat))
		for(j in seq_len(nrow(mat))) {
			group = factor(annotation)
			data = mat[j, ]
			data = data + rnorm(length(data), 0, 0.01)
			df = data.frame(data = data, group = group)
			p[j] = oneway.test(data ~ group, data = df)$p.value
		}
		p = p.adjust(p, method = adj_method)
		l = p < cutoff

		gr = gr[l]
		mat = as.matrix(mcols(gr))
			
		qqcat("@{nrow(mat)} rows in `gr` after filtering by oneway-ANOVA test (FDR < @{cutoff}).\n")
	}
	
	ogr = gr
	# don't need that much rows
	nr = nrow(mat)
	if(nr > 5000) {
		l = sort(sample(nr, 5000))
		gr = gr[l]
		mat = as.matrix(mcols(gr))
	}

	### draw heatmap ###

	### cluster columns ###########
	if(cluster_cols == "subgroup") {
		type = unique(annotation)
		mm = NULL
		for(t in type) {
			m1 = mat[, annotation == t, drop = FALSE]
			hc = hclust(dist(t(m1)))
			m1 = m1[, hc$order, drop = FALSE]
			cn = c(colnames(mm), colnames(m1))
			mm = cbind(mm, m1)
			colnames(mm) = cn
		}
		mat = mm
		cluster_cols = FALSE
	} else if(cluster_cols == "all") {
		cluster_cols = TRUE
	} else {
		cluster_cols = FALSE
	}

	if(is.null(ha))	ha = HeatmapAnnotation(df = data.frame(anno = annotation), col = list(anno = annotation_color))
	ht_list = Heatmap(mat, col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), top_annotation = ha,
		cluster_columns = cluster_cols, show_row_dend = FALSE,
		heatmap_legend_param = list(title = "methylation"), ...)
	
	# make a copy of `gr`
	gr2 = gr
	mcols(gr2) = NULL

	## width of gr
	w = width(gr2)
	if(length(unique(w)) != 1) {
		ht_list = ht_list + Heatmap(w, name = "width", col = colorRamp2(quantile(w, c(0, 0.8)), c("black", "white")),
			heatmap_legend_param = list(title = "width"), width = unit(5, "mm"))
	}

	# calcualte distance to closest tss
	if(!is.null(txdb)) {
		gene = genes(txdb)
		tss = promoters(gene, upstream = 0, downstream = 1)
		dst = nearest(gr2, tss, select = "arbitrary")
		ht_list = ht_list + Heatmap(dst, name = "tss_dist", col = colorRamp2(c(0, 10000), c("black", "white")),
			heatmap_legend_param = list(title = "tss_dist"), width = unit(5, "mm"))
	}

	# annotate to other gf
	if(!is.null(gf_list)) {
		for(i in seq_along(gf_list)) {
			if(!inherits(gf_list[[i]], "GRanges")) {
				gf_list[[i]] = GRanges(seqnames = gf_list[[i]][[1]],
					                   ranges = IRanges(gf_list[[i]][[2]], gf_list[[i]][[3]]))
			}
		}

		gf_name = names(gf_list)
		for(gn in gf_name) {
			gr2 = annotate_to_genomic_features(gr2, gf_list[[gn]], gn, prefix = "", type = gf_type)
		}

		if(gf_type == "percent") {
			col = colorRamp2(c(0, 1), c("white", "orange"))
		} else {
			col = colorRamp2(c(0, 0.5), c("white", "orange"))
		}
		ht_list = ht_list + Heatmap(as.data.frame(mcols(gr2)), col = col, cluster_columns = FALSE,
			heatmap_legend_param = list(title = "anno_gf"), width = unit(5*length(gf_list), "mm"))
	}
	
	draw(ht_list, column_title = qq("@{title}\n@{length(ogr)}/@{nn} rows"))
	
	return(invisible(ogr))
}


# == title
# calculate average methylation value in a list of regions
#
# == param
# -sample_id  a list of sample IDs
# -gf_list a list of genomic features in `GenomicRanges::GRanges` class
# -average    whether to calcualte average methylation in a interval? if not,
#             the function will randomly sample some CpG sites in the interval
# -p if average is FALSE, the probability to pick cpg sites
# -chromosome chromosome name
# -filter_fun filtering function on sites in each intersection
#
# == value
# a list of ``GRanges`` objects.
get_mean_methylation_in_genomic_features = function(sample_id, gf_list, average = TRUE, p = 0.001,
	chromosome = paste0("chr", 1:22),
	filter_fun = function(s) length(s) >= 10 && any(diff(s) < 50)) {
	
	# initialize the mat_list
	# mat_list has the same name as `gf_list`
	
	gr_list = rep(list(GRanges(seqnames = character(0), ranges = IRanges(start = integer(0), end = integer(0)))), 
		          length(gf_list))
	names(gr_list) = names(gf_list)

	for(chr in chromosome) {
		methylation_hooks$set(chr)
		meth_mat = methylation_hooks$meth(col_index = sample_id)
		meth_gr = methylation_hooks$GRanges()
		meth_site = methylation_hooks$site()

		for(i in seq_along(gf_list)) {
			qqcat("overlapping to @{names(gf_list)[i]} on @{chr}\n")
			mtch = as.matrix(findOverlaps(gf_list[[i]], meth_gr))  ## <- test memory usage

			if(average){
				l = tapply(mtch[,2], mtch[,1], function(i) {
						x = meth_site[i]
						filter_fun(x)
					})
				mean_meth = tapply(mtch[,2], mtch[,1], function(i) colMeans(meth_mat[i, , drop = FALSE]))
				mean_meth = mean_meth[l]

				ncpg = tapply(mtch[,2], mtch[,1], length)

				mean_meth_mat = matrix(unlist(mean_meth), nrow = length(mean_meth), byrow = TRUE)
				rownames(mean_meth_mat) = names(mean_meth); colnames(mean_meth_mat) = sample_id
				ind = as.integer(names(mean_meth))
				gr = gf_list[[i]][ind]
				mcols(gr) = cbind(as.data.frame(mean_meth_mat), ncpg = ncpg[names(mean_meth)])
			} else {
				l = unique(mtch[, 2])
				gr = meth_gr[l]
				mcols(gr) = as.data.frame(meth_mat[l, , drop = FALSE])
				if(length(gr) > 10000) gr = gr[sample(c(TRUE, FALSE), length(gr), prob = c(p, 1-p), replace = TRUE)]
				gc(verbose = FALSE)
			}
			gr_list[[i]] = c(gr_list[[i]], gr)
		}
	}

	return2(gr_list)
}
