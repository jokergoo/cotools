
cr_enriched_at_tss = function(cr, txdb) {
	
	gene = genes(txdb)
	tx = transcripts(txdb)
	tx_list = transcriptsBy(txdb, "gene")

	gene = gene[names(gene) %in% cr$gene_id]
	tx = tx[tx$tx_name %in% cr$nearest_tx_tss]
	names(tx) = tx$tx_name
	tx$gene_id = structure(cr$gene_id, names = cr$nearest_tx_tss)[names(tx)]

	gene_tss = promoters(gene, upstream = 1, downstream = 0)
	tx_tss = promoters(tx, upstream = 1, downstream = 0)

	# align to gene tss
	mat_pos = normalizeToMatrix(cr[cr$corr > 0], gene_tss, mapping_column = "gene_id",
        extend = c(10000, 20000), mean_mode = "absolute", w = 50)
	mat_neg = normalizeToMatrix(cr[cr$corr < 0], gene_tss, mapping_column = "gene_id",
        extend = c(10000, 20000), mean_mode = "absolute", w = 50)
	x_pos = colSums(mat_pos)
	x_neg = colSums(mat_neg)


	tx_list = tx_list[names(gene)]
	tx_list = lapply(names(tx_list), function(gi) {
		gr = tx_list[[gi]]
		gr = gr[gr$tx_name != gi]
		gr$gene_id = rep(gi, length(gr))
		gr
	})

	tx2 = do.call("c", tx_list)

	mat_tx = normalizeToMatrix(tx2, gene_tss, mapping_column = "gene_id",
		extend = c(10000, 20000), mean_mode = "w0", w = 50)
	mat_nearest_tx = normalizeToMatrix(tx, gene_tss, mapping_column = "gene_id",
		extend = c(10000, 20000), mean_mode = "w0", w = 50)
	
	par(mar = c(4, 4, 4, 4))
	plot(x_pos, ylim = c(0, max(c(x_pos, x_neg))), type = "l", col = "red", axes = FALSE, xlab = 'pos relateive to gene tss', ylab = "coverage")
	lines(x_neg, col = "green")
	abline(v = length(x_pos)/3, lty = 2, col = "grey")
	axis(side = 2)
	axis(side = 1, at = seq(1, length(x_pos), length = 7), label = c(-10000, -5000, 0, 5000, 10000, 15000, 20000))

	par(new = TRUE)
	plot(colSums(mat_tx), type = "l", col = "blue", axes = FALSE, ann = FALSE)
	lines(colSums(mat_nearest_tx), col = "orange")
	axis(side = 4)

	legend("topright", lty = 1, col = c("green", "red", "blue"), legend = c("neg_cr", "pos_cr", "tx"))
	par(new = FALSE)

	# align to tx tss
	mat_pos = normalizeToMatrix(cr[cr$corr > 0], tx_tss, mapping_column = "nearest_tx_tss",
        extend = c(10000, 20000), mean_mode = "absolute", w = 50)
	mat_neg = normalizeToMatrix(cr[cr$corr < 0], tx_tss, mapping_column = "nearest_tx_tss",
        extend = c(10000, 20000), mean_mode = "absolute", w = 50)
	x_pos = colSums(mat_pos)
	x_neg = colSums(mat_neg)

	par(mar = c(4, 4, 4, 4))
	plot(x_pos, ylim = c(0, max(c(x_pos, x_neg))), type = "l", col = "red", axes = FALSE, xlab = 'pos relateive to tx tss', ylab = "coverage")
	lines(x_neg, col = "green")
	abline(v = length(x_pos)/3, lty = 2, col = "grey")
	axis(side = 2)
	axis(side = 1, at = seq(1, length(x_pos), length = 7), label = c(-10000, -5000, 0, 5000, 10000, 15000, 20000))

	legend("topright", lty = 1, col = c("green", "red"), legend = c("neg_cr", "pos_cr"))
}

enriched_heatmap_list_on_gene = function(cr, cgi, txdb, expr, hm_list = NULL, hm_name = NULL, on = "tss", by = "gene", ...) {

	sample_id = attr(cr, "sample_id")
	factor = attr(cr, "factor")

	if(on == "tss") {
		if(by == "gene") {
			qqcat("extracting gene tss\n")
			gene = genes(txdb)
			tss = promoters(gene, upstream = 1, downstream = 0)
			tss = tss[names(tss) %in% cr$gene_id]
			mapping_column = "gene_id"
		} else {
			qqcat("extracting nearest tx tss\n")
			tx = transcripts(txdb)
			nearest_tx = tx[tx$tx_name %in% cr_filtered$nearest_tx_tss]
			names(nearest_tx) = nearest_tx$tx_name
			nearest_tx$gene_id = structure(cr_filtered$gene_id, names = cr_filtered$nearest_tx_tss)[names(nearest_tx)]
			nearest_tx_tss = promoters(nearest_tx, upstream = 1, downstream = 0)
			tss = nearest_tx_tss
			mapping_column = "nearest_tx_tss"
		}

		target = tss
		target_ratio = 0.1
		axis_name = c("-5KB", "TSS", "5KB")
	} else {
		qqcat("extracting gene body\n")
		gene = genes(txdb)
		gene = gene[names(gene) %in% cr$gene_id]

		target = gene
		target_ratio = 0.6
		mapping_column = "gene_id"
		axis_name = c("-5KB", "TSS", "TES", "5KB")
	}

	target = sort(target)

	mat = normalizeToMatrix(cr, target, mapping_column = mapping_column,
	        extend = 5000, mean_mode = "absolute", w = 50, target_ratio = target_ratio)

	l = rowSums(mat) > 0
	mat = mat[l, ]
	target = target[l]

	mat_cgi = normalizeToMatrix(cgi, target,
	        extend = 5000, mean_mode = "absolute", w = 50, target_ratio = target_ratio)

	if(!missing(hm_list)) {
		lt = enrich_with_histone_mark(target, hm_list, sample_id, factor, return_arr = TRUE, target_ratio = target_ratio)

		arr = lt[[1]]
		# detect regions that histone marks correlate to expression
		expr2 = expr[names(target), sample]
		cor_mat = matrix(nr = nrow(expr2), nc = ncol(mat))
		cor_p_mat = cor_mat

		counter = set_counter(nrow(cor_mat))
		for(i in seq_len(nrow(cor_mat))) {
			counter()
		    for(j in seq_len(ncol(cor_mat))) {
		        x = cor(arr[i, j, ], expr2[i, ], method = "spearman")
		        cor_mat[i, j] = x
		        cor_p_mat[i, j] = cor.test(arr[i, j, ], expr2[i, ], method = "spearman")$p.value
		    }
		}
		cat("\n")
		cor_mat[is.na(cor_mat)] = 0
		cor_p_mat = p.adjust(cor_p_mat, method = "BH")
		l1 = cor_p_mat < 0.01 & cor_mat > 0
		cor_mat[l1] = 1
		l2 = cor_p_mat < 0.01 & cor_mat < 0
		cor_mat[l2] = -1 
		cor_mat[!(l1 | l2)] = 0
		cor_mat = copyAttr(mat, cor_mat)
	}


	hist_mat_list = lt[[2]]

	meth_mat_list = enrich_with_methylation(target, sample_id, factor, target_ratio = target_ratio)

	base_expr = rowMeans(expr[names(target), , drop = FALSE])

	expr2 = t(apply(expr[names(target), , drop = FALSE], 1, function(x) {
			q = quantile(x, 0.9)
			x[x > q] = q
			scale(x)
		}))
	rel_expr = tapply(seq_along(factor), factor, function(ind) rowMeans(expr2[, ind]))
	rel_expr = do.call("cbind", rel_expr)

	if(all(cr$corr > 0)) {
		cr_col = "red"
		cr_name = "pos_cr"
	} else if(all(cr$corr < 0)) {
		cr_col = "darkgreen"
		cr_name = "neg_cr"
	} else {
		cr_col = "blue"
		cr_name = "cr"
	}

	qqcat("making heatmap...\n")
	ht_list = EnrichedHeatmap(mat, col = c("white", cr_col), name = cr_name, cluster_rows = TRUE, show_row_hclust = FALSE,
	              top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = cr_col))), 
	              top_annotation_height = unit(2, "cm"), column_title = "meth <cor> expr", axis_name = axis_name)
	if(!is.null(hm_list)) {
	    ht_list = ht_list + EnrichedHeatmap(cor_mat, col = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red")), name = qq("@{hm_name}"), 
	          	column_title = qq("@{hm_name} <cor> expr"))
	}

	ht_list = ht_list +
	          EnrichedHeatmap(mat_cgi, col = c("white", "darkorange"), name = "CGI",
	          	  top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "darkorange"))), 
	              top_annotation_height = unit(2, "cm"), column_title = "CGI", axis_name = axis_name) +
	          Heatmap(width(gene[names(target)]), name = "gene_len", col = colorRamp2(c(2000, 10000), c("black", "white")), 
	          	  width = unit(5, "mm")) +
	          Heatmap(rel_expr, name = "rel_expr", width = unit(12, "mm"), show_row_names = FALSE, cluster_columns = FALSE) +
	          Heatmap(base_expr, name = "base_e", width = unit(5, "mm"), show_row_names = FALSE)

	if(!null(hm_list)) {
		for(type in names(hist_mat_list)) {
			ht_list = ht_list + EnrichedHeatmap(hist_mat_list[[type]], col = c("white", "red"), name = qq("hist_@{type}"),
				column_title = qq("@{hm_name}_@{type}"), axis_name = axis_name, show_heatmap_legend = FALSE)
		}
	}
	for(type in names(meth_mat_list)) {
		ht_list = ht_list + EnrichedHeatmap(meth_mat_list[[type]], col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
			name = qq("meth_@{type}"), column_title = qq("meth_@{type}"), axis_name = axis_name, show_heatmap_legend = FALSE)
	}

	draw(ht_list, heatmap_legend_side = "bottom", ...)

}


enriched_heatmap_list_on_tss_cgi = function(cr, cgi, txdb, expr, hm_list = NULL, hm_name = NULL, by = "gene", ...) {

	sample_id = attr(cr, "sample_id")
	factor = attr(cr, "factor")

	gene = genes(txdb)

	strand(cr) = strand(gene[cr$gene_id])

	tss = promoters(gene, upstream = 1, downstream = 0)
	tss = tss[names(tss) %in% cr$gene_id]   # ???

	dist = distanceToNearest(cgi, tss)
	l = dist@elementMetadata@listData$distance < 5000
	mtch = as.matrix(dist)[l, ]

	cgi = cgi[mtch[,1]]
	cgi$nearest_tss = names(tss[mtch[,2]])

	cgi_extend = cgi; start(cgi_extend) = start(cgi) - 5000; end(cgi_extend) = end(cgi) + 5000
	mtch = as.matrix(findOverlaps(cr, cgi_extend))
	cgi2 = cgi[unique(mtch[, 2])]
	cgi_extend = cgi_extend[unique(mtch[, 2])]

	mat3 = normalizeToMatrix(cr[unique(mtch[, 1])], cgi2, extend = 5000, w = 50)
	l = rowSums(mat3) > 0
	mat3 = mat3[l, ]
	cgi2 = cgi2[l]

	n_tss = countOverlaps(cgi_extend, tss)

	dist = distanceToNearest(cgi2, tss)

	strd = as.vector(strand(gene[cgi2$nearest_tss]))

	if(all(cr$corr) > 0) {
		cr_col = "red"
		cr_name = "pos_cr"
	} else if(all(cr$corr) < 0) {
		cr_col = "green"
		cr_name = "neg_cr"
	} else {
		cr_col = "blue"
		cr_name = "cr"
	}

	ht_list = EnrichedHeatmap(mat3, col = c("white", cr_col), cluster_rows = TRUE, show_row_hclust = FALSE, split = strd,
	                top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = c("darkgreen", "red")))), 
	                top_annotation_height = unit(2, "cm"), column_title = qq("@{cr_name} ~ cgi"), axis_name = c("-5kb", "start", "end", "5kb"))+
	            Heatmap(width(cgi2[l]), name = "cgi_width", col = colorRamp2(c(0, 1000), c("orange", "white")),
	                width = unit(5, "mm")) +
	            Heatmap(strd, name = "strand", col = c("+" = "red", "-" = "darkgreen"),
	                width = unit(5, "mm")) +
	            Heamtap(n_tss, name = "nearby_tss", col = c("white", "black"), width = unit(5, "mm"))


	if(!is.null(hm_list)) {
		hist_mat_list = enrich_with_histone_mark(cgi2, hm_list, sample_id, factor)

		for(type in names(hist_mat_list)) {
			ht_list = ht_list + EnrichedHeatmap(hist_mat_list[[type]], col = c("white", "red"), name = qq("hist_@{type}"),
				column_title = type, axis_name = c("-5kb", "start", "end", "5kb"), show_heatmap_legend = FALSE)
		}
	}

	meth_mat_list = enrich_with_methylation(cgi2, sample_id, factor)
	for(type in names(meth_mat_list)) {
		ht_list = ht_list + EnrichedHeatmap(meth_mat_list[[type]], col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
			name = qq("meth_@{type}"), column_title = type, axis_name = c("-5kb", "start", "end", "5kb"), show_heatmap_legend = FALSE)
	}
	draw(ht_list, heatmap_legend_side = "bottom", ...)
}

enrich_with_histone_mark = function(target, hm_list, sample_id, factor, target_ratio = 0.1, return_arr = FALSE) {

	target_name = deparse(substitute(target))
	# 200 is number of windows (5000 + 5000)/50

	if(all(width(target) <= 1)) {
		w = 200
	} else {
		w = ceiling(200/(1-target_ratio)*target_ratio)
	}

	arr = array(dim = c(length(target), w, length(hm_list)))
	sample = names(hm_list)
	for(i in seq_along(sample)) {
		qqcat("@{sample[i]}: normalize histone modifications to @{target_name}.\n")
	    tm = normalizeToMatrix(hm_list[[i]], target, value_column = "density", extend = 5000, mean_mode = "w0", w = 50, target_ratio = target_ratio)
	    arr[, , i] = tm
	}

	########### matrix for each subgroup
	hist_mat_list = list()
	for(type in unique(factor)) {
	    l = sample %in% sample_id[factor]
	    hist_mat_list[[type]] = apply(arr[, , l], c(1, 2), mean)
	    hist_mat_list[[type]] = copyAttr(tm, hist_mat_list[[type]])
	}

	if(return_arr) {
		return(list(arr = arr, list = hist_mat_list))
	} else {
		return(hist_mat_list)
	}
}

enrich_with_methylation = function(target, sample_id, factor, target_ratio = 0.1) {

	# extrace methylation value which in [-5k, 5k] from TSS, and calculate mean methylation in each subgroup
	target_extend = GRanges(seqnames = seqnames(target), ranges = IRanges(start(target)-5000, end(target)+5000))
	# process raw methylation data
	meth_gr = GRanges()
	for(chr in sort(unique(seqnames(target)))) {
	    cat(chr, "\n")
	    methylation_hooks$set(chr)
	    gr = methylation_hooks$GRanges()
	    mtch = as.matrix(findOverlaps(gr, target_extend))
	    ind = unique(mtch[, 1])
	    mm = do.call("cbind", lapply(unique(SAMPLE$type), function(s) rowMeans(methylation_hooks$meth(row_index = ind, col_index = SAMPLE$type == s))))
	    gr = gr[ind]
	    mcols(gr) = mm
	    meth_gr = c(meth_gr, gr)
	}
	mm = mcols(meth_gr)
	colnames(mm) = unique(factor)
	mcols(meth_gr) = mm
	meth_mat_list = list()
	for(type in unique(factor)) {
		qqcat("@{type}: normalize methylation to target.\n")
	    meth_mat_list[[type]] = normalizeToMatrix(meth_gr, target, value_column = type, extend = 5000, w = 50, 
	    	mean_mode = "absolute", empty_value = 0.5, target_ratio = target_ratio)
	}

	return(meth_mat_list)
}

