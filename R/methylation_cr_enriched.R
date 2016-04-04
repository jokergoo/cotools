
# == title
# enrichment of cr on tss
#
# == param
# -cr cr
# -txdb txdb
#
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
	# lines(colSums(mat_nearest_tx), col = "orange")
	axis(side = 4)
	mtext("coverage for tx", side = 4, line = 3)

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

# == title
# Enriched heatmap
# 
# == param
# -cr cr
# -cgi CpG Island
# -txdb txdb
# -expr expr
# -hm_list a list of histome marks
# -hm_name names for hm
# -on tss or body
# -by gene or tx
# -hm_cor_p_cutoff cutoff for the correlation between hm and expression
# -show_expr whether show heatmap of expression
# -... pass to
#
enriched_heatmap_list_on_gene = function(cr, cgi, txdb, expr, hm_list = NULL, hm_name = NULL, on = "tss", by = "gene", 
	hm_cor_p_cutoff = 0.05, show_expr = TRUE, ...) {

	sample_id = attr(cr, "sample_id")
	factor = attr(cr, "factor")
	if(is.null(factor)) factor = rep("foo", length(sample_id))
	gene = genes(txdb)

	if(on == "tss") {
		if(by == "gene") {
			qqcat("extracting gene tss\n")
			tss = promoters(gene, upstream = 1, downstream = 0)
			tss = tss[names(tss) %in% cr$gene_id]
			mapping_column = "gene_id"
		} else {
			qqcat("extracting nearest tx tss\n")
			tx = transcripts(txdb)
			l = !(tx$tx_name %in% gene$gene_id)
			tx = tx[l]
			names(tx) = tx$tx_name
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
		if(by != "gene") stop("'by' should be 'gene' if 'on' is not 'tss'.")
		qqcat("extracting gene body\n")
		gene = gene[names(gene) %in% cr$gene_id]

		target = gene
		target_ratio = 0.6
		mapping_column = "gene_id"
		axis_name = c("-5KB", "TSS", "TES", "5KB")
	}

	target = sort(target)

	mat = normalizeToMatrix(cr, target, mapping_column = mapping_column,
	        extend = 5000, mean_mode = "absolute", w = 50, target_ratio = target_ratio, trim = 0)

	l = rowSums(mat) > 0
	mat = mat[l, ]
	target = target[l]

	mat_cgi = normalizeToMatrix(cgi, target,
	        extend = 5000, mean_mode = "absolute", w = 50, target_ratio = target_ratio)

	meth_mat_list = enrich_with_methylation(target, sample_id, factor, target_ratio = target_ratio)

	if(length(hm_list) > 0) {
		lt = enrich_with_histone_mark(target, hm_list, sample_id, factor, return_arr = TRUE, target_ratio = target_ratio)
		sample = names(hm_list)
		arr = lt[[1]]
		if(length(hm_list) >= 5) {
			# detect regions that histone marks correlate to expression
			expr2 = expr[target$gene_id, intersect(colnames(expr), sample)]
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
			l1 = cor_p_mat < hm_cor_p_cutoff & cor_mat > 0
			cor_mat[l1] = 1
			l2 = cor_p_mat < hm_cor_p_cutoff & cor_mat < 0
			cor_mat[l2] = -1 
			cor_mat[!(l1 | l2)] = 0
			cor_mat = copyAttr(mat, cor_mat)
		}

		hist_mat_list = lt[[2]]
	}

	base_expr = rowMeans(expr[target$gene_id, , drop = FALSE])
	base_expr_label = ifelse(base_expr > mean(base_expr, trim = 0.1), "high", "low")

	expr2 = t(apply(expr[target$gene_id, , drop = FALSE], 1, function(x) {
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
	ht_list = EnrichedHeatmap(mat, col = c("white", cr_col), name = cr_name, cluster_rows = TRUE, show_row_dend = FALSE,
	              top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = cr_col))), 
	              top_annotation_height = unit(2, "cm"), column_title = "meth <> expr", axis_name = axis_name)
	gap = unit(1, "cm")
	if(length(hm_list) >= 5) {
	    ht_list = ht_list + EnrichedHeatmap(cor_mat, col = colorRamp2(c(-1, 0, 1), c("darkgreen", "white", "red")), name = qq("@{hm_name}"), 
	          	column_title = qq("@{hm_name} <> expr\n@{length(hm_list)} samples"))
	    gap = unit.c(gap, unit(3, "mm"))
	}

	gl = width(gene[target$gene_id])
	gl[gl > quantile(gl, 0.99)] = quantile(gl, 0.99)
	ht_list = ht_list +
	          EnrichedHeatmap(mat_cgi, col = c("white", "darkorange"), name = "CGI",
	          	  top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = "darkorange"))), 
	              top_annotation_height = unit(2, "cm"), column_title = "CGI", axis_name = axis_name) +
	          {if(by == "gene") rowAnnotation(gene_len = row_anno_barplot(gl, axis = TRUE), width = unit(1, "cm"))
	           else rowAnnotation(tx_len = row_anno_barplot(width(tx[names(target)]), axis = TRUE), width = unit(1, "cm"))}
	gap = unit.c(gap, unit(1, "cm"), unit(3, "mm"))

	if(show_expr)   {
		ht_list = ht_list + Heatmap(rel_expr, name = "rel_expr", width = unit(20, "mm"), show_row_names = FALSE, cluster_columns = FALSE)
		gap = unit.c(gap, unit(3, "mm"))
	}
	ht_list = ht_list + Heatmap(base_expr, name = "base_e", width = unit(5, "mm"), show_row_names = FALSE)
	gap = unit.c(gap, unit(3, "mm"))

	if(length(hm_list) > 0) {
		ymin = min(sapply(hist_mat_list, function(mat) min(colMeans(mat, na.rm = TRUE))), na.rm = TRUE)
		ymax = max(sapply(hist_mat_list, function(mat) max(colMeans(mat, na.rm = TRUE))), na.rm = TRUE)
		if(ymax > ymin) {
			for(type in names(hist_mat_list)) {
				ht_list = ht_list + EnrichedHeatmap(hist_mat_list[[type]], col = c("white", "purple"), name = qq("hist_@{type}"),
					column_title = qq("@{hm_name}_@{type}"), axis_name = axis_name, show_heatmap_legend = type == names(hist_mat_list)[1],
					heatmap_legend_param = list(title = "hm_density"),
					top_annotation = HeatmapAnnotation(lines1 = anno_enriched(ylim = c(ymin, ymax), gp = gpar(col = "purple"))))
				gap = unit.c(gap, unit(1, "cm"))
			}
		}
	}
	for(type in names(meth_mat_list)) {
		ymin = min(sapply(meth_mat_list, function(mat) min(colMeans(mat, na.rm = TRUE))), na.rm = TRUE)
		ymax = max(sapply(meth_mat_list, function(mat) max(colMeans(mat, na.rm = TRUE))), na.rm = TRUE)
		ht_list = ht_list + EnrichedHeatmap(meth_mat_list[[type]], col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
			name = qq("meth_@{type}"), column_title = qq("meth_@{type}"), axis_name = axis_name, show_heatmap_legend = type == names(meth_mat_list)[1],
				heatmap_legend_param = list(title = "methylation"),
				top_annotation = HeatmapAnnotation(lines1 = anno_enriched(ylim = c(ymin, ymax), gp = gpar(col = "red"))))
		gap = unit.c(gap, unit(1, "cm"))
	}

	draw(ht_list, gap = gap, heatmap_legend_side = "bottom", padding = unit(c(2, 2, 2, 15), "mm"), ...)
	if(by == "gene") {
		decorate_annotation("gene_len", {
			grid.text("gene_len", unit(0.5, "npc"), unit(1, "npc") +unit(2, "mm"), just = "left", rot = 90)
		})
	} else {
		decorate_annotation("tx_len", {
			grid.text("tx_len", unit(0.5, "npc"), unit(1, "npc") +unit(2, "mm"), just = "left", rot = 90)
		})
	}
}


# == title
# Enriched heatmap on TSS-CGI
#
# == param
# -cr cr
# -cgi CpG Islands
# -txdb txdb
# -expr expr
# -hm_list hm_list
# -hm_name hm_name
# -by by
# -... pass to
#
enriched_heatmap_list_on_tss_cgi = function(cr, cgi, txdb, expr, hm_list = NULL, hm_name = NULL, by = "gene", ...) {

	sample_id = attr(cr, "sample_id")
	factor = attr(cr, "factor")
	if(is.null(factor)) factor = rep("foo", length(sample_id))
	gene = genes(txdb)

	if(by == "gene") {
		qqcat("extracting gene tss\n")
		tss = promoters(gene, upstream = 1, downstream = 0)
		tss = tss[names(tss) %in% cr$gene_id]
		expr = expr[names(tss), , drop = FALSE]
	} else {
		qqcat("extracting nearest tx tss\n")
		tx_list = transcriptsBy(txdb, "gene")
		gene_name = names(tx_list)
		tx = unlist(tx_list)
		tx$gi = rep(gene_name, times = sapply(tx_list, length))
		tx = tx[tx$tx_name != tx$gi]
		tx = tx[tx$gi %in% cr$gene_id]
		tss = promoters(tx, upstream = 1, downstream = 0)
		expr = expr[tx$gi, , drop = FALSE]
	}

	strand(cr) = strand(gene[cr$gene_id])

	dist = distanceToNearest(cgi, tss)
	l = dist@elementMetadata@listData$distance < 5000
	mtch = as.matrix(dist)[l, ]

	cgi = cgi[mtch[,1]]
	cgi$nearest_tss = names(tss[mtch[,2]])

	cgi_1kb = cgi; start(cgi_1kb) = start(cgi) - 1000; end(cgi_1kb) = end(cgi) + 1000
	mtch = as.matrix(findOverlaps(cgi_1kb, tss))
	high_expressed_genes = tapply(mtch[, 2], mtch[, 1], function(ind) {
			if(length(ind) == 1) {
				return(ind)
			} else {
				which.max(rowMeans(expr[ind, , drop = FALSE]))[1]
			}
		})
	cgi$nearest_expressed_tss = cgi$nearest_tss
	cgi$nearest_expressed_tss[as.numeric(names(high_expressed_genes))] = names(tss)[high_expressed_genes]

	cgi_extend = cgi; start(cgi_extend) = start(cgi) - 5000; end(cgi_extend) = end(cgi) + 5000
	mtch = as.matrix(findOverlaps(cr, cgi_extend))
	# cgi gene should be the same as cr gene
	l = cr[mtch[, 1]]$gene_id == cgi_extend[mtch[,2]]$nearest_expressed_tss
	mtch = mtch[l, , drop = FALSE]
	cgi2 = cgi[unique(mtch[, 2])]
	cgi_extend = cgi_extend[unique(mtch[, 2])]

	qqcat("There are @{length(cgi2)} CGIs that exists within TSS +- 5kb\n")

	mat3 = normalizeToMatrix(cr[unique(mtch[, 1])], cgi2, extend = 5000, w = 50, trim = 0)
	l = rowSums(mat3) > 0
	mat3 = mat3[l, ]
	cgi2 = cgi2[l]

	n_tss = countOverlaps(cgi_extend, tss)

	dist = distanceToNearest(cgi2, tss)

	strd = as.vector(strand(gene[cgi2$nearest_expressed_tss]))
	strd = factor(strd, levels = c("-", "+"))

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

	cgi_width = width(cgi2[l])
	cgi_width[cgi_width > quantile(cgi_width, 0.99)] = quantile(cgi_width, 0.99)

	ht_list = EnrichedHeatmap(mat3, col = c("white", cr_col), cluster_rows = TRUE, show_row_dend = FALSE, split = strd,
	                top_annotation = HeatmapAnnotation(lines1 = anno_enriched(gp = gpar(col = c("darkgreen", "red")))), 
	                top_annotation_height = unit(2, "cm"), column_title = qq("@{cr_name} ~ cgi"), axis_name = c("-5kb", "start", "end", "5kb"))+
	            rowAnnotation(cgi_width = row_anno_barplot(cgi_width, axis = TRUE),
	                width = unit(1, "cm")) +
	            Heatmap(strd, name = "strand", col = c("+" = "red", "-" = "darkgreen"),
	                width = unit(5, "mm")) +
	            rowAnnotation(nearby_n_tss = row_anno_barplot(n_tss, axis = TRUE), width = unit(1, "cm"))


	if(length(hm_list) > 0) {
		hist_mat_list = enrich_with_histone_mark(cgi2, hm_list, sample_id, factor)
		ymin = min(sapply(hist_mat_list, function(mat) min(colMeans(mat, na.rm = TRUE))))
		ymax = max(sapply(hist_mat_list, function(mat) max(colMeans(mat, na.rm = TRUE))))
		if(ymax > ymin) {
			for(type in names(hist_mat_list)) {
				ht_list = ht_list + EnrichedHeatmap(hist_mat_list[[type]], col = c("white", "purple"), name = qq("hist_@{type}"),
					column_title = type, axis_name = c("-5kb", "start", "end", "5kb"), show_heatmap_legend = type == names(hist_mat_list)[1],
					heatmap_legend_param = list(title = "hm_density"),
					top_annotation = HeatmapAnnotation(lines1 = anno_enriched(ylim = c(ymin, ymax), gp = gpar(col = "purple"))))
			}
		}
	}

	meth_mat_list = enrich_with_methylation(cgi2, sample_id, factor)
	for(type in names(meth_mat_list)) {
		ymin = min(sapply(meth_mat_list, function(mat) min(colMeans(mat, na.rm = TRUE))))
		ymax = max(sapply(meth_mat_list, function(mat) max(colMeans(mat, na.rm = TRUE))))
		ht_list = ht_list + EnrichedHeatmap(meth_mat_list[[type]], col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
			name = qq("meth_@{type}"), column_title = type, axis_name = c("-5kb", "start", "end", "5kb"), show_heatmap_legend = type == names(meth_mat_list)[1],
				heatmap_legend_param = list(title = "methylation"),
				top_annotation = HeatmapAnnotation(lines1 = anno_enriched(ylim = c(ymin, ymax), gp = gpar(col = "red"))))
	}
	draw(ht_list, heatmap_legend_side = "bottom", ...)
	decorate_annotation("cgi_width", {
		grid.text("cgi_width", unit(0.5, "npc"), unit(1, "npc") + unit(2, "mm"), just = "left", rot = 90)
	})
	decorate_annotation("nearby_n_tss", {
		grid.text("nearby_tss", unit(0.5, "npc"), unit(1, "npc") + unit(2, "mm"), just = "left", rot = 90)
	})
}

# == title
# Enriched heatmap on a certain gf
#
# == param
# -cr cr
# -gf genomic features
# -hm_list hm_list
# -hm_name hm_name
# -... pass to
#
enriched_heatmap_list_on_genomic_features = function(cr, gf, hm_list = NULL, hm_name = NULL, ...) {

	sample_id = attr(cr, "sample_id")
	factor = attr(cr, "factor")
	if(is.null(factor)) factor = rep("foo", length(sample_id))

	gf_extend = gf; start(gf_extend) = start(gf) - 5000; end(gf_extend) = end(gf) + 5000
	mtch = as.matrix(findOverlaps(cr, gf_extend))

	gf2 = gf[unique(mtch[, 2])]
	
	mat3 = normalizeToMatrix(cr[unique(mtch[, 1])], gf2, extend = 5000, w = 50, trim = 0)
	l = rowSums(mat3) > 0
	mat3 = mat3[l, ]
	gf2 = gf2[l]

	qqcat("There are @{length(gf2)} gfs that exists\n")

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

	gf_width = width(gf2[l])
	gf_width[gf_width > quantile(gf_width, 0.99)] = quantile(gf_width, 0.99)

	ht_list = EnrichedHeatmap(mat3, col = c("white", cr_col), cluster_rows = TRUE, show_row_dend = FALSE,
	                top_annotation = HeatmapAnnotation(lines1 = anno_enriched()), 
	                top_annotation_height = unit(2, "cm"), column_title = qq("@{cr_name} ~ gf"), axis_name = c("-5kb", "start", "end", "5kb"))


	if(length(hm_list) > 0) {
		hist_mat_list = enrich_with_histone_mark(gf2, hm_list, sample_id, factor)
		ymin = min(sapply(hist_mat_list, function(mat) min(colMeans(mat, na.rm = TRUE))))
		ymax = max(sapply(hist_mat_list, function(mat) max(colMeans(mat, na.rm = TRUE))))
		if(ymax > ylim) {
			for(type in names(hist_mat_list)) {
				ht_list = ht_list + EnrichedHeatmap(hist_mat_list[[type]], col = c("white", "purple"), name = qq("hist_@{type}"),
					column_title = type, axis_name = c("-5kb", "start", "end", "5kb"), show_heatmap_legend = type == names(hist_mat_list)[1],
					heatmap_legend_param = list(title = "hm_density"),
					top_annotation = HeatmapAnnotation(lines1 = anno_enriched(ylim = c(ymin, ymax), gp = gpar(col = "purple"))))
			}
		}
	}

	meth_mat_list = enrich_with_methylation(gf2, sample_id, factor)
	for(type in names(meth_mat_list)) {
		ymin = min(sapply(meth_mat_list, function(mat) min(colMeans(mat, na.rm = TRUE))))
		ymax = max(sapply(meth_mat_list, function(mat) max(colMeans(mat, na.rm = TRUE))))
		ht_list = ht_list + EnrichedHeatmap(meth_mat_list[[type]], col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), 
			name = qq("meth_@{type}"), column_title = type, axis_name = c("-5kb", "start", "end", "5kb"), show_heatmap_legend = type == names(meth_mat_list)[1],
				heatmap_legend_param = list(title = "methylation"),
				top_annotation = HeatmapAnnotation(lines1 = anno_enriched(ylim = c(ymin, ymax), gp = gpar(col = "red"))))
	}
	draw(ht_list, heatmap_legend_side = "bottom", ...)
}


enrich_with_histone_mark = function(target, hm_list, sample_id, factor, target_ratio = 0.1, return_arr = FALSE, 
	value_column = "density", extend = 5000, mean_mode = "w0", w = 50, ...) {

	target_name = deparse(substitute(target))
	# 200 is number of windows (5000 + 5000)/50

	if(is.null(factor)) factor = rep("foo", length(sample_id))

	sample = names(hm_list)
	flag = 0
	for(i in seq_along(sample)) {
		qqcat("@{sample[i]}: normalize histone modifications to @{target_name}.\n")
	    tm = normalizeToMatrix(hm_list[[i]], target, target_ratio = target_ratio, 
	    	value_column = value_column, extend = extend, mean_mode = mean_mode, w = w, trim = c(0, 0.01), ...)
	    if(!flag) {
	    	arr = array(dim = c(length(target), dim(tm)[2], length(hm_list)))
	    	flag = 1
	    }
	    arr[, , i] = tm
	}

	########### matrix for each subgroup
	hist_mat_list = list()
	for(type in unique(factor)) {
	    l = sample %in% sample_id[factor == type]
	    if(sum(l) == 0) {
	    	hist_mat_list[[type]] = NULL
	    } else {
	    	hist_mat_list[[type]] = apply(arr[, , l,drop = FALSE], c(1, 2), mean, na.rm = TRUE)
	    	hist_mat_list[[type]] = copyAttr(tm, hist_mat_list[[type]])
	    }
	}

	if(return_arr) {
		return(list(arr = arr, list = hist_mat_list))
	} else {
		return(hist_mat_list)
	}
}

enrich_with_methylation = function(target, sample_id, factor, target_ratio = 0.1,
	extend = 5000, w = 50, mean_mode = "absolute", empty_value = NA, smooth = TRUE, ...) {

	# extrace methylation value which in [-5k, 5k] from TSS, and calculate mean methylation in each subgroup
	target_extend = GRanges(seqnames = seqnames(target), ranges = IRanges(start(target)-extend, end(target)+extend))
	# process raw methylation data
	meth_gr = GRanges()
	for(chr in sort(unique(as.vector(seqnames(target))))) {
	    cat(chr, "\n")
	    methylation_hooks$set(chr)
	    gr = methylation_hooks$GRanges()
	    mtch = as.matrix(findOverlaps(gr, target_extend))
	    ind = unique(mtch[, 1])
	    mm = do.call("cbind", lapply(unique(factor), function(s) rowMeans(methylation_hooks$meth(row_index = ind, col_index = sample_id[factor == s]), na.rm = TRUE)))
	    gr = gr[ind]
	    mcols(gr) = mm
	    meth_gr = c(meth_gr, gr)
	}
	rm(gr)
	gc(verbose = FALSE)
	mm = mcols(meth_gr)
	colnames(mm) = unique(factor)
	mcols(meth_gr) = mm
	meth_mat_list = list()
	for(type in unique(factor)) {
		qqcat("@{type}: normalize methylation to target.\n")
	    mat = normalizeToMatrix(meth_gr, target, value_column = type, target_ratio = target_ratio, extend = extend, w = w, mean_mode = mean_mode, empty_value = empty_value, 
	    	smooth = smooth, ...)
	    mat[mat > 1] = 1
	    mat[mat < 0] = 0
	    meth_mat_list[[type]] = mat
	}

	return(meth_mat_list)
}

