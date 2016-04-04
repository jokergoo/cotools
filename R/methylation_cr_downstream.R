
# == title
# QC on crs
#
# == param
# -chromosome chromosome
# -template template
# 
cr_qc = function(chromosome = paste0("chr", 1:22), template) {
	
	cutoff = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
	diameter = c(0, 0.1, 0.2, 0.3, 0.4)
	
	n = matrix(0, nr = length(cutoff), nc = 2)
	n = rep(list(n), length(diameter))
	w = n

	for(k in seq_along(chromosome)) {
	    chr = chromosome[k]
	    qqcat("reading @{chr}...\n")
	    cr = readRDS(qq(template))
	    
	    l = cr$corr > 0; l[is.na(l)] = FALSE
	    pos_cr = cr[l]
	    l = cr$corr < 0; l[is.na(l)] = FALSE
	    neg_cr = cr[l]
	    
	    for(i in seq_along(cutoff)) {
	    	for(j in seq_along(diameter)) {
		        n[[j]][i, 1] = n[[j]][i, 1] + sum(neg_cr$meth_diameter >= diameter[j] & neg_cr$corr <= -cutoff[i])
		        n[[j]][i, 2] = n[[j]][i, 2] + sum(pos_cr$meth_diameter >= diameter[j] & pos_cr$corr >= cutoff[i])
		        w[[j]][i, 1] = w[[j]][i, 1] + sum(width(neg_cr[neg_cr$meth_diameter >= diameter[j] & neg_cr$corr <= -cutoff[i]]))
		        w[[j]][i, 2] = w[[j]][i, 2] + sum(width(pos_cr[pos_cr$meth_diameter >= diameter[j] & pos_cr$corr >= cutoff[i]]))
		    }
	    }
	}

	r1 = max(sapply(n, function(x) x[,1]/x[,2]))
	r2 = max(sapply(w, function(x) x[,1]/x[,2]))

	par(mfrow = c(length(diameter), 2), mar = c(4, 4, 2, 1), xpd = NA)
	for(i in seq_along(n)) {
		plot(n[[i]][, 1]/n[[i]][, 2], ylim = c(0, r1), type = "b", axes = FALSE, ylab = "#neg/#pos", xlab = "abs(correlation)")
		axis(side = 1, at = seq_along(cutoff), labels = paste0(">",cutoff))
		axis(side = 2); box()
		text(length(cutoff)/2, r1, paste0("diameter > ", diameter[i]), adj = c(0.5, 1))
		text(seq_along(n[[i]][, 1]), n[[i]][, 1]/n[[i]][, 2], qq("@{n[[i]][, 1]}\n@{n[[i]][, 2]}", collapse=FALSE), srt = 45, col = "red")

		plot(w[[i]][, 1]/w[[i]][, 2], ylim = c(0, r2), type = "b", axes = FALSE, ylab = "width(neg)/width(pos)", xlab = "abs(correlation)")
		axis(side = 1, at = seq_along(cutoff), labels = paste0(">",cutoff))
		axis(side = 2); box()
		text(length(cutoff)/2, r2, paste0("diameter > ", diameter[i]), adj = c(0.5, 1))
		text(seq_along(w[[i]][, 1]), w[[i]][, 1]/w[[i]][, 2], qq("@{w[[i]][, 1]}\n@{w[[i]][, 2]}", collapse=FALSE), srt = 45, col = "red")
	}

	mat_n = matrix(nrow = length(diameter), ncol = length(cutoff))
	rownames(mat_n) = diameter
	colnames(mat_n) = cutoff
	for(i in seq_along(n)) {
		mat_n[i, ] = n[[i]][, 1]/n[[i]][,2 ]
	}

	mat_w = matrix(nrow = length(diameter), ncol = length(cutoff))
	rownames(mat_w) = diameter
	colnames(mat_w) = cutoff
	for(i in seq_along(n)) {
		mat_w[i, ] = w[[i]][, 1]/w[[i]][,2 ]
	}

	# Heatmap(mat_n, cluster_rows = FALSE, cluster_columns = FALSE, name = "#neg/#pos")
	# Heatmap(mat_w, cluster_rows = FALSE, cluster_columns = FALSE, name = "w(neg)/w(pos)")

	return(invisible(list(n = n, w = w)))
}

# == title
# enrichment between cr and genomic feature
#
# == param
# -cr cr
# -gf_list gf_list
# -species species
# -chromosome chromosomes
# 
cr_overlap_to_genomic_features = function(cr, gf_list, species = NULL, chromosome = paste0("chr", 1:22)) {

	neg_cr = cr[cr$corr < 0]
	pos_cr = cr[cr$corr > 0]

	neg_cr_reduced = reduce(neg_cr)
	pos_cr_reduced = reduce(pos_cr)
	pct_mat = matrix(0, nrow = length(gf_list), ncol = 3)
	rownames(pct_mat) = names(gf_list)
	colnames(pct_mat) = c("genome", "neg_cr", "pos_cr")
	overlap_mat = pct_mat
	chr_sum_len = sum(read.chromInfo(species = species)$chr.len[chromosome])
	for(i in seq_along(gf_list)) {
	    gr = gf_list[[i]]
	    gr = reduce(gr)
	    gr = gr[seqnames(gr) %in% chromosome]
	    pct_mat[i, 1] = sum(width(gr))/chr_sum_len
	    overlap_mat[i, 1] = sum(width(gr))
	    mtch = as.matrix(findOverlaps(neg_cr_reduced, gr))
	    pct_mat[i, 2] = sum(width(pintersect(neg_cr_reduced[mtch[,1]], gr[mtch[,2]])))/sum(width(neg_cr_reduced))
	    overlap_mat[i, 2] = sum(width(pintersect(neg_cr_reduced[mtch[,1]], gr[mtch[,2]])))
	    mtch = as.matrix(findOverlaps(pos_cr_reduced, gr))
	    pct_mat[i, 3] = sum(width(pintersect(pos_cr_reduced[mtch[,1]], gr[mtch[,2]])))/sum(width(pos_cr_reduced))
	    overlap_mat[i, 3] = sum(width(pintersect(pos_cr_reduced[mtch[,1]], gr[mtch[,2]])))
	}

	par(mar = c(8, 4, 4, 1))
	mat = log2(t(pct_mat[, 2:3]/pct_mat[,1]))
	colnames(mat) = NULL
	pos = barplot(mat, beside = TRUE, ylim = range(mat)*1.1, axes = FALSE, ann = FALSE, col = c("green", "red"), 
	    ylab = "fold change", main = expression(frac("width(intersect(gf, cr))", "width(cr)")/frac("width(intersect(gf, genome))", "width(genome)")))
	axis(side = 2, at = seq(-3, 3), labels = c(0.125, 0.25, 0.5, 1, 2, 4, 8))
	par(xpd = NA)
	text(colMeans(pos), min(mat)*1.2, names(gf_list), srt = 90, adj = c(1, 0.5))
	box()
	legend("bottomleft", pch = 15, col = c("green", "red"), legend = c("neg_cr", "pos_cr"))

	par(mar = c(8, 4, 4, 1))
	mat = t(overlap_mat[, 2:3])
	colnames(mat) = NULL
	pos = barplot(mat, beside = TRUE, ylim = range(mat)*1.1, ann = FALSE, col = c("green", "red"), 
	    ylab = "bp", main = "Intersection between gf and cr")
	par(xpd = NA)
	text(colMeans(pos), 0, names(gf_list), srt = 90, adj = c(1.2, 0.5))
	box()
	legend("topright", pch = 15, col = c("green", "red"), legend = c("neg_cr", "pos_cr"))

	return(invisible(list(pct = pct_mat, overlap = overlap_mat)))
}

# == title
# Hilbert curve visualization of CRs
#
# == param
# -cr cr
# -template template
# -txdb txdb
# -chromosome chromosome
# -merge_chr whether merge chromsomes into one plot
#
cr_hilbert = function(cr, template, txdb, chromosome = paste0("chr", 1:22), merge_chr = TRUE) {

	gene = genes(txdb)
	
	chr_len = read.chromInfo()$chr.len

	if(!missing(cr)) {
		cm = ColorMapping(levels = c("neg", "pos"), colors = c("darkgreen", "red"))
		lgd = color_mapping_legend(cm, title = "type", plot = FALSE)
		cr = cr[!is.na(cr$corr)]
		if(merge_chr) {
			hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = "CR for all chromosomes", legend = lgd)
		    hc_layer(hc, gene, col = "#F0F0FF")
		    hc_layer(hc, cr, col = ifelse(cr$corr > 0, "red", "darkgreen"), mean_mode = "absolute")
		    hc_map(hc, add = TRUE, fill = NA, border = "grey")
		    hc_map(hc, title = "map for all chromosomes")
		} else {
			# grid.newpage()
			# pushViewport(viewport(layout = grid.layout(nr = 4, nc = 6)))
			for(i in seq_along(chromosome)) {
			    chr = chromosome[i]
			    cat(chr, "\n")
			    cr_subset = cr[seqnames(cr) == chr]
			    gene_subset = gene[seqnames(gene) == chr]
			    # pushViewport(viewport(layout.pos.row = ceiling(i/6), layout.pos.col = i - (ceiling(i/6)-1)*6))
			    hc = HilbertCurve(s = 1, e = max(chr_len), mode = "pixel", level = 10, title = chr, legend = lgd)
			    hc_layer(hc, ranges(reduce(gene_subset)), col = "#F0F0FF")
			    hc_layer(hc, ranges(cr_subset), col = ifelse(cr_subset$corr > 0, "red", "darkgreen"), mean_mode = "absolute")
			    # upViewport()
			}
			# upViewport()
		}
	}

	if(!missing(template)) {
		## all cr windows
		col_fun = colorRamp2(c(-1, 0, 1), c("green", "white", "red"))
		cm = ColorMapping(col_fun = col_fun)
		lgd = color_mapping_legend(cm, title = "type", plot = FALSE)
		if(merge_chr) {
			cr = GRanges()
			for(i in seq_along(chromosome)) {
				chr = chromosome[i]
				cr = c(cr, readRDS(qq(template)))
			}
			cr = cr[!is.na(cr$corr)]
			hc = GenomicHilbertCurve(chr = chromosome, mode = "pixel", level = 10, title = "cr for all chromosomes", legend = lgd)
		    hc_layer(hc, cr, col = col_fun(cr$corr), mean_mode = "absolute")
		    hc_map(hc, add = TRUE, fill = NA, border = "grey")
		   	hc_map(hc, title = "map for all chromosomes")
		} else {
			# grid.newpage()
			# pushViewport(viewport(layout = grid.layout(nr = 4, nc = 6)))
			for(i in seq_along(chromosome)) {
			    chr = chromosome[i]
			    cat(chr, "\n")
			    cr = readRDS(qq(template))
			    cr = cr[!is.na(cr$corr)]
			    # pushViewport(viewport(layout.pos.row = ceiling(i/6), layout.pos.col = i - (ceiling(i/6)-1)*6))
			    hc = HilbertCurve(s = 1, e = max(chr_len), mode = "pixel", level = 10, title = chr, legend = lgd)
			    hc_layer(hc, ranges(cr), col = col_fun(cr$corr), mean_mode = "absolute")
			    # upViewport()
			}
			# upViewport()
		}
	}
}

# == title
# compare methylation between smoothed and raw methylation data
#
# == param
# -cr cr
# -chr chromosome
# -start start position
# -end end position
# -x cr with smoothed methylation
# -x2 cr with raw methylation
#
compare_meth = function(cr, chr, start, end, x = NULL, x2 = NULL) {
	
	sample_id = attr(cr, "sample_id")
	factor = attr(cr, "factor")
	col = attr(cr, "col")
	n_sample = length(sample_id)
	if(is.null(col)) col = sample(n_sample, n_sample)
	
	methylation_hooks$set(chr)

	ind = extract_sites(start, end, methylation_hooks$site())
	meth = methylation_hooks$meth(row_index = ind, col_index = sample_id)
	raw = methylation_hooks$raw(row_index = ind, col_index = sample_id)
	cov = methylation_hooks$coverage(row_index = ind, col_index = sample_id)
	site = extract_sites(start, end, methylation_hooks$site(), index = FALSE)

	par(mfrow = c(5 + (!is.null(x)) + (!is.null(x2)), 1), mar = c(1, 4, 1, 1))

	matplot(site, meth, type = "l", col = col[factor], lty = 1, 
	    ylab = "smoothed meth", xlab = NULL)
	legend("bottomleft", lty = 1, col = col, legend = names(col))

	if(!is.null(x)) plot((start(x)+end(x))/2, x$corr, xlim = range(site), ylab = "corr\n(from smoothed)", xlab = NULL, type = "l")

	matplot(site, raw, type = "l", col = col[factor], lty = 1, 
	    ylab = "raw meth", xlab = NULL)

	if(!is.null(x2)) plot((start(x2)+end(x2))/2, x2$corr, xlim = range(site), ylab = "corr\n(from raw)", xlab = NULL, type = "l")

	plot(site, xlim = range(site), ylim = c(0, 1), type = "n", 
	    ylab = "raw meth (cov >= 5)", xlab = NULL)
	for(i in seq_len(ncol(raw))) {
	    xx = raw[, i]
	    yy = cov[, i]
	    lines(site[yy >= 5], xx[yy >= 5], col = col[factor[i]])
	}
	plot(site, xlim = range(site), ylim = c(0, 1), type = "n", 
	    ylab = "raw meth (cov >= 10)", xlab = NULL)
	for(i in seq_len(ncol(raw))) {
	    xx = raw[, i]
	    yy = cov[, i]
	    lines(site[yy >= 10], xx[yy >= 10], col = col[factor[i]])
	}
	par(mar = c(4, 4, 1, 1))
	plot(site, rowMeans(cov), type = 'h', ylab = "CpG coverage", xlab = "CpG sites")

}

