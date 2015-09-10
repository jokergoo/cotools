

# == title
# correlated regions in extended gene model
#
# == param
# -site CpG sites
# -meth methylation matrix corresponding to ``site``
# -expr expression for current gene
# -chr chromosome
# -type negative or positive
# -p_cutoff cutoff for p-value for correlation test
# -other_filter other filter on methylation and expression
#
# == detail
# wrapper on `GenomicRegions::genomic_regions_finder`
correlated_regions_per_gene = function(site, meth, cov, expr, chr, cov_cutoff = 3, min_dp = 4,
	cor_method = "spearman", window_size = 5, factor = NULL, max_width = 10000) {

	if(ncol(meth) != length(expr)) {
		stop("number of columsn of `meth` should be same as length of `expr`.\n")
	}

	index = seq(1, length(site), by = window_size)
	
	i = seq_len(length(index) - 1)
	ir = IRanges(site[index[i]], site[index[i+1]-1])

	m = lapply(index[i], function(x) {
		ind = x+0:(window_size-1)
		meth_m = meth[ind, , drop = FALSE]
		cov_m = cov[ind, , drop = FALSE]
		sapply(seq_len(ncol(meth_m)), function(i) {
			xm = meth_m[, i]
			ym = cov_m[, i]
			mean(xm[ym >= cov_cutoff])
		})
	})
	m = do.call('rbind', m)
	colnames(m) = paste0("mean_meth_", colnames(meth))
	corr = apply(m, 1, function(x) {
		l = !is.na(x)
		if(sum(l) < min_dp) return(NA)
		cor(x[l], expr[l], method = cor_method)
	})
	corr_p = suppressWarnings(apply(m, 1, function(x) {
		l = !is.na(x)
		if(sum(l) < min_dp) return(NA)
		cor.test(x[l], expr[l], method = cor_method)$p.value
	}))

	if(!is.null(factor)) {
		factor = as.vector(factor)
		meth_anova = apply(m, 1, function(x) {
			l = !is.na(x)
			data = data.frame(value = x[l], class = factor[l], stringsAsFactors = FALSE)
			if(length(unique(data$class)) < 2) return(NA)
			if(any(table(data$class) < 2)) return(NA)
			oneway.test(value ~ class, data = data)$p.value
		})
		meth_diameter = apply(m, 1, function(x) {
			l = !is.na(x)
			if(any(table(factor[l]) < 2)) return(NA)
			diameter(as.vector(tapply(x[l], factor[l], mean)))
		})
	}
	gr = GRanges(seqnames = rep(chr, length(ir)),
		    ranges = ir)
	if(is.null(factor)) {
		df = DataFrame(n = window_size,
		    m, # mean methylation
		    corr = corr,
		    corr_p = corr_p)	
	} else {
		df = DataFrame(n = window_size,
		    m, # mean methylation
		    corr = corr,
		    corr_p = corr_p,
		    meth_anova = meth_anova,
		    meth_diameter = meth_diameter)
	}
	mcols(gr) = df

	return(gr)
}

# == title
# correlated regions
#
# == param
# -sample_id sample id
# -expr expression matrix
# -txdb ``transcritpDb`` object
# -chr chromosome
# -extend extension of gene model, both upstream and downstream
# -cov_filter function to filter on coverage
# -type type of correlation
# -p_cutoff cutoff for correlation test
# -factor if specified, only regions show differential methylation are keeped
# -mean_diff mean different between groups
#
# == detail
# based on `correlated_regions_per_gene`
correlated_regions = function(sample_id, expr, txdb, chr, extend = 50000,
	cov_filter = function(x) sum(x > 0) > length(x)/2,
	cor_method = "spearman", factor = NULL, window_size = 5, max_width = 10000,
	raw_meth = FALSE, cov_cutoff = 3, min_dp = 4) {

	qqcat("extracting gene model (extend = @{extend}, chr = @{chr})...\n")
	gene = genes(txdb)
	tx_list = transcriptsBy(txdb, by = "gene")

	gene = gene[seqnames(gene) == chr]

	g = intersect(rownames(expr), names(gene))
	expr = expr[g, , drop = FALSE]
	gene = gene[g]
	tx_list = tx_list[g]

	genemodel = gene
	start(genemodel) = start(genemodel) - extend
	end(genemodel) = end(genemodel) + extend
	s = start(genemodel)
	start(genemodel) = ifelse(s > 0, s, 1)

	expr = expr[, sample_id, drop = FALSE]

	all_gi = rownames(expr)
	n_gene = length(all_gi)

	methylation_hooks$set(chr)
	site = methylation_hooks$site()
	if(raw_meth) {
		meth = methylation_hooks$raw(col_index = sample_id)
	} else {
		meth = methylation_hooks$meth(col_index = sample_id)
	}
	cov = methylation_hooks$coverage(col_index = sample_id)

	if(!is.null(cov_filter)) {
		
		l = apply(cov, 1, cov_filter)
		site = site[l]
		meth = meth[l, , drop = FALSE]
		cov = cov[l, , drop = FALSE]
	}

	op = qq.options("cat_prefix")

	if(!raw_meth) cov_cutoff = 0
	
	res = GRanges()
	for(i in seq_len(n_gene)) {
			
		# current gene name
		gi = all_gi[i]

		# expression of current gene
		e = expr[gi, ]

		if(is.function(op)) {
			qq.options(cat_prefix = function(x) {
				qq("@{op()} [@{chr}:@{gi}, @{i}/@{n_gene}]")
			})
		} else {
			qq.options(cat_prefix = qq("@{op} [@{chr}:@{gi}, @{i}/@{n_gene}]"))
		}

		# if gene has low expression in many samples
		if(all(e == 0) || sd(e) == 0)  {
			qqcat("@{gi} has zero expression in all samples, skip\n")
			next
		}

		start = start(genemodel[gi])
		end = end(genemodel[gi])
		gm_site_index = extract_sites(start, end, site, index = TRUE)
		gm_site = site[gm_site_index]
		gm_meth = meth[gm_site_index, sample_id, drop = FALSE]
		gm_cov = cov[gm_site_index, sample_id, drop = FALSE]

		if(length(gm_site) < 10) {
			qqcat("@{gi} has too few cpg sites, skip\n")
			next
		}

		qqcat("...\n")
		gr = correlated_regions_per_gene(gm_site, gm_meth, gm_cov, e, cov_cutoff = cov_cutoff, chr = chr,
			factor = factor, cor_method = cor_method, window_size = window_size, min_dp = min_dp,
			max_width = max_width)
		gr$gene_id = rep(gi, length(gr))

		## distance to gene tss
		tss = promoters(gene[gi], upstream = 1, downstream = 0)
		gene_tss_dist = as.data.frame(distanceToNearest(gr, tss))[, 3]
		if(as.vector(strand(gene[gi])) == "+") {
			gene_tss_dist = ifelse(end(gr) < start(tss), -gene_tss_dist, gene_tss_dist)
		} else if(as.vector(strand(gene[gi])) == "-") {
			gene_tss_dist = ifelse(start(gr) > end(tss), -gene_tss_dist, gene_tss_dist)
		}
		gr$gene_tss_dist = gene_tss_dist

		## distance to tx tss
		tx = tx_list[[gi]]
		tx = tx[tx$tx_name != gi]
		tss = promoters(tx, upstream = 1, downstream = 0)
		dist = as.data.frame(distanceToNearest(gr, tss))
		tx_tss_dist = dist[, 3]
		if(as.vector(strand(gene[gi])) == "+") {
			tx_tss_dist = ifelse(end(gr) < start(tss[dist[,2]]), -tx_tss_dist, tx_tss_dist)
		} else if(as.vector(strand(gene[gi])) == "-") {
			tx_tss_dist = ifelse(start(gr) > end(tss[dist[,2]]), -tx_tss_dist, tx_tss_dist)
		}
		gr$tx_tss_dist = tx_tss_dist
		gr$nearest_tx_tss = tss[dist[,2]]$tx_name

		res = c(res, gr)

	}
	qq.options(cat_prefix = op)
	
	attr(res, "factor") = factor
	attr(res, "col") = col
	attr(res, "cor_method") = cor_method
	attr(res, "extend") = extend
	attr(res, "window_size") = window_size
	attr(res, "sample_id") = sample_id
	attr(res, "cor_method") = cor_method
	attr(res, "cov_filter") = cov_filter
	attr(res, "raw_meth") = raw_meth
	attr(res, "cov_filter") = cov_filter
	attr(res, "min_dp") = min_dp

	return(res)
}

filter_correlated_regions = function(chromosome = paste0("chr", c(1:22, "X")), template, 
	cutoff = 0.05, adj_method = "BH", meth_diameter_cutoff = 0.25,
	anova_cutoff = 0.05) {

	if(length(cutoff) == 1) cutoff = rep(cutoff, 2)

	cat("calculate fdr...\n")
	corr_p = NULL
	meth_anova = NULL
	meth_diameter = NULL
	chr_name = NULL
	corr = NULL
	for(chr in chromosome) {
		qqcat("reading cr for @{chr}\n")
		cr = readRDS(qq(template))
		has_anova = FALSE
		if("meth_anova" %in% colnames(mcols(cr))) {
			has_anova = TRUE
		}

		if(has_anova) {
			cr = cr[(!is.na(cr$corr)) & (!is.na(cr$corr_p)) & (!is.na(cr$meth_anova)) & (!is.na(cr$meth_diameter))]
		} else {
			cr = cr[(!is.na(cr$corr)) & (!is.na(cr$corr_p))]
		}
	    
	    corr_p = c(corr_p, cr$corr_p)
	    corr = c(corr, cr$corr)
	   	if(has_anova) {
	   		meth_anova = c(meth_anova, cr$meth_anova)
	   		meth_diameter = c(meth_diameter, cr$meth_diameter)
	   	}
	    chr_name = c(chr_name, rep(chr, length(cr)))
	}

	if(has_anova) {
		anova_fdr = p.adjust(meth_anova, method = adj_method)
		l = anova_fdr <= cutoff[1] & meth_diameter >= meth_diameter_cutoff
		qqcat("filter out @{sum(!l)}/@{length(l)} by differential methylation.\n")
		corr_fdr = rep(Inf, length(corr_p))
		corr_fdr[l] = p.adjust(corr_p[l], method = adj_method)
		l = l & ifelse(corr > 0, corr_fdr <= cutoff[1], corr_fdr <= cutoff[2])
	} else {
		corr_fdr = p.adjust(corr_p, method = adj_method)
		l = ifelse(corr > 0, corr_fdr <= cutoff[1], corr_fdr <= cutoff[2])
	}

	cat("filter by fdr...\n")
	cr2 = GRanges()
	for(chr in chromosome) {
		qqcat("reading cr for @{chr}\n")
		cr = readRDS(qq(template))

		if(has_anova) {
			cr = cr[(!is.na(cr$corr)) & (!is.na(cr$corr_p)) & (!is.na(cr$meth_anova)) & (!is.na(cr$meth_diameter))]
		} else {
			cr = cr[(!is.na(cr$corr)) & (!is.na(cr$corr_p))]
		}

		lo = chr_name == chr
		cr$corr_fdr = corr_fdr[lo]
		if(has_anova) {
			cr$meth_anova_fdr = anova_fdr[lo]
		}
		cr2 = suppressWarnings(c(cr2, cr[l[lo]]))
	}

	attr(cr2, "factor") = attr(cr, "factor")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "extend") = attr(cr, "extend")
	attr(cr2, "window_size") = attr(cr, "window_size")
	attr(cr2, "sample_id") = attr(cr, "sample_id")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "cov_filter") = attr(cr, "cov_filter")
	attr(cr2, "cov_cutoff") = attr(cr, "cov_cutoff")
	attr(cr2, "raw_meth") = attr(cr, "raw_meth")
	attr(cr2, "min_dp") = attr(cr, "min_dp")

	cr2
}

reduce_cr = function(cr, expr, txdb, max_gap = 1000) {

	sample_id = attr(cr, "sample_id")
	cor_method = attr(cr, "cor_method")
	raw_meth = attr(cr2, "raw_meth")

	if(raw_meth) warning("only smoothed meth is supported")

	gene = genes(txdb)
	tx_list = transcriptsBy(txdb, by = "gene")

	reduce_cr_by_gene = function(cr, e, max_gap = 1000) {
		if(length(cr) == 0) return(GRanges())
		gr = reduce(cr, min.gapwidth = max_gap, with.revmap = TRUE)
		n = sapply(mcols(gr)[, 1], function(ind) sum(cr$n[ind]))
		mat1 = as.matrix(mcols(cr)[, paste0("mean_meth_", sample_id)])
		mat = lapply(mcols(gr)[, 1], function(ind) colMeans(mat1[ind, , drop =FALSE]))
		mat = do.call("rbind", mat)
		colnames(mat) = paste0("mean_meth_", sample_id)
		corr = apply(mat, 1, function(x) cor(x, e, method = cor_method))
		mcols(gr) = cbind(n, mat, corr)
		gr
	}

	cr_list = split(cr, cr$gene_id)
	i = 0
	res = lapply(names(cr_list), function(gi) {

		qqcat("reducing cr on @{gi}, @{i <<- i+1}/@{length(cr_list)}...\n")
		cr = cr_list[[gi]]
		neg_cr = cr[cr$corr < 0]
		pos_cr = cr[cr$corr > 0]
		gr = c(reduce_cr_by_gene(neg_cr, expr[gi, sample_id], max_gap = max_gap),
			   reduce_cr_by_gene(pos_cr, expr[gi, sample_id], max_gap = max_gap))
		gr = sort(gr)
		gr$gene_id = rep(gi, length(gr))

		## distance to gene tss
		tss = promoters(gene[gi], upstream = 1, downstream = 0)
		gene_tss_dist = as.data.frame(distanceToNearest(gr, tss))[, 3]
		if(as.vector(strand(gene[gi])) == "+") {
			gene_tss_dist = ifelse(end(gr) < start(tss), -gene_tss_dist, gene_tss_dist)
		} else if(as.vector(strand(gene[gi])) == "-") {
			gene_tss_dist = ifelse(start(gr) > end(tss), -gene_tss_dist, gene_tss_dist)
		}
		gr$gene_tss_dist = gene_tss_dist

		## distance to tx tss
		tx = tx_list[[gi]]
		tx = tx[tx$tx_name != gi]
		tss = promoters(tx, upstream = 1, downstream = 0)
		dist = as.data.frame(distanceToNearest(gr, tss))
		tx_tss_dist = dist[, 3]
		if(as.vector(strand(gene[gi])) == "+") {
			tx_tss_dist = ifelse(end(gr) < start(tss[dist[,2]]), -tx_tss_dist, tx_tss_dist)
		} else if(as.vector(strand(gene[gi])) == "-") {
			tx_tss_dist = ifelse(start(gr) > end(tss[dist[,2]]), -tx_tss_dist, tx_tss_dist)
		}
		gr$tx_tss_dist = tx_tss_dist
		gr$nearest_tx_tss = tss[dist[,2]]$tx_name

		gr
	})

	cr2 = do.call("c", res)

	attr(cr2, "factor") = attr(cr, "factor")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "extend") = attr(cr, "extend")
	attr(cr2, "window_size") = attr(cr, "window_size")
	attr(cr2, "sample_id") = attr(cr, "sample_id")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "cov_filter") = attr(cr, "cov_filter")

	cr2
}
