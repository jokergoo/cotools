

# == title
# correlated regions in extended gene model
#
# == param
# -site CpG sites
# -meth methylation matrix corresponding to ``site``
# -cov coverage
# -expr expression for current gene
# -chr chromosome
# -cov_cutoff cutoff for coverage
# -min_dp minimal number of non-NA values for calculating correlations
# -cor_method method for calcualting correlations
# -window_size how many CpG sites in a window
# -factor subtype
# -max_width maximum width of a window
#
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
			mean(xm[ym >= cov_cutoff], na.rm = TRUE)
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
		if(length(unique(factor)) == 1) factor = NULL
	}

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
	if(nrow(m) == 1) {
		meth_IQR = iqr(m)
	} else {
		meth_IQR = rowIQRs(m, na.rm = TRUE)
	}
	gr = GRanges(seqnames = rep(chr, length(ir)),
		    ranges = ir)
	if(is.null(factor)) {
		df = DataFrame(n = window_size,
		    m, # mean methylation
		    corr = corr,
		    corr_p = corr_p,
		    meth_IQR = meth_IQR)	
	} else {
		df = DataFrame(n = window_size,
		    m, # mean methylation
		    corr = corr,
		    corr_p = corr_p,
		    meth_IQR = meth_IQR,
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
# -cor_method method to calculate correlation
# -factor subtype
# -window_size how many CpGs in a window
# -max_width maximum width of a window
# -raw_meth whether use raw methylation value (unsmoothed)
# -cov_cutoff cutoff for coverage
# -min_dp minimal non-NA values for calculating correlations
# -col color for subtypes
#
# == detail
# based on `correlated_regions_per_gene`
correlated_regions = function(sample_id, expr, txdb, chr, extend = 50000,
	cov_filter = function(x) sum(x > 0, na.rm = TRUE) > length(x)/2,
	cor_method = "spearman", factor = NULL, window_size = 5, max_width = 10000,
	raw_meth = FALSE, cov_cutoff = 3, min_dp = 4, col = NULL) {

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
		if(any(is.na(l))) {
			stop("`cov_filter` generates `NA`, check it.")
		}
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

# == title
# filter correlation regions
# 
# == param
# -chromosome chromosomes
# -template template to find cr files
# -cutoff cutoff of adjusted p-values
# -adj_method method for calculating adjusted p-values
# -meth_diameter_cutoff cutoff for diameters
# -meth_IQR_cutoff cutoff for IQR, if there is no subtype information, IQR is used to remove less variable methylation
# -anova_cutoff cutoff for ANOVA test
#
filter_correlated_regions = function(chromosome = paste0("chr", 1:22), template, 
	cutoff = 0.05, adj_method = "BH", meth_diameter_cutoff = 0.25, meth_IQR_cutoff = 0.25,
	anova_cutoff = 0.05) {

	if(length(cutoff) == 1) cutoff = rep(cutoff, 2)

	cat("calculate fdr...\n")
	corr_p = NULL
	meth_anova = NULL
	meth_diameter = NULL
	meth_IQR = NULL
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
	   	} else {
	   		meth_mat = as.matrix(mcols(cr)[, grep("^mean_meth", colnames(mcols(cr)))])
	   		meth_IQR = c(meth_IQR, rowIQRs(meth_mat, na.rm = TRUE))
	   		# meth_IQR = c(meth_IQR, cr$meth_IQR)
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
		l = ifelse(corr > 0, corr_fdr <= cutoff[1], corr_fdr <= cutoff[2]) & meth_IQR >= meth_IQR_cutoff & !is.na(meth_IQR)
	}

	l = l & !is.na(corr)

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
		if(sum(l[lo])) {
			cr2 = suppressWarnings(c(cr2, cr[l[lo]]))
		}
	}

	attr(cr2, "factor") = attr(cr, "factor")
	attr(cr2, "col") = attr(cr, "col")
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

# == title
# plot that helps to choose a gap
reduce_cr_gap_test = function(cr) {
	neg_cr = cr[cr$corr < 0]
	neg_cr_list = split(as.data.frame(neg_cr), neg_cr$gene_id)
	neg_cr_rainfall = do.call("rbind", lapply(neg_cr_list, rainfallTransform))
	neg_cr_rainfall$ratio = neg_cr_rainfall$dist/(neg_cr_rainfall$end - neg_cr_rainfall$start + 1)

	x = neg_cr_rainfall$dist

	par(mfrow = c(2, 3))
	for(i in c(10, 100, 500, 1000, 2000)) {
		plot(density(log10(x[x >= i])), axes = FALSE, main = qq("neg_cr, dist >= @{i}bp"))
		axis(side = 1, at = 1:10, labels = 10^(1:10))
		axis(side = 2)
		box()
	}

	pos_cr = cr[cr$corr > 0]
	pos_cr_list = split(as.data.frame(pos_cr), pos_cr$gene_id)
	pos_cr_rainfall = do.call("rbind", lapply(pos_cr_list, rainfallTransform))
	pos_cr_rainfall$ratio = pos_cr_rainfall$dist/(pos_cr_rainfall$end - pos_cr_rainfall$start + 1)

	x = pos_cr_rainfall$dist

	par(mfrow = c(2, 3))
	for(i in c(10, 100, 500, 1000, 2000)) {
		plot(density(log10(x[x >= i])), axes = FALSE, main = qq("pos_cr, dist >= @{i}bp"))
		axis(side = 1, at = 1:10, labels = 10^(1:10))
		axis(side = 2)
		box()
	}
}

# == title
# refuce cr regions
#
# == param
# -cr cr 
# -expe expression
# -txdb txdb
# -max_gap maximum gap
# -gap gap
# -mc.cores number of cores
#
# == detail
# pos_CR and neg_CR are reduced separatedly
reduce_cr = function(cr, expr, txdb, max_gap = 1000, gap = 1.0, mc.cores = 1) {

	sample_id = attr(cr, "sample_id")
	cor_method = attr(cr, "cor_method")
	raw_meth = attr(cr, "raw_meth")

	if(raw_meth) warning("only smoothed meth is supported")

	qqcat("extracting gene and tx models.\n")
	gene = genes(txdb)
	tx_list = transcriptsBy(txdb, by = "gene")

	reduce_cr_by_gene = function(cr, e, max_gap = 1000) {
		if(length(cr) == 0) return(GRanges())
		n = cr$n
		meth_mat = as.matrix(mcols(cr)[, paste0("mean_meth_", sample_id)]) * n

		mcols(cr) = cbind(n = n, meth_mat)
		gr = reduce2(cr, max_gap = max_gap, gap = gap)

		n = gr$n
		meth_mat = as.matrix(mcols(gr)[, paste0("mean_meth_", sample_id)]) / n
	
		corr = apply(meth_mat, 1, function(x) cor(x, e, method = cor_method))
		mcols(gr) = cbind(n = n, meth_mat, corr = corr)
		return(gr)
	}

	cr_list = split(cr, cr$gene_id)
	i = 0
	res = mclapply(names(cr_list), function(gi) {

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
	}, mc.cores = mc.cores)

	cr2 = do.call("c", res)

	attr(cr2, "factor") = attr(cr, "factor")
	attr(cr2, "col") = attr(cr, "col")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "extend") = attr(cr, "extend")
	attr(cr2, "window_size") = attr(cr, "window_size")
	attr(cr2, "sample_id") = attr(cr, "sample_id")
	attr(cr2, "cor_method") = attr(cr, "cor_method")
	attr(cr2, "cov_filter") = attr(cr, "cov_filter")

	cr2
}

# == title
# ad subtype specificity columns
#
# == param
# -cr cr
# -cutoff cutoff for ANOVA test
#
# == details
# 1 is defined as the methylation is higher than all other subtypes and the difference is significant.
# -1 is defined as the methylation is lower than all other subtypes and the difference is significant.
# All the others are defined as 0.
add_subtype_specificity = function(cr, cutoff = 0.05, suffix = "_ss") {
	
	factor = attr(cr, "factor")
	if(length(unique(factor)) <= 1) {
		warning("no grouping settings.")
		return(cr)
	}

	level = unique(factor)
	n_level = length(level)
	sample_id = attr(cr, "sample_id")
	subtype_ss = matrix(nrow = length(cr), ncol = n_level)
	colnames(subtype_ss) = level

	meth_mat = mcols(cr)
	meth_mat = meth_mat[, grep("^mean_meth_", colnames(meth_mat))]
	meth_mat = as.matrix(meth_mat)
	counter = set_counter(length(cr))
	for(i in seq_len(length(cr))) {
		x = meth_mat[i, paste0("mean_meth_", sample_id)]
		# pairwise t-test, t-value and p-value
		mat_t = matrix(nrow = n_level, ncol = n_level)
		rownames(mat_t) = level
		colnames(mat_t) = level
		mat_p = mat_t

		for(i1 in 2:n_level) {
			for(i2 in 1:(i1-1)) {
				x1 = x[factor == level[i1]]
				x2 = x[factor == level[i2]]
				test = t.test(x1, x2)
				mat_t[i1, i2] = test$statistic
				mat_p[i1, i2] = test$p.value
				mat_t[i2, i1] = -mat_t[i1, i2]
				mat_p[i2, i1] = mat_p[i1, i2]	
			}
		}

		ss = apply(mat_t, 1, function(x) {
			x = x[!is.na(x)]
			if(all(x > 0)){
				return(1)
			} else if(all(x < 0)) {
				return(-1)
			} else {
				return(0)
			}
		})

		l = apply(mat_p, 1, function(x) {
			x = x[!is.na(x)]
			all(x < cutoff)
		})
		ss[!l] = 0
		subtype_ss[i, ] = ss

		counter()
	}

	colnames(subtype_ss) = paste0(colnames(subtype_ss), suffix)

	mcols(cr) = cbind(as.data.frame(mcols(cr)), subtype_ss)
	return(cr)
}
