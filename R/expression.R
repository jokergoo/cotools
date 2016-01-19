### functions for expression analysis

# == title
# read row count
#
# == param
# -sample_id sample ids
# -template path for count files
#
read_count = function(sample_id, template) {
	
	count = NULL

	sid = sample_id[1]
	qqcat("[@{sid}] reading count data.\n")
	d = read.table(qq(template, envir = list(sid = sid)), sep = "\t", row.names = 1)
	nr = nrow(d)
	d = d[0:4 - nr, , drop = FALSE]
	count = cbind(count, d[[1]])
	gene_id = rownames(d)
	
	for(sid in sample_id[-1]) {
		qqcat("[@{sid}] reading count data.\n")
		d = read.table(qq(template, envir = list(sid = sid)), sep = "\t", row.names = 1)
		nr = nrow(d)
		d = d[0:4 - nr, , drop = FALSE]
		d = d[gene_id, ,drop = FALSE]
		count = cbind(count, d[[1]])
	}
	colnames(count) = sample_id
	rownames(count) = gene_id
	return(count)
}

# == title
# normalize raw count
#
# == param
# -count count matrix, rownames are used to map genes in ``txdb``
# -txdb transcriptome, used if normalized by gene length
# -method rpkm, voom, tpm, tc, md, deseq2, tmm
# -param additional parameters for different normalization methods
# -gene_length_type if normalize to gene length, type
#
# == detail
# useful links:
# http://bib.oxfordjournals.org/content/14/6/671.full.pdf+html
# https://www.biostars.org/p/56919/
# http://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
#
# == value
# normalized expression matrix
normalize_count = function(count, txdb = NULL, method = "rpkm", 
	param = NULL, gene_length_type = c("exon", "gene")) {

	method = tolower(method)[1]
	gene_length_type = match.arg(gene_length_type)[1]
	default_param = list(
		"varianceStabilize" = 1,
		"normalizeTOGeneLength" = 0
	)
	if(!is.null(param)) {
		default_param[names(param)] = param
		default_param$varianceStabilize = as.numeric(default_param$varianceStabilize)
		default_param$normalizeTOGeneLength = as.numeric(default_param$normalizeTOGeneLength)
	}
	param = default_param
	
	if(method == "rpkm" || param$normalizeTOGeneLength) {
		gene_length = get_gene_length(txdb, by = gene_length_type)

		if(length(setdiff(rownames(count), names(gene_length))) != 0) {
			stop("make sure the GTF file that you did counting and imported as TranscriptDb is the same file.\n")
		}
		gene_length = gene_length[rownames(count)]
	}
	qqcat("normalizing raw count by @{method}\n")
	
	if(method == "rpkm") {
		all_count = colSums(count)
		
		rpkm = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(rpkm) = rownames(count)
		colnames(rpkm) = colnames(count)
		for(i in seq_len(nrow(count))) {
			rpkm[i, ] = count[i, ] / gene_length[i]
		}
			
		for(j in seq_len(ncol(count))) {
			rpkm[, j] = rpkm[, j] / all_count[j]
		}
		
		rpkm = 10^9 * rpkm
		return(rpkm)
	} else if(method == "voom") {
		require(limma)
		expr = voom(count, normalize.method = param$normalize.method)$E
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else if(method == "tpm") {
		
		tpm = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(tpm) = rownames(count)
		colnames(tpm) = colnames(count)
		for(i in seq_len(nrow(count))) {
			tpm[i, ] = count[i, ] / gene_length[i]
		}
		
		sum_ratio = colSums(tpm)
		for(j in seq_len(ncol(count))) {
			tpm[, j] = tpm[, j] / sum_ratio[j]
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				tpm[i, ] = tpm[i, ] / gene_length[i] * 1e3
			}
		}
		tpm = tpm * 1e6
		return(tpm)
	} else if(method == "tc") {
		all_count = colSums(count)
		mean_all_count = mean(all_count)
		expr = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(expr) = rownames(count)
		colnames(expr) = colnames(count)
		for(j in seq_len(ncol(count))) {
			expr[, j] = count[, j] / all_count[j] * mean_all_count
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else if(method == "med") {
		median_count = apply(count, 2, function(x) median(x[x > 0]))
		mean_median_count = median(median_count)
		expr = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(expr) = rownames(count)
		colnames(expr) = colnames(count)
		for(j in seq_len(ncol(count))) {
			expr[, j] = count[, j] / median_count[j] * mean_median_count
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else if(method == "deseq2") {
		require("DESeq2")
		df = data.frame(sth = rep("foo", ncol(count)))
		rownames(df) = colnames(count)
		cds = DESeqDataSetFromMatrix(countData = count, colData = df, design = ~ 1)
		cds = estimateSizeFactors(cds)
		cts = counts(cds, normalized=TRUE)
		if(param$varianceStabilize) {
			cds = estimateDispersions(cds)
			vsd = getVarianceStabilizedData(cds)
			expr = as.data.frame(vsd)
		} else {
			expr = as.data.frame(cts)
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else if(method == "tmm") {
		require(edgeR)
		sizes = calcNormFactors(DGEList(counts = count))
		expr = matrix(0, nrow = nrow(count), ncol = ncol(count))
		rownames(expr) = rownames(count)
		colnames(expr) = colnames(count)
		for(j in seq_len(ncol(count))) {
			expr[, j] = count[, j] / sizes$samples$norm.factors[j]
		}
		if(param$normalizeTOGeneLength) {
			for(i in seq_len(nrow(count))) {
				expr[i, ] = expr[i, ] / gene_length[i] * 1e3
			}
		}
		return(expr)
	} else {
		stop(qq("@{method} is not supported.\n"))
	}
}


# == title
# distribution of expression
#
# == param
# -expr expressed matrix
# -log log2 or log10
# -label type of the expression
# -main title
#
# == detail
# it generates two figures, only make plot in [q0, q90]
plot_expr_distribution = function(expr, log = "none", label = NULL, main = NULL) { 
	
	sample_id = colnames(expr)
	if(is.null(label)) label = deparse(substitute(expr))
	
	if(log == "log2") {
		expr = log2(expr + 1)
		label2 = qq("log2(@{label} + 1)")
	} else if(log == "log10") {
		expr = log10(expr + 1)
		label2 = qq("log10(@{label} + 1)")
	} else {
		label2 = label
	}
	
	q90 = quantile(as.matrix(expr), 0.9)
	expr = apply(expr, 2, function(x) {
		x[x > q90] = q90
		x
	})

	op = par(no.readonly = TRUE)
	on.exit(par(op))
	
	par(mar = c(12, 4, 4 ,2), xpd = NA)
	density_heatplot(expr, draw.quantiles = TRUE)
	
	aa = seq(0, 50, by = 2)
	aa = aa[aa <= max(expr)]
	axis(side = 2, at = aa, labels = 2^aa-1)
	par(las = 3)
	axis(side = 1, at = seq_along(sample_id), labels = sample_id)
	mtext(label2, side = 2, line = 2) 
	title(qq("density heatmap for @{label} (< q90) distribution\n@{main}"))

	par(xpd = FALSE)
	den_x = matrix(nrow = 512, ncol = dim(expr)[2])
	den_y = matrix(nrow = 512, ncol = dim(expr)[2])
	for(i in seq_len(dim(expr)[2])) {
		den  = density(expr[, i])
		den_x[, i] = den$x
		den_y[, i] = den$y
	}
	matplot(den_x, den_y, type = "l", xlab = label2, ylab = "density", main = qq("density distribution of @{label}"))
}

