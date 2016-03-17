
# == title
# Gviz plot for a gene model with tracks
#
# == param
# -gi gene id
# -expe expression matrix
# -txdb transcriptDb object
# -gene_start start of gene
# -gene_end end of the gene
# -tx_list a list of tx
# -species species
# -gf_list a list of gf
# -hm_list a list of hm
# -symbol gene symbol
#
cr_gviz = function(cr, gi, expr, txdb, gene_start = NULL, gene_end = NULL, tx_list = NULL,
	species = "hg19", gf_list = NULL, hm_list = NULL, symbol = NULL) {

	sample_id = attr(cr, "sample_id")
	extend = attr(cr, "extend")
	window_size = attr(cr, "window_size")
	cor_method = attr(cr, "cor_method")
	factor = attr(cr, "factor")
	cov_filter = attr(cr, "cov_filter")
	raw_meth = attr(cr, "raw_meth")
	cov_cutoff = attr(cr, "cov_cutoff")
	min_dp = attr(cr, "min_dp")

	if(is.null(raw_meth)) raw_meth = FALSE
	if(is.null(cov_cutoff)) cov_cutoff = 0
	if(is.null(min_dp)) min_dp = 5
	if(!raw_meth) cov_cutoff = 0
	
	if(!gi %in% cr$gene_id) {
		stop(paste0("cannot find ", gi, "in cr.\n"))
	}

	chr = as.vector(seqnames(cr[cr$gene_id == gi]))[1]
	cr = cr[cr$gene_id == gi]

	if(is.null(methylation_hooks$obj)) methylation_hooks$set(chr)
	if(attr(methylation_hooks$obj, "chr") != chr) methylation_hooks$set(chr)
	
	e = expr[gi, sample_id]

	if(is.null(gene_start) || is.null(gene_end)) {
		gene = genes(txdb)
		gene_start = start(gene[gi])
		gene_end = end(gene[gi])
	}

	gene_start = gene_start - extend
	gene_end = gene_end + extend

	site = methylation_hooks$site()

	gm_site_index = extract_sites(gene_start, gene_end, site, index = TRUE)
	gm_site = site[gm_site_index]
	gm_meth = methylation_hooks$meth(row_index = gm_site_index, col_index = sample_id)
	gm_cov = methylation_hooks$coverage(row_index = gm_site_index, col_index = sample_id)

	if(!is.null(cov_filter)) {
		l = apply(gm_cov, 1, cov_filter)
		gm_site = gm_site[l]
		gm_meth = gm_meth[l, , drop = FALSE]
		gm_cov = gm_cov[l, , drop = FALSE]
	}

	qqcat("rescan on @{gi} to calculate corr in @{window_size} bp cpg window...\n")
	gr = correlated_regions_per_gene(gm_site, gm_meth, gm_cov, e, chr = chr,
		factor = factor, cor_method = cor_method, window_size = window_size,
		cov_cutoff = cov_cutoff, min_dp = min_dp)

	qqcat("add transcripts to givz tracks...\n")
	options(Gviz.ucscUrl="http://genome-euro.ucsc.edu/cgi-bin/")
	trackList = list()
	trackList = pushTrackList(trackList, GenomeAxisTrack())
	trackList = pushTrackList(trackList, IdeogramTrack(genome = species, chromosome = chr))
	grtrack = GeneRegionTrack(txdb, chromosome = chr, start = gene_start, end = gene_end, name="Gene\nmodel", showId = TRUE, rotate.title = TRUE)
	if(!is.null(tx_list)) {
		sg = symbol(grtrack)
		sg[sg %in% tx_list] = paste0("[", sg[sg %in% tx_list], "]")
		symbol(grtrack) = sg
	}
	trackList = pushTrackList(trackList, grtrack)

	## correlation track
	qqcat("add correlation line to givz tracks...\n")
	trackList = pushTrackList(trackList, DataTrack(name = qq("Correlation\nCpG window = @{window_size}"),
								range = gr,
								genome = species,
								data = gr$corr,
								type = c("l", "g"),
								ylim = c(-1, 1)))

	qqcat("add cr to givz tracks...\n")
	pos_cr = cr[cr$corr > 0]
	if(length(pos_cr))
		trackList = pushTrackList(trackList, constructAnnotationTrack(reduce(pos_cr), chr, name = "POS_CR", fill = "red", col = NA, rotate.title = TRUE, start = gene_start, end = gene_end))
	neg_cr = cr[cr$corr < 0]
	if(length(neg_cr))
		trackList = pushTrackList(trackList, constructAnnotationTrack(reduce(neg_cr), chr, name = "NEG_CR", fill = "green", col = NA, rotate.title = TRUE, start = gene_start, end = gene_end))

	qqcat("add methylation to givz tracks...\n")
	meth_mat = as.matrix(mcols(gr)[, paste0("mean_meth_", sample_id)])
	colnames(meth_mat) = NULL
	for(t in unique(factor)) {
		mat = meth_mat[, factor == t]
		trackList = pushTrackList(trackList, DataTrack(name = t,
									start = start(gr),
									end = end(gr),
									chromosome = seqnames(gr),
									genome = species,
									data = t(mat),
									type = "heatmap",
									showSampleNames = FALSE,
									gradient = c("blue", "white", "red"),
									size = 3,
									col = NA))
	}
	
	### CpG density per 1000bp
	qqcat("add cpg density to givz tracks...\n")
	segment = seq(gm_site[1], gm_site[length(gm_site)], by = 500)
	start = segment[-length(segment)]
	end = segment[-1]-1
	num = sapply(seq_along(start), function(i) sum(gm_site >= start[i] & gm_site <= end[i]))
	trackList = pushTrackList(trackList, DataTrack(name = "#CpG\nper 500bp",
		                            start = start,
		                            end = end,
		                            chromosome = rep(chr, length(start)),
									genome = species,
									data = num,
									col = "black",
									type = "l",
									rotate.title = TRUE,
									size = 2))
	
	qqcat("add other genomic features to givz tracks...\n")
	gf_name = names(gf_list)
	for(i in seq_along(gf_list)) {
		trackList = pushTrackList(trackList, constructAnnotationTrack(gf_list[[i]], chr, name = gf_name[i], rotate.title = TRUE, start = gene_start, end = gene_end))
	}

	if(!is.null(hm_list)) {
		hm_list2 = lapply(hm_list, function(gr) {
			gr = gr[seqnames(gr) == chr]
			l = start(gr) > gene_start & end(gr) < gene_end
			gr[l]
		})

		hm_merged = GRanges()
		seqinfo(hm_merged) = seqinfo(hm_list[[1]])
		for(i in seq_along(hm_list2)) {
			if(length(hm_list2[[i]])) hm_merged = c(hm_merged, hm_list2[[i]])
		}
		if(length(hm_merged) > 0) {
			segments = as(coverage(hm_merged), "GRanges")
			# covert to matrix
			hm_mat = matrix(0, nrow = length(hm_list), ncol = length(segments))
			rownames(hm_mat) = names(hm_list)
			for(i in seq_along(hm_list2)) {
				mtch = as.matrix(findOverlaps(segments, hm_list2[[i]]))
				hm_mat[i, mtch[, 1]] = hm_list2[[i]][mtch[, 2]]$density
			}
			segments = c(segments, GRanges(chr, ranges = IRanges(gene_end - 100, gene_end), score = 0))
			
			for(t in unique(factor)) {
				mat = hm_mat[rownames(hm_mat) %in% sample_id[factor == t], , drop = FALSE]
				mat = cbind(mat, rep(0, nrow(mat)))
				mat[1, ncol(mat)] = max(hm_mat)
				trackList = pushTrackList(trackList, DataTrack(name = t,
											start = start(segments),
											end = end(segments),
											chromosome = seqnames(segments),
											genome = species,
											data = mat,
											type = "heatmap",
											showSampleNames = TRUE,
											gradient = c("white", "purple"),
											size = 2,
											col = NA))
			}
		}
	}

	qqcat("draw gviz plot...\n")
	plotTracks(trackList, from = gene_start, to = gene_end, chromosome = chr, main = paste(gi, symbol))

	#grid.text(paste(gf_name, collapse = "\n"), x = unit(4, "mm"), y = unit(4, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 8))
		
	rm(list = ls())
	gc()

	return(invisible(NULL))
}

pushTrackList = function(trackList, track) {
	if(!is.null(track)) {
		trackList[[length(trackList) + 1]] = track
	}
	return(trackList)
}

constructAnnotationTrack = function(gr, chr, name = NULL, genome = "hg19", start = 0, end = Inf, ...) {
	gr2 = gr[seqnames(gr) == chr]
	gr2 = gr2[end(gr2) > start & start(gr2) < end]

	if(length(gr2)) {

		AnnotationTrack(name = name,
		                start = start(gr2),
		                end = end(gr2),
		                chromosome = seqnames(gr2),
		                genome = genome, 
		                stacking = "dense",
		                showTitle = TRUE, 
		                height = unit(5, "mm"),
		                ...)
	} else {
		NULL
	}
}

