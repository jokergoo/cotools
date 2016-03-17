library(GetoptLong)
NCORES = 1
GetoptLong(c("NCORES=i", "number of cores"))

setwd(dirname(get_scriptname()))

source("hipo16_head.R")


plot_var_heatmap = function(sample, mat, ml = NULL, main = NULL, test = FALSE) {
	if(is.null(ml)) {
		qqcat("no dist column\n")
	}
	nn = nrow(mat)
	l = apply(mat, 1, function(x) all(!is.na(x)))
	mat = mat[l, , drop = FALSE]
	for(i in seq_along(ml)) {
		ml[[i]] = ml[[i]][l, , drop = FALSE]
	}
	#dist = dist[l]
	qqcat("removed @{length(l)-sum(l)}/@{length(l)} regions\n")
	
	
	# l = apply(mat[, sample$type != "xx", drop = FALSE], 1, function(x) diameter(tapply(x, sample$type[sample$type != "xx"], mean)) > 0.2)
	# mat = mat[l, ,drop = FALSE]; 
	# #dist = dist[l]
	for(i in seq_along(ml)) {
		ml[[i]] = ml[[i]][l, , drop = FALSE]
	}

	if(test) {
		# l = logical(nrow(mat))
		# for(j in seq_len(nrow(mat))) {
		# 	group = factor(sample$type[sample$type != "xx"])
		# 	data = mat[j, sample$id[sample$type != "xx"]]
		# 	data = data + rnorm(length(data), 0, 0.01)
		# 	df = data.frame(data = data, group = group)
		# 	p1 = oneway.test(data ~ group, data = df)$p.value
		# 	l[j] = p1 < 0.01
		# }
		# l[is.na(l)] = FALSE
		l = order(rowSds(mat), decreasing = TRUE)[1:5000]

		mat = mat[l, ,drop = FALSE]; 
		#dist = dist[l]
		for(i in seq_along(ml)) {
			ml[[i]] = ml[[i]][l, , drop = FALSE]
		}
	}
	xx = as.logical(ml[[1]])
	dim(xx) = dim(ml[[1]])
	colnames(ml[[1]]) = paste(colnames(ml[[1]]), colSums(xx), sep = "_")

	colSide = data.frame(subtype = sample$type)
	rownames(colSide) = sample$id
	
	col.fun = generate_col_fun(c(0, 0.5, 1), c("blue", "white", "red"))
	col.fun2 = generate_col_fun(c(0, 1), c("white", "black"))
	col.fun3 = generate_col_fun(c(0, 1000), c("black", "white"))
	mat_list = ml
	

		nx = dim(mat)[1]
		
		if(nx > 10000) {
			l = sample(nx, 10000)
			mat = mat[l, ,drop = FALSE]; 
			for(i in seq_along(mat_list)) {
				mat_list[[i]] = mat_list[[i]][l, , drop = FALSE]
			}
		}
		type = unique(sample$type[sample$type != "xx"])
		mm = NULL
		for(t in type) {
			m1 = mat[, sample$type == t, drop = FALSE]
			hc = hclust(dist(t(m1)))
			m1 = m1[, hc$order, drop = FALSE]
			cn = c(colnames(mm), colnames(m1))
			mm = cbind(mm, m1)
			colnames(mm) = cn
		}
		
		pheatmap_modified(mm, cluster_cols = FALSE,
			mat_list = mat_list, col.fun = col.fun, col.fun_list = list(col.fun2),
			cluster_rows = TRUE, annotation = colSide, border_color = NA, show_rownames=FALSE, annotation_legend = TRUE, main = qq("@{main}, n=@{nx}, percentage=@{sprintf('%.4f', nx/nn)}"))
		
		pheatmap_modified(mm, cluster_cols = TRUE,
			mat_list = mat_list, col.fun = col.fun, col.fun_list = list(col.fun2),
			cluster_rows = TRUE, annotation = colSide, border_color = NA, show_rownames=FALSE, annotation_legend = TRUE, main = qq("@{main}, n=@{nx}, percentage=@{sprintf('%.4f', nx/nn)}"))
					
}


GENOMIC_FEATURE_LIST = list(genomic_window_1kb = "bed/genomic_window_1kb.bed")

gr_list_1 = lapply(GENOMIC_FEATURE_LIST, function(file) {
	qqcat("reading @{file}\n");
	x = read.table(file, sep = "\t", stringsAsFactors = FALSE)
	x = x[x[[1]] %in% CHROMOSOME, , drop = FALSE]
	bedtools("sort -V -k1,1 -k2,2 `x`")
})


gr_list_name = names(gr_list_1)

# scatter or density plot for mean methylation in regions
# NOTE: rows in mat_list are same as in gr, so there may be NAs
# IF group = TRUE, rows in `mat_list` is same as unique(gr$name)
mat_list = get_avg_methyrate_in_regions_mc(SAMPLE$id, gr_list_1, chromosome = CHROMOSOME, ncores = NCORES)
save(mat_list, file = "RData/mat_list_genomic_window_background.RData")

pdf(file = qq("@{DIR_IMAGE}/heatmap_of_methylation_in_genomic_window_1kb.pdf"), width = 12, height = 16)
for(i in seq_along(mat_list)) {
	qqcat("drawing heatmap for @{gr_list_name[i]}\n")

	t = gr_list_1[[i]]
	t = annotate_to_other_region(t, "bed/gene.bed", name = "gene")
	# t = annotate_to_other_region(t, "bed/intergenic.bed", name = "intergenic")
	# t = annotate_to_other_region(t, "bed/exon.bed", name = "exon")
	# t = annotate_to_other_region(t, "bed/intron.bed", name = "intron")
	# t = annotate_to_other_region(t, "bed/tss_2k.bed", name = "tss_2k")
	# t = annotate_to_other_region(t, "bed/ucsc/cpgIslandExt.bed", name = "CGI")
	# t = annotate_to_other_region(t, "bed/repeats_DNA.bed", name = "repeats_DNA")
	# t = annotate_to_other_region(t, "bed/repeats_LINE.bed", name = "repeats_LINE")
	# t = annotate_to_other_region(t, "bed/repeats_Low_complexity.bed", name = "repeats_Low_complexity")
	# t = annotate_to_other_region(t, "bed/repeats_LTR.bed", name = "repeats_LTR")
	# t = annotate_to_other_region(t, "bed/repeats_Simple_repeat.bed", name = "repeats_Simple_repeat")
	# t = annotate_to_other_region(t, "bed/repeats_SINE.bed", name = " repeats_SINE")
	
	plot_var_heatmap(SAMPLE, mat_list[[i]], list(as.matrix(t[, -(1:3)])), main = qq("@{gr_list_name[i]}, most variable"), test = TRUE)
	plot_var_heatmap(SAMPLE[SAMPLE$type != "IDH", ], mat_list[[i]][, SAMPLE$type != "IDH"], list(as.matrix(t[, -(1:3)])), main = qq("CiMP-, @{gr_list_name[i]}, most variable"), test = TRUE)

}
dev.off()

gr_list_1[[1]] = bedtools("bedtools intersect -a `gr_list_1[[1]]` -b bed/ucsc/cpgIslandExt.bed -wa -f 0.25 | uniq")
mat_list = get_avg_methyrate_in_regions_mc(SAMPLE$id, gr_list_1, chromosome = CHROMOSOME, ncores = NCORES)

pdf(file = qq("@{DIR_IMAGE}/heatmap_of_methylation_in_genomic_window_1kb_CiMP.pdf"), width = 12, height = 16)
for(i in seq_along(mat_list)) {
	qqcat("drawing heatmap for @{gr_list_name[i]}\n")

	t = gr_list_1[[i]]
	t = annotate_to_other_region(t, "bed/gene.bed", name = "gene")
	# t = annotate_to_other_region(t, "bed/intergenic.bed", name = "intergenic")
	# t = annotate_to_other_region(t, "bed/exon.bed", name = "exon")
	# t = annotate_to_other_region(t, "bed/intron.bed", name = "intron")
	# t = annotate_to_other_region(t, "bed/tss_2k.bed", name = "tss_2k")
	# t = annotate_to_other_region(t, "bed/ucsc/cpgIslandExt.bed", name = "CGI")
	# t = annotate_to_other_region(t, "bed/repeats_DNA.bed", name = "repeats_DNA")
	# t = annotate_to_other_region(t, "bed/repeats_LINE.bed", name = "repeats_LINE")
	# t = annotate_to_other_region(t, "bed/repeats_Low_complexity.bed", name = "repeats_Low_complexity")
	# t = annotate_to_other_region(t, "bed/repeats_LTR.bed", name = "repeats_LTR")
	# t = annotate_to_other_region(t, "bed/repeats_Simple_repeat.bed", name = "repeats_Simple_repeat")
	# t = annotate_to_other_region(t, "bed/repeats_SINE.bed", name = " repeats_SINE")
	
	plot_var_heatmap(SAMPLE, mat_list[[i]], list(as.matrix(t[, -(1:3)])), main = qq("@{gr_list_name[i]}, most variable"), test = TRUE)
	#plot_var_heatmap(SAMPLE[SAMPLE$type != "IDH", ], mat_list[[i]][, SAMPLE$type != "IDH"], list(as.matrix(t[, -(1:3)])), main = qq("CiMP-, @{gr_list_name[i]}, most variable"), test = TRUE)

}
dev.off()
