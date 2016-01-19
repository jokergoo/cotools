
source("~/project/development/cotools/script/load_all.R")
source("~/project/development/cotools/test/head.R")

# chromosome = c("chr21", "chr22")

setwd("/home/guz/project/analysis/hipo16_new/figure_prepare/")

pdf("general_methylation_distribution.pdf", width = 10, height = 10)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
	annotation_color = SAMPLE_COLOR, chromosome = chromosome)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
	annotation_color = SAMPLE_COLOR, chromosome = chromosome, by_chr = TRUE)
dev.off()

pdf("general_methylation_distribution_cgi_and_shores.pdf", width = 10, height = 10)
extended_cgi = GENOMIC_FEATURE_LIST$cgi
start(extended_cgi) = start(extended_cgi) - 2000
end(extended_cgi) = end(extended_cgi) + 2000
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
	annotation_color = SAMPLE_COLOR, chromosome = chromosome, background = extended_cgi, p = 0.01)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
	annotation_color = SAMPLE_COLOR, chromosome = chromosome, by_chr = TRUE, background = extended_cgi, p = 0.01)
dev.off()

pdf("general_methylation_distribution_no_cgi_and_shores.pdf", width = 10, height = 10)
chromInfo = getChromInfoFromUCSC(genome)
chromInfo = chromInfo[chromInfo$chrom %in% chromosome, ]
chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
complement = setdiff(chromGr, extended_cgi)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
	annotation_color = SAMPLE_COLOR, chromosome = chromosome, background = complement)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
	annotation_color = SAMPLE_COLOR, chromosome = chromosome, by_chr = TRUE, background = complement)
dev.off()

# maybe only for QC
#plot_coverage_and_methylation_on_genome("AK100", chromosome = c("chr21", "chr22"))

pdf("multiple_samples_methylation_on_genome.pdf", width = 20, height = 20)
plot_multiple_samples_methylation_on_genome(SAMPLE$id[od], annotation = SAMPLE$type[od], chromosome = chromosome)
dev.off()

## only for QC
wgbs_qcplot("AK100", chromosome = chromosome)


gr_list = get_mean_methylation_in_genomic_features(SAMPLE$id, chromosome = chromosome, gf_list = GENOMIC_FEATURE_LIST)

pdf("heatmap_diff_methylation_in_genomic_features.pdf", width = 16, height = 16)
for(i in seq_along(gr_list)) {
	heatmap_diff_methylation_in_genomic_features(gr_list[[i]], annotation = SAMPLE$type, 
		annotation_color = SAMPLE_COLOR, title = names(gr_list)[i], txdb = txdb, gf_list = GENOMIC_FEATURE_LIST)
}
dev.off()


####################
## split genome into 1kb window
chromInfo = getChromInfoFromUCSC(genome)
chromInfo = chromInfo[chromInfo$chrom %in% chromosome, ]
chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))

genome_1kb_window = makeWindows(chromGr, w = 1000)
gr_list_genome = get_mean_methylation_in_genomic_features(SAMPLE$id, chromosome = chromosome, gf_list = list(genome_1kb = genome_1kb_window))


# res1 = genomic_regions_correlation(gr_genome, GENOMIC_FEATURE_LIST, chromosome = chromosome, 
# 	background = gr_list_genome[[1]], nperm = 50, stat_fun = genomicCorr.sintersect)


####################
## overlap to CGI
cgi_1kb_window = makeWindows(extended_cgi, w = 1000, short.keep = TRUE)
cgi_1kb_window = cgi_1kb_window[width(cgi_1kb_window) > 500]
gr_list_cgi = get_mean_methylation_in_genomic_features(SAMPLE$id, chromosome = chromosome, gf_list = list(cgi_1kb_window = cgi_1kb_window))

# res2 = genomic_regions_correlation(gr_cgi, GENOMIC_FEATURE_LIST, chromosome = chromosome, 
# 	background = gr_list_cgi[[1]], nperm = 50, stat_fun = genomicCorr.sintersect)

genome_1kb_window = makeWindows(setdiff(chromGr, extended_cgi), w = 1000)
gr_list_genome = get_mean_methylation_in_genomic_features(SAMPLE$id, chromosome = chromosome, gf_list = list(genome_1kb = genome_1kb_window))



pdf("heatmap_diff_methylation_1kb_window.pdf", width = 14, height = 14)
gr_genome = heatmap_diff_methylation_in_genomic_features(gr_list_genome[[1]], annotation = SAMPLE$type, 
		annotation_color = SAMPLE_COLOR, title = "genome 1kb window")
gr_cgi = heatmap_diff_methylation_in_genomic_features(gr_list_cgi[[1]], annotation = SAMPLE$type, 
		annotation_color = SAMPLE_COLOR, title = "cgi 1kb window")
dev.off()



## for the genomic 1kb window, just use the top 5000 most variable windows

mat = as.matrix(mcols(gr_list_genome[[1]]))
od = order(rowVars(mat), decreasing = TRUE)[1:round(nrow(mat)*0.05)]
mat = mat[od, ]
gr = gr_list_genome[[1]][od]

ct = cor_cols(t(mat), abs_cutoff = seq(0.5, 0.9, by = 0.025), mc = 2)

nn = do.call("rbind", lapply(seq(1000, 10000, by = 1000), function(x) apply(ct, 2, function(y) sum(y > x))))
matplot(y = nn, x = seq(1000, 10000, by = 1000), type = "l")


ind = which(ct[, "0.8"] >= 1000)
mcols(gr) = NULL

ind = sample(ind, 1000)
mat = mat[ind, ]
gr = gr[ind, ]

gr = annotate_to_genomic_features(gr, list(GENOMIC_FEATURE_LIST$cgi, extended_cgi), name = c("cgi", "shore"))

anno = ifelse(gr$overlap_to_shore > 0, ifelse(gr$overlap_to_cgi > 1 - gr$overlap_to_shore, "CGI", "Shore"), "Others")

ha = HeatmapAnnotation(subtype = SAMPLE$type, col = list("subtype" = SAMPLE_COLOR))
ht = Heatmap(mat, name = "methylation", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), top_annotation = ha, split = anno) + 
	Heatmap(anno, name = "anno", col = c("CGI" = "orange", "Shore" = "green", "Others" = "blue"))


pdf("methylation_classification_wgbs.pdf", width = 14, height = 14)
draw(ht)
dev.off()



mat = log2(expression$rpkm[, SAMPLE$id] + 1)
gt2 = ifelse(gt == "protein_coding", "protein_coding", "others")
od = order(rowVars(mat), decreasing = TRUE)[1:round(nrow(mat)*0.2)]
mat = mat[od, ]
ct = cor_cols(t(mat), abs_cutoff = seq(0.4, 0.8, by = 0.025), method = "spearman")

nn = do.call("rbind", lapply(seq(0, 1000, by = 20), function(x) apply(ct, 2, function(y) sum(y > x))))
matplot(y = nn, x = seq(0, 1000, by = 20), type = "l")
l = ct[, "0.6"] >= 100


l = order(rowVars(mat), decreasing = TRUE)[1:5000]

dist = function(x, y) {
	q1 = quantile(x, c(0.1, 0.9))
	q2 = quantile(y, c(0.1, 0.9))
	l = x > q1[1] & x < q1[2] & y > q2[1] & y < q2[2]
	sqrt(sum((x[l] - y[l])^2))
}
ht = Heatmap(t(scale(t(mat[l, ]))), name = "log2(rpkm + 1)", top_annotation = ha, 
	show_row_names = FALSE, clustering_distance_columns = dist) +
Heatmap(gt2[rownames(mat[l, ])], col = c("protein_coding" = "red", "others" = "grey"), width = unit(5, "mm"), show_row_names = FALSE)
pdf("expression_classification_rnaseq.pdf", width = 14, height = 14)
draw(ht)
dev.off()





## for the genomic 1kb window, just use the top 5000 most variable windows

mat = as.matrix(mcols(gr_list_cgi[[1]])); rownames(mat) = seq_len(nrow(mat))
mat = mat[, colnames(mat) != "ncpg"]
od = order(rowVars(mat), decreasing = TRUE)[1:round(nrow(mat)*0.2)]
mat = mat[od, ]

ct_cgi = cor_cols(t(mat), abs_cutoff = seq(0.5, 0.9, by = 0.025), mc = 2)

ind = which(ct_cgi[, "0.6"] >= 100)
mcols(gr) = NULL
mat_cgi = mat[ind, ]
gr_cgi = gr[ind, ]

# genome
mat = as.matrix(mcols(gr_list_genome[[1]]))
mat = mat[, colnames(mat) != "ncpg"]
od = order(rowVars(mat), decreasing = TRUE)[1:round(nrow(mat)*0.02)]
mat = mat[od, ]
gr = gr_list_genome[[1]][od]

ct_genome = cor_cols(t(mat), abs_cutoff = seq(0.5, 0.9, by = 0.025), mc = 2)

ind = which(ct_genome[, "0.6"] >= 100)
mcols(gr) = NULL
mat_genome = mat[ind, ]
gr_genome = gr[ind, ]






nn = do.call("rbind", lapply(seq(1000, 10000, by = 1000), function(x) apply(ct, 2, function(y) sum(y > x))))
matplot(y = nn, x = seq(1000, 10000, by = 1000), type = "l")






ind = sample(seq_len(nrow(mat_cgi)), 5000)

mat2 = mat_cgi[ind, ]
gr2 = gr_cgi[ind]



library(ConsensusClusterPlus)
mat = as.matrix(mcols(gr_list_cgi[[1]])); rownames(mat) = seq_len(nrow(mat))
mat = mat[, colnames(mat) != "ncpg"]

set.seed(12345)
get_class = function(mat) {
	od = order(rowVars(mat), decreasing = TRUE)[1:round(nrow(mat)*0.2)]
	mat2 = mat[od, ]

	ct_cgi = cor_cols(t(mat2), abs_cutoff = 0.6, mc = 2)

	ind = which(ct_cgi[, "0.6"] >= 500)
	mat2 = mat2[ind, ]
	if(nrow(mat2) > 5000) {
		l = sample(seq_len(nrow(mat2)), 5000)
	} else {
		l = rep(TRUE, nrow(mat2))
	}
	res = ConsensusClusterPlus(mat2[l, ], maxK = 4, 
		clusterAlg = "km", distance = "euclidean", reps = 1000, verbose = TRUE)

	list(class = lapply(res[-1], function(x) x$consensusClass), row_index = rownames(mat2))
}

gr = gr_list_cgi[[1]]
gr2 = annotate_to_genomic_features(gr, list(GENOMIC_FEATURE_LIST$cgi, extended_cgi), name = c("cgi", "shore"))
anno = ifelse(gr2$overlap_to_shore > 0, ifelse(gr2$overlap_to_cgi > 1 - gr2$overlap_to_shore, "CGI", "Shore"), "Others")

cl = get_class(mat[anno == "CGI", ])
class = cl$class[[1]]
row_index = cl$row_index

tb = table(class); i = as.numeric(names(tb)[which.max(tb)])
cl = get_class(mat[anno == "CGI", class == i])
class[class == i] = cl$class[[1]] + 2
row_index = union(row_index, cl$row_index)

tb = table(class); i = as.numeric(names(tb)[which.max(tb)])
cl = get_class(mat[anno == "CGI", class == i])
class[class == i] = cl$class[[1]] + 4
row_index = union(row_index, cl$row_index)

n = length(row_index)
if(length(row_index) > 5000) row_index = sample(row_index, 5000)

gr = gr_list_cgi[[1]]
gr2 = annotate_to_genomic_features(gr[as.numeric(row_index)], list(GENOMIC_FEATURE_LIST$cgi, extended_cgi), name = c("cgi", "shore"))
anno = ifelse(gr2$overlap_to_shore > 0, ifelse(gr2$overlap_to_cgi > 1 - gr2$overlap_to_shore, "CGI", "Shore"), "Others")


m = NULL
type = NULL
class2 = NULL
for(i in unique(class)) {
	dend = as.dendrogram(hclust(dist(t(mat[row_index, class == i]))))
	dend = stats:::reorder.dendrogram(dend, colMeans(mat[row_index, class == i]))
	col_order = order.dendrogram(dend)
	m = cbind(m, mat[row_index, class == i][, col_order])
	type = c(type, SAMPLE$type[class == i][col_order])
	class2 = c(class2, class[class == i][col_order])
}
ha = HeatmapAnnotation(subtype = type, consensusCL = class2, 
	col = list("subtype" = SAMPLE_COLOR, consensusCL = structure(2:5, names = unique(class2))))
ht = Heatmap(m, name = "methylation", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), show_row_names = FALSE,
	top_annotation = ha, split = anno, cluster_columns = FALSE, column_title = qq("@{n} 1kb windows")) + 
	Heatmap(anno, name = "anno", col = c("CGI" = "orange", "Shore" = "green", "Others" = "blue"))


pdf("methylation_classification_wgbs.pdf", width = 14, height = 14)
draw(ht)
dev.off()


#######################################################
## show simply use top 5K most variable rows may contain rows which only contains random noise

ind = order(rowVars(mat), decreasing = TRUE)[1:5000]
ha = HeatmapAnnotation(subtype = SAMPLE$type, 
	col = list("subtype" = SAMPLE_COLOR))
ht = Heatmap(mat[ind, ], name = "methylation", col = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red")), show_row_names = FALSE,
	top_annotation = ha, split = anno, cluster_columns = TRUE)


pdf("methylation_classification_wgbs_top_5k_var.pdf", width = 14, height = 14)
draw(ht)
dev.off()

