
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

pdf("general_methylation_distribution_cgi.pdf", width = 10, height = 10)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
  annotation_color = SAMPLE_COLOR, chromosome = chromosome, background = GENOMIC_FEATURE_LIST$cgi, p = 0.01)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
  annotation_color = SAMPLE_COLOR, chromosome = chromosome, by_chr = TRUE, background = GENOMIC_FEATURE_LIST$cgi, p = 0.01)
dev.off()


pdf("general_methylation_distribution_cgi_shores.pdf", width = 10, height = 10)
extended_cgi = GENOMIC_FEATURE_LIST$cgi
start(extended_cgi) = start(extended_cgi) - 2000
end(extended_cgi) = end(extended_cgi) + 2000
shore = setdiff(extended_cgi, GENOMIC_FEATURE_LIST$cgi)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
  annotation_color = SAMPLE_COLOR, chromosome = chromosome, background = shore, p = 0.01)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
  annotation_color = SAMPLE_COLOR, chromosome = chromosome, by_chr = TRUE, background = shore, p = 0.01)
dev.off()

pdf("general_methylation_distribution_neither_cgi_nor_shores.pdf", width = 10, height = 10)
chromInfo = getChromInfoFromUCSC(genome)
chromInfo = chromInfo[chromInfo$chrom %in% chromosome, ]
chromGr = GRanges(chromInfo$chrom, ranges = IRanges(rep(1, nrow(chromInfo)), chromInfo$length))
complement = setdiff(chromGr, extended_cgi)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
  annotation_color = SAMPLE_COLOR, chromosome = chromosome, background = complement)
global_methylation_distribution(sample_id = SAMPLE$id, annotation = SAMPLE$type, 
  annotation_color = SAMPLE_COLOR, chromosome = chromosome, by_chr = TRUE, background = complement)
dev.off()


###########################################################
# differnetial methylation in different genomic features
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
#   background = gr_list_genome[[1]], nperm = 50, stat_fun = genomicCorr.sintersect)


####################
## overlap to CGI
cgi_1kb_window = makeWindows(GENOMIC_FEATURE_LIST$cgi, w = 1000, short.keep = TRUE)
cgi_1kb_window = cgi_1kb_window[width(cgi_1kb_window) > 500]
gr_list_cgi = get_mean_methylation_in_genomic_features(SAMPLE$id, chromosome = chromosome, gf_list = list(cgi_1kb_window = cgi_1kb_window))

shore_1kb_window = makeWindows(shore, w = 1000, short.keep = TRUE)
shore_1kb_window = shore_1kb_window[width(shore_1kb_window) > 500]
gr_list_shore = get_mean_methylation_in_genomic_features(SAMPLE$id, chromosome = chromosome, gf_list = list(shore_1kb_window = shore_1kb_window))

# res2 = genomic_regions_correlation(gr_cgi, GENOMIC_FEATURE_LIST, chromosome = chromosome, 
#   background = gr_list_cgi[[1]], nperm = 50, stat_fun = genomicCorr.sintersect)

genome_1kb_window = makeWindows(setdiff(chromGr, extended_cgi), w = 1000)
genome_1kb_window = genome_1kb_window[width(genome_1kb_window) > 500]
gr_list_genome = get_mean_methylation_in_genomic_features(SAMPLE$id, chromosome = chromosome, gf_list = list(genome_1kb = genome_1kb_window))


save(gr_list_cgi, gr_list_shore, gr_list_genome, file = "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/genome_windows.RData")

pdf("heatmap_diff_methylation_1kb_window.pdf", width = 14, height = 14)
gr_diff_genome = heatmap_diff_methylation_in_genomic_features(gr_list_genome[[1]], annotation = SAMPLE$type, 
    annotation_color = SAMPLE_COLOR, title = "genome 1kb window")
gr_diff_cgi = heatmap_diff_methylation_in_genomic_features(gr_list_cgi[[1]], annotation = SAMPLE$type, 
    annotation_color = SAMPLE_COLOR, title = "cgi 1kb window")
gr_diff_shore = heatmap_diff_methylation_in_genomic_features(gr_list_shore[[1]], annotation = SAMPLE$type, 
    annotation_color = SAMPLE_COLOR, title = "shore 1kb window")
dev.off()


MR_list = list(genome_1kb_window = gr_list_genome[[1]],
             cgi_1kb_window = gr_list_cgi[[1]],
             shore_1kb_window = gr_list_shore[[1]])
res = genomic_regions_correlation(MR_list, GENOMIC_FEATURE_LIST[c("gene", "intron", "intergenic", "repeats_SINE", "repeats_LINE", "dnase", "exon", "tfbs", "tss_2k", "enhancer")], chromosome = chromosome, nperm = 2)
pdf(qq("genome_1kb_window_correlation.pdf"), width = 6, height = 6)
ht = Heatmap(res$stat, name = "jaccard", column_title = qq("jaccard coefficient"), cluster_columns = FALSE)
draw(ht)
dev.off()


################################ only for testing ########################
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

##############################################################
## let's start from here

cgi_1kb_window = makeWindows(extended_cgi, w = 1000, short.keep = TRUE)
cgi_1kb_window = cgi_1kb_window[width(cgi_1kb_window) > 500]
gr_list_cgi = get_mean_methylation_in_genomic_features(SAMPLE$id, chromosome = chromosome, gf_list = list(cgi_1kb_window = cgi_1kb_window))


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

l = anno == "CGI"
cl = get_class(mat[l, ])
class = cl$class[[1]]
row_index = cl$row_index

tb = table(class); i = as.numeric(names(tb)[which.max(tb)])
cl = get_class(mat[l, class == i])
class[class == i] = cl$class[[1]] + 2
row_index = union(row_index, cl$row_index)

tb = table(class); i = as.numeric(names(tb)[which.max(tb)])
cl = get_class(mat[l, class == i])
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


pdf("methylation_classification_wgbs_cgi.pdf", width = 14, height = 14)
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



############################
## classify by methylation


methylation_subtype_classfication = function(gr, p, corr, k, name) {
  gr_1kb_window = makeWindows(gr, w = 1000, short.keep = TRUE)
  gr_1kb_window = gr_1kb_window[width(gr_1kb_window) > 500]
  gr_list_cgi = get_mean_methylation_in_genomic_features(SAMPLE$id, chromosome = chromosome, gf_list = list(gr_1kb_window = gr_1kb_window))

  mat = as.matrix(mcols(gr_list_cgi[[1]])); rownames(mat) = seq_len(nrow(mat))
  mat = mat[, colnames(mat) != "ncpg"]

  set.seed(12345)
  get_class = function(mat, p, corr, k) {
    od = order(rowVars(mat), decreasing = TRUE)[1:round(nrow(mat)*p)]
    mat2 = mat[od, ]

    ct_cgi = cor_cols(t(mat2), abs_cutoff = corr, mc = 2)

    ind = which(ct_cgi[, as.character(corr)] >= k)
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

  cl = get_class(mat, p = p, corr = corr, k = k)
  class = cl$class[[1]]
  row_index = cl$row_index

  tb = table(class); i = as.numeric(names(tb)[which.max(tb)])
  cl = get_class(mat[, class == i], p = p, corr = corr, k = k)
  class[class == i] = cl$class[[1]] + 2
  row_index = union(row_index, cl$row_index)

  tb = table(class); i = as.numeric(names(tb)[which.max(tb)])
  cl = get_class(mat[, class == i], p = p, corr = corr, k = k)
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


  pdf(qq("methylation_classification_wgbs_@{name}.pdf"), width = 14, height = 14)
  draw(ht)
  dev.off()
}


methylation_subtype_classfication(GENOMIC_FEATURE_LIST$cgi, p = 0.2, corr = 0.6, k = 500, name = "cgi")
methylation_subtype_classfication(setdiff(extended_cgi, GENOMIC_FEATURE_LIST$cgi), p = 0.2, corr = 0.6, k = 500, name = "cgi_shores")
methylation_subtype_classfication(setdiff(chromGr, extended_cgi), p = 0.02, corr = 0.8, k = 1000, name = "neither_cgi_nor_shores")
