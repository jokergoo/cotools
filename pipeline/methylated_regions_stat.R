 
source("~/project/development/cotools/script/load_all.R")
source("~/project/development/cotools/pipeline/head/head.R")


expr = log2(expr + 1)

MR_stat = function(MR_list, name) {
  MR_list = lapply(MR_list[SAMPLE$id], function(gr) {
    gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
    gr[seqnames(gr) %in% paste0("chr", 1:22)]
  })

  pdf(qq("@{output_dir}/@{name}_statistics.pdf"))
  basic_genomic_regions_stat(MR_list, annotation = SAMPLE$type, annotation_color = SAMPLE_COLOR)
  basic_genomic_regions_stat(MR_list, type = "number", annotation = SAMPLE$type, annotation_color = SAMPLE_COLOR)
  basic_genomic_regions_stat(MR_list, type = "median_width", annotation = SAMPLE$type, annotation_color = SAMPLE_COLOR)
  dev.off()
}

load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/LMR_real_volker.RData"); MR_stat(LMR_list, "LMR")
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/UMR_filter_volker.RData"); MR_stat(UMR_list, "UMR")
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/PMD_filter_volker.RData"); MR_stat(PMD_list, "PMD")
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/DMV_filter_volker.RData"); MR_stat(DMV_list, "DMV")

### enrichment to other ...

MR_corr = function(MR_list, name) {
  MR_list = lapply(MR_list[SAMPLE$id], function(gr) {
    gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
    gr[seqnames(gr) %in% paste0("chr", 1:22)]
  })

  res = genomic_regions_correlation(MR_list, GENOMIC_FEATURE_LIST, chromosome = chromosome, nperm = 2)
  pdf(qq("@{output_dir}/@{name}_correlation.pdf"), width = 10, height = 5)
  ha = HeatmapAnnotation(subtype = SAMPLE$type, age = SAMPLE$age, col = list(subtype = SAMPLE_COLOR, age = AGE_COL_FUN))
  # ht = Heatmap(res$foldChange, top_annotation = ha, column_title = qq("@{name}, fold change, stat: jaccard coefficient"))
  # draw(ht)
  ht = Heatmap(res$stat, top_annotation = ha, column_title = qq("@{name}, jaccard coefficient"), cluster_columns = FALSE)
  draw(ht)
  dev.off()
}

load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/LMR_real_volker.RData"); MR_corr(LMR_list, "LMR")
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/UMR_filter_volker.RData"); MR_corr(UMR_list, "UMR")
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/PMD_filter_volker.RData"); MR_corr(PMD_list, "PMD")
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/DMV_filter_volker.RData"); MR_corr(DMV_list, "DMV")


### LMR/UMR overlap with CGI


MR_overlap_to_cgi = function(MR_list) {
  MR_list = lapply(MR_list[SAMPLE$id], function(gr) {
    gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
    gr[seqnames(gr) %in% paste0("chr", 1:22)]
  })

  sapply(MR_list, function(gr) {
    mtch = as.matrix(findOverlaps(gr, GENOMIC_FEATURE_LIST$cgi))
    length(unique(mtch[, 1]))/length(gr)
  })
}

load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/LMR_real_volker.RData"); LMR_cgi = MR_overlap_to_cgi(LMR_list)
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/UMR_filter_volker.RData"); UMR_cgi = MR_overlap_to_cgi(UMR_list)

ha = HeatmapAnnotation(barplot = function(index) {
  n = length(index)
  pushViewport(viewport(xscale = c(0, n), yscale = c(0, 1)))
  grid.rect(index - 0.5, 0, width = 0.4, height = LMR_cgi[index], default.unit = "native", just = c("left", "bottom"), gp = gpar(fill = "blue"))
  grid.rect(index - 0.5, 0, width = 0.4, height = UMR_cgi[index], default.unit = "native", just = c("right", "bottom"), gp = gpar(fill = "green"))
  grid.yaxis()
  upViewport()
}, subtype = SAMPLE$type, col = list(subtype = SAMPLE_COLOR),
annotation_height = unit(c(9, 0.5), "cm"),
show_legend = FALSE)

mat = matrix(nrow = 0, ncol = length(DMV_genes))
colnames(mat) = SAMPLE$id
ht = Heatmap(mat, top_annotation = ha, column_title = "Percent of LMR/UMR that overlap with CGI")
cm = ColorMapping(level = names(SAMPLE_COLOR), color = SAMPLE_COLOR)
lgd1 = color_mapping_legend(cm, title = "Subtype", plot = FALSE)
cm = ColorMapping(level = c("LMR", "UMR"), color = c(UMR = "green", LMR = "blue"))
lgd2 = color_mapping_legend(cm, title = "Methylated regions", plot = FALSE)
pdf(qq("@{output_dir}/LMR_UMR_overlap_with_CGI.pdf"), width = 10, height = 6)
draw(ht, padding = unit(c(10, 20, 5, 5), "mm"), heatmap_legend_list = list(lgd2, lgd1))
decorate_annotation("barplot", {
  grid.text("percent", x = unit(-15, "mm"), rot = 90, just = "bottom")
})
dev.off()
