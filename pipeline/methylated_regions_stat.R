 
source("~/project/development/cotools/script/load_all.R")
source("~/project/development/cotools/test/head.R")

setwd("/home/guz/project/analysis/hipo16_new/figure_prepare/")



MR_stat = function(MR_list, name) {
	MR_list = lapply(MR_list[SAMPLE$id], function(gr) {
		gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
		gr[seqnames(gr) %in% paste0("chr", 1:22)]
	})

	pdf(qq("@{name}_statistics.pdf"))
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
	pdf(qq("@{name}_correlation.pdf"), width = 10, height = 5)
	ha = HeatmapAnnotation(subtype = SAMPLE$type, col = list(subtype = SAMPLE_COLOR))
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

load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/LMR_real_volker.RData"); MR_overlap_to_cgi(LMR_list)
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/UMR_filter_volker.RData"); MR_overlap_to_cgi(UMR_list)
