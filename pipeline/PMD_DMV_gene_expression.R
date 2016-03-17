 
source("~/project/development/cotools/script/load_all.R")
source("~/project/development/cotools/pipeline/head/head.R")

# chromosome = c("chr21", "chr22")

#setwd(output_dir)

genes = genes(txdb)
gt = extract_field_from_gencode(gencode_gtf_file, level = "gene", primary_key = "gene_id", field = "gene_type")

genes = genes[intersect(names(gt[gt == "protein_coding"]), rownames(expr))]

### DMV
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/DMV_filter_volker.RData")
DMV_list = lapply(DMV_list[SAMPLE$id], function(gr) {
		gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
		gr[seqnames(gr) %in% paste0("chr", 1:22)]
	})


DMV_genes = lapply(DMV_list, function(gr) {
	p = percentOverlaps(genes, gr)
	p = p[p >= 0.25]
	names(p)
})

non_DMV_genes = lapply(DMV_list, function(gr) {
	p = percentOverlaps(genes, gr)
	p = p[p == 0]
	names(p)
})

DMV_genes_expr = lapply(seq_along(DMV_genes), function(i) expr[DMV_genes[[i]], i])
non_DMV_genes_expr = lapply(seq_along(non_DMV_genes), function(i) expr[non_DMV_genes[[i]], i])

ha = HeatmapAnnotation(DMV_genes = anno_boxplot(DMV_genes_expr, ylim = c(0, 9), outline = FALSE,
										axis = TRUE, gp = gpar(fill = SAMPLE_COLOR[SAMPLE$type])),
					   non_DMV_genes = anno_boxplot(non_DMV_genes_expr, ylim = c(0, 9), outline = FALSE,
					   					axis = TRUE, gp = gpar(fill = SAMPLE_COLOR[SAMPLE$type])),
					   number = anno_barplot(sapply(DMV_genes, length), axis = TRUE,
					   					gp = gpar(fill = SAMPLE_COLOR[SAMPLE$type])),
					   percent = anno_barplot(sapply(DMV_genes, function(x) length(x)/length(genes)), axis = TRUE,
					   					gp = gpar(fill = SAMPLE_COLOR[SAMPLE$type])),
					   age = SAMPLE$age,
					   col = list(age = AGE_COL_FUN),
					   annotation_height = c(8, 8, 8, 8, 1))

mat = matrix(nrow = 0, ncol = length(DMV_genes))
colnames(mat) = SAMPLE$id
ht = Heatmap(mat, top_annotation = ha, top_annotation_height = unit(20, "cm"), 
	column_title = "Expresion of genes that are covered by DMV or not")
cm = ColorMapping(level = names(SAMPLE_COLOR), color = SAMPLE_COLOR)
lgd = color_mapping_legend(cm, title = "Subtype", plot = FALSE)
pdf(qq("@{output_dir}/DMV_genes_expression.pdf"), width = 10, height = 10)
draw(ht, padding = unit(c(10, 20, 5, 5), "mm"), heatmap_legend_list = list(lgd))
decorate_annotation("DMV_genes", {
	grid.text("DMV genes", x = unit(-15, "mm"), rot = 90, just = "bottom")
})
decorate_annotation("non_DMV_genes", {
	grid.text("none DMV genes", x = unit(-15, "mm"), rot = 90, just = "bottom")
})
decorate_annotation("number", {
	grid.text("number", x = unit(-15, "mm"), rot = 90, just = "bottom")
})
decorate_annotation("percent", {
	grid.text("percent", x = unit(-15, "mm"), rot = 90, just = "bottom")
})
decorate_annotation("age", {
	grid.text("age", x = unit(-2, "mm"), just = "right")
})
dev.off()



##  PMD
load("/icgc/dkfzlsdf/analysis/B080/guz/hipo16_gbm/RData/PMD_filter_volker.RData")
PMD_list = lapply(PMD_list[SAMPLE$id], function(gr) {
		gr = GRanges(seqnames = gr[[1]], ranges = IRanges(gr[[2]], gr[[3]]))
		gr[seqnames(gr) %in% paste0("chr", 1:22)]
	})


PMD_genes = lapply(PMD_list, function(gr) {
	p = percentOverlaps(genes, gr)
	p = p[p >= 0.5]
	names(p)
})

non_PMD_genes = lapply(PMD_list, function(gr) {
	p = percentOverlaps(genes, gr)
	p = p[p >= 0.5]
	setdiff(names(genes), names(p))
})

PMD_genes_expr = lapply(seq_along(PMD_genes), function(i) expr[PMD_genes[[i]], i])
non_PMD_genes_expr = lapply(seq_along(non_PMD_genes), function(i) expr[non_PMD_genes[[i]], i])

ha = HeatmapAnnotation(PMD_genes = anno_boxplot(PMD_genes_expr, ylim = c(0, 8), outline = FALSE,
										axis = TRUE, gp = gpar(fill = SAMPLE_COLOR[SAMPLE$type])),
					   non_PMD_genes = anno_boxplot(non_PMD_genes_expr, ylim = c(0, 8), outline = FALSE,
					   					axis = TRUE, gp = gpar(fill = SAMPLE_COLOR[SAMPLE$type])),
					   number = anno_barplot(sapply(PMD_genes, length), axis = TRUE,
					   					gp = gpar(fill = SAMPLE_COLOR[SAMPLE$type])),
					   percent = anno_barplot(sapply(PMD_genes, function(x) length(x)/length(genes)), axis = TRUE,
					   					gp = gpar(fill = SAMPLE_COLOR[SAMPLE$type])),
					   age = SAMPLE$age,
					   col = list(age = AGE_COL_FUN),
					   annotation_height = c(8, 8, 8, 8, 1))

mat = matrix(nrow = 0, ncol = length(PMD_genes))
colnames(mat) = SAMPLE$id
ht = Heatmap(mat, top_annotation = ha, top_annotation_height = unit(20, "cm"), 
	column_title = "Expresion of genes that are covered by PMD or not")
cm = ColorMapping(level = names(SAMPLE_COLOR), color = SAMPLE_COLOR)
lgd = color_mapping_legend(cm, title = "Subtype", plot = FALSE)
pdf(qq("@{output_dir}/PMD_genes_expression.pdf"), width = 10, height = 10)
draw(ht, padding = unit(c(10, 20, 5, 5), "mm"), heatmap_legend_list = list(lgd))
decorate_annotation("PMD_genes", {
	grid.text("PMD genes", x = unit(-15, "mm"), rot = 90, just = "bottom")
})
decorate_annotation("non_PMD_genes", {
	grid.text("none PMD genes", x = unit(-15, "mm"), rot = 90, just = "bottom")
})
decorate_annotation("number", {
	grid.text("number", x = unit(-15, "mm"), rot = 90, just = "bottom")
})
decorate_annotation("percent", {
	grid.text("percent", x = unit(-15, "mm"), rot = 90, just = "bottom")
})
decorate_annotation("age", {
	grid.text("age", x = unit(-2, "mm"), just = "right")
})
dev.off()

