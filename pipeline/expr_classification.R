 
suppressPackageStartupMessages(library(GetoptLong))

head = "~/project/development/cotools/pipeline/head/head.R"
GetoptLong(c("head=s", "head R script"))

source("~/project/development/cotools/script/load_all.R")
source(head)


#setwd(output_dir)

gene = genes(txdb)
chr = as.vector(seqnames(gene))
names(chr) = names(gene)
gt = extract_field_from_gencode(gencode_gtf_file, level = "gene", primary_key = "gene_id", field = "gene_type")
gt[gt != "protein_coding"] = "others"


count = as.matrix(expression$count[, SAMPLE$id])
rpkm = as.matrix(expression$rpkm[, SAMPLE$id])
expr = log2(rpkm + 1)
l = apply(count, 1, function(x) sum(x > 0) > length(x)/2)
expr = expr[l, ]

#expr = expr[chr[rownames(expr)] %in% paste0("chr", 1:22), ]

all_mean = apply(expr, 1, mean, trim = 0.1)  # maybe median?
m_10 = quantile(all_mean, 0.1)
m_50 = quantile(all_mean, 0.5)
cv = apply(expr, 1, function(x) {
    if(mean(x) < m_50) return(-Inf)
    if(sd(x) == 0) return(-Inf)
    sd(x)/(mean(x) + m_10)  # maybe mad(x)/median(x)
  })

ind = order(cv, decreasing = TRUE)[1:2000]

mat = expr[ind, ]
set.seed(123)
pdf(NULL)
res = ConsensusClusterPlus(mat, maxK = 6, 
    clusterAlg = "hc", distance = "spearman", reps = 1000, verbose = TRUE)
dev.off()

class = res[[4]]$consensusClass
pdf(qq("@{output_dir}/expression_classification_rnaseq.pdf"), width = 10, height = 10)
ha = HeatmapAnnotation(subtype = SAMPLE$type, age = SAMPLE$age, class = class, 
  col = list(subtype = SAMPLE_COLOR, age = AGE_COL_FUN, class = structure(2:5, names = unique(class))))
ht = Heatmap(mat, col = colorRamp2(quantile(mat, c(0, 0.5, 0.9)), c("blue", "white", "red")), 
  top_annotation = ha, show_row_names = FALSE, cluster_columns = res[[4]]$consensusTree) +
  Heatmap(gt[rownames(mat)], col = c("protein_coding" = "red", "others" = "grey"), 
    show_row_names = FALSE, name = "type")
draw(ht)
decorate_annotation("subtype", {grid.text("subtype", unit(2, "mm") + unit(1, "npc"), just = "left")})
decorate_annotation("age", {grid.text("age", unit(2, "mm") + unit(1, "npc"), just = "left")})
dev.off()
