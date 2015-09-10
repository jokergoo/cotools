library(GetoptLong)



histone_mark = "H3K4me3"
by = "gene"
on = "tss"
which = "neg"

GetoptLong(c("histone_mark=s", "",
	         "by=s", "gene|tx",
	         "on=s", "tss|body",
	         "which=s", "neg|pos"))

setwd("~/project/development/cotools/R")
for(f in dir(pattern = "\\.R$")) source(f)
source("~/project/development/cotools/test/head.R")


cr_filtered = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/correlated_region/cr_filtered_fdr_0.05.rds")

############ histone marks ###############

sample = dir("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/", pattern = qq("AK.*@{histone_mark}"))
sample = gsub(qq("_@{histone_mark}"), "", sample)

hm_list = list()
for(i in seq_along(sample)) {
    cat(sample[i], "\n")
    df = read.table(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/@{sample[i]}_@{histone_mark}/sicer/@{sample[i]}_@{histone_mark}.mkdup-W200-G200-islands-summary-FDR0.01"))
    gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), density = df[[7]])
    hm_list[[i]] = gr
}
names(hm_list) = sample

if(which == "neg") {
	cr = cr_filtered[cr_filtered$corr < 0]
} else {
	cr = cr_filtered[cr_filtered$corr > 0]
}

ht_global_opt(heatmap_column_title_gp = gpar(fontsize = 10))
pdf(qq("heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf"), width = 18, height = 10)
enriched_heatmap_list_on_gene(cr, GENOMIC_FEATURE_LIST$cgi, txdb, log2(expr+1), hm_list, 
    hm_name = histone_mark, on = on, by = by)
dev.off()

system(qq("convert heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.png"))
