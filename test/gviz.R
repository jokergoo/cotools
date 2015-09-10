library(GetoptLong)

GetoptLong(c("chr=s", ""))

setwd("~/project/development/cotools/R")
for(f in dir(pattern = "\\.R$")) source(f)
source("~/project/development/cotools/test/head.R")

cr_filtered = readRDS("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/correlated_region/cr_filtered_fdr_0.05.rds")

cr_filtered = cr_filtered[seqnames(cr_filtered) == chr]

txdb2 = subset_txdb(txdb, chr)

tx_list = transcriptsBy(txdb2, by = "gene")
gene = genes(txdb2)

i = 0
for(gi in unique(cr_filtered$gene_id)) {
    qqcat("@{i = i+1}/@{length(unique(cr_filtered$gene_id))} genes finished.\n")
    gg = gene[gi]
    pdf(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/correlated_region/image/gviz_@{chr}_@{gi}.pdf"), width = 16, height = 12)
    cr_gviz(cr_filtered, gi, expr, txdb2, gene_start = start(gg)[1], gene_end = end(gg)[1], tx_list = tx_list[[gi]]$tx_name, gf_list = GENOMIC_FEATURE_LIST[c("cgi", "dnase", "tfbs")])
    dev.off()
}
