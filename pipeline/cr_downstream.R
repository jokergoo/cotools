 

source("~/project/development/cotools/script/load_all.R")
source("~/project/development/cotools/pipeline/head/head.R")

cr_filtered = readRDS(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/rds/cr_filtered_fdr_0.01.rds"))
cr_reduced = readRDS(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/rds/cr_reduced_fdr_0.01_gap_1kb.rds"))

attr(cr_filtered, "col") = SAMPLE_COLOR
attr(cr_reduced, "col") = SAMPLE_COLOR

### some per chromosome barplots
neg_cr = cr_filtered[cr_filtered$corr < 0]
pos_cr = cr_filtered[cr_filtered$corr > 0]

n1 = length(neg_cr) * 5
n2 = length(pos_cr) * 5
w1 = sum(width(neg_cr))
w2 = sum(width(pos_cr))

pdf(qq("@{output_dir}/cr_number.pdf"), width = 4, height = 6)
par(mfrow = c(2, 1), mar = c(3, 4, 1, 1))
barplot(c("neg_cr" = n1, "pos_cr" = n2), beside = TRUE, col = c("green", "red"), 
    ylab = "#CpG", axes = FALSE)
axis(side = 2, at = c(0, 2, 4, 6, 8)*100000, labels = c("0K", "200K", "400K", "600K", "800K"))
barplot(c("neg_cr" = w1, "pos_cr" = w2), beside = TRUE, col = c("green", "red"), 
    ylab = "sum(width(cr))", axes = FALSE)
axis(side = 2, at = c(0, 10, 20, 30)*1000000, labels = c("0MB", "10MB", "20MB", "30MB"))
dev.off()

#######
## test
gi = "ENSG00000060982.10"  # BCAT1, neg_cr at two alternative tx start site, 25.05M, 25.1M
cr = readRDS(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/rds/chr12_cr.rds"))
x = cr[cr$gene_id == "ENSG00000060982.10"]
attr(cr_filtered, "col") = attr(cr, "col")

cr = readRDS(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/rds/chr12_cr_raw_meth.rds"))
x2 = cr[cr$gene_id == "ENSG00000060982.10"]

methylation_hooks$set("chr12")
s = start(x)
e = end(x)
cov = sapply(seq_along(s), function(i) {
        ind = extract_sites(s[i], e[i], methylation_hooks$site())
        mean(methylation_hooks$coverage(row_index = ind))
    })

pdf(qq("@{output_dir}/BCAT1_correlation_compare.pdf"))
plot(x$corr, x2$corr, col = ifelse(cov > 5, "red", "black"), pch = ifelse(x$meth_anova < 0.01, 16, 1),
xlab = "corr calculated from smoothed meth", 
    ylab = "corr calculated from raw meth", 
    main = "Compare correlation (to expression) for BCAT1\n(50kb extended, every 5CpG window)")
legend("topleft", pch = c(1, 1, 16), col = c("red", "black", "black"), legend = c("cov > 5", "cov <= 5", "meth_anova < 0.01"))
text(0.2, -0.6, qq("corr_all = @{l = !(is.na(x$corr) | is.na(x2$corr)); sprintf('%.2f', cor(x$corr[l], x2$corr[l]))}\ncorr_diff_meth = @{l = x$meth_anova < 0.01; sprintf('%.2f', cor(x$corr[l], x2$corr[l]))}"), adj = c(0, 0.5))
dev.off()

pdf(qq("@{output_dir}/BCAT1_meth_plot.pdf"), width = 10, height = 8)
compare_meth(cr_filtered, "chr12", 25050000, 25060000, x, x2)
dev.off()


# #cpg and sum_width changing with cutoff
pdf(qq("@{output_dir}/cr_qc.pdf"), width = 8, height = 8)
foo = cr_qc(template = paste0(co_opt$wd, "/rds/@{chr}_cr.rds"))
dev.off()


pdf(qq("@{output_dir}/hilbert_sig.pdf"), width = 8, height = 8)
cr_hilbert(cr = cr_filtered, txdb = txdb, merge_chr = TRUE)
cr_hilbert(cr = cr_filtered, txdb = txdb, merge_chr = FALSE)
dev.off()


pdf(qq("@{output_dir}/hilbert_all.pdf"), width = 18, height = 12)
cr_hilbert(template = paste0(co_opt$wd, "/rds/@{chr}_cr.rds"), txdb = txdb, merge_chr = TRUE)
cr_hilbert(template = paste0(co_opt$wd, "/rds/@{chr}_cr.rds"), txdb = txdb, merge_chr = FALSE)
dev.off()



pdf(qq("@{output_dir}/cr_overlap.pdf"), width = 8, height = 5)
cr_overlap_to_genomic_features(cr_filtered, GENOMIC_FEATURE_LIST)
dev.off()



## how cr enriched at tss
pdf(qq("@{output_dir}/cr_tss.pdf"), width = 10, height = 4)
cr_enriched_at_tss(cr_filtered, txdb)
dev.off()



## test a cgi
pdf("@{output_dir}/LRRC3_meth_plot.pdf", width = 10, height = 8)
compare_meth(cr_filtered, "chr21", 45874000, 45878000)
dev.off()


### scatterplot
pdf(qq("@{output_dir}/cr_scatter.pdf"), width = 8, height = 8)
cr_scatterplot_me(cr_filtered, expr, gi = "ENSG00000060982.10", text_column = c("corr_p", "meth_anova", "meth_diameter", "gene_tss_dist"))
dev.off()




######## part 3: complex heatmaps ###############

make_enriched_heatmap = function(histone_mark, by, on, which) {
	sample = dir("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/", pattern = qq("AK.*@{histone_mark}"))
	sample = gsub(qq("_@{histone_mark}"), "", sample)
	sample = intersect(SAMPLE$id, sample)

	hm_list = list()
	for(i in seq_along(sample)) {
	    cat(sample[i], "\n")
	    gr = get_hm(sample[i], histone_mark)
	    hm_list[[i]] = gr
	}
	names(hm_list) = sample


	if(which == "neg") {
	    cr = cr_filtered[cr_filtered$corr < 0]
	} else {
	    cr = cr_filtered[cr_filtered$corr > 0]
	}

	ht_global_opt(heatmap_column_title_gp = gpar(fontsize = 10))
	pdf(qq("@{output_dir}/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf"), width = 18, height = 10)
	enriched_heatmap_list_on_gene(cr, GENOMIC_FEATURE_LIST$cgi, txdb, log2(expr+1), hm_list, 
	    hm_name = histone_mark, on = on, by = by)
	dev.off()

	system(qq("convert @{output_dir}/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf @{output_dir}/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.png"))

	pdf(qq("@{output_dir}/heatmap_cgi_@{by}_@{which}_@{histone_mark}.pdf"), width = 12, height = 10)
	enriched_heatmap_list_on_tss_cgi(cr, GENOMIC_FEATURE_LIST$cgi, txdb, log2(expr+1), hm_list, 
	    hm_name = histone_mark, by = "gene")
	dev.off()

	system(qq("convert @{output_dir}/heatmap_cgi_@{by}_@{which}_@{histone_mark}.pdf @{output_dir}/heatmap_cgi_@{by}_@{which}_@{histone_mark}.png"))

}

for(hm in c("H3K4me3", "H3K27Ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K9me3")) {
	make_enriched_heatmap(hm, by = "gene", on = "tss", which = "neg")
	make_enriched_heatmap(hm, by = "gene", on = "tss", which = "pos")
	make_enriched_heatmap(hm, by = "tx", on = "tss", which = "neg")
	make_enriched_heatmap(hm, by = "tx", on = "tss", which = "pos")
	make_enriched_heatmap(hm, by = "gene", on = "body", which = "neg")
	make_enriched_heatmap(hm, by = "gene", on = "body", which = "pos")
}


######################## gviz ########################
sample = dir("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/", pattern = qq("AK.*H3K4me3"))
sample = gsub(qq("_H3K4me3"), "", sample)
sample = intersect(SAMPLE$id, sample)

hm_list = list()
for(i in seq_along(sample)) {
    cat(sample[i], "\n")
    gr = get_hm(sample[i], "H3K4me3")
    hm_list[[i]] = gr
}
names(hm_list) = sample

load("/icgc/dkfzlsdf/analysis/B080/guz/gencode/gencode_v19_transcript_merged.RData")
gn = sapply(gene_annotation$gtf, function(x) x$name)
tx_list = transcriptsBy(txdb, by = "gene")
gi = "ENSG00000060982.10"
for(chr in unique(as.character(seqnames(cr_filtered)))) {
	cat(chr, "...\n")
    cr_subset = cr_filtered[seqnames(cr_filtered) == chr]
    for(gi in unique(cr_subset$gene_id)) {
        pdf(qq("@{output_dir}/gviz/gviz_@{chr}_@{gi}_@{gn[gi]}.pdf"), width = 16, height = 12)
        cr_gviz(cr_filtered, gi, expr, txdb, tx_list = tx_list[[gi]]$tx_name, gf_list = GENOMIC_FEATURE_LIST[c("cgi", "tfbs")], 
        	hm_list = hm_list, symbol = gn[gi])
        dev.off()
    }
}
##############################################


## methylation around enhancers


## enrichment to histone mark
CR_corr = function(histone_mark, ...) {

	sample = dir("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/", pattern = qq("AK.*@{histone_mark}"))
	sample = gsub(qq("_@{histone_mark}"), "", sample)
	sample = intersect(SAMPLE$id, sample)

	hm_list = list()
	for(i in seq_along(sample)) {
	    cat(sample[i], "\n")
	    gr = get_hm(sample[i], histone_mark)
	    hm_list[[i]] = gr
	}
	names(hm_list) = sample

	res = genomic_regions_correlation(hm_list, cr_list, chromosome = chromosome, nperm = 0)
	ha = HeatmapAnnotation(subtype = SAMPLE[names(hm_list), "type"], col = list(subtype = SAMPLE_COLOR))
	# ht = Heatmap(res$foldChange, top_annotation = ha, column_title = qq("@{name}, fold change, stat: jaccard coefficient"))
	# draw(ht)
	ht = Heatmap(res$stat, top_annotation = ha, column_title = qq("@{histone_mark}, jaccard coefficient"), 
		cluster_columns = FALSE, ...)
	draw(ht)
}

neg_cr = cr_filtered[cr_filtered$corr < 0]
pos_cr = cr_filtered[cr_filtered$corr > 0]
cr_list = list(neg_cr = neg_cr, pos_cr = pos_cr)
pdf(qq("@{output_dir}/cr_enriched_to_histone_marks.pdf"), width = 12, height = 4)
for(hm in c("H3K4me3", "H3K27Ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K9me3")) {
	CR_corr(hm, cluster_rows = FALSE)
}
dev.off()

label = ifelse(cr_filtered$corr < 0, "neg_cr", "pos_cr")
label2 = paste(cr_filtered$IDH_ss, cr_filtered$MES_ss, cr_filtered$PDGFRA_ss, cr_filtered$RTK_II_ss, sep = "|")
label_merged = paste(label, label2, sep = "_")
cr_list = split(cr_filtered, label_merged)
pdf(qq("@{output_dir}/cr_subtype_enriched_to_histone_marks_subtype_specific_cr.pdf"), width = 12, height = 16)
for(hm in c("H3K4me3", "H3K27Ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K9me3")) {
	CR_corr(hm, split = ifelse(grepl("neg", names(cr_list)), "neg", "pos"), cluster_rows = FALSE)
}
dev.off()


## enrichment to enhancers
# active enhancers K4me1+k27ac - k4me3 (>1kb)
# posied enhancers K4me1 - k4me3 (>1kb). 

CR_corr_to_enhancer = function(type = "active", ...) {

	sample = dir("/icgc/dkfzlsdf/analysis/B080/wangqi/gbm_enhancer/data", pattern = qq("AK.*_enhancer.bed"))
	sample = gsub(qq("_enhancer.bed"), "", sample)
	sample = intersect(SAMPLE$id, sample)

	hm_list = list()
	for(i in seq_along(sample)) {
	    cat(sample[i], "\n")
	    df = read.table(qq("/icgc/dkfzlsdf/analysis/B080/wangqi/gbm_enhancer/data/@{sample[i]}_enhancer.bed"), stringsAsFactors = FALSE)
	    df = df[df[ncol(df)] == type, ]
	    hm_list[[i]] = GRanges(df[[1]], ranges = IRanges(df[[2]], df[[3]]))
	}
	names(hm_list) = sample

	res = genomic_regions_correlation(hm_list, cr_list, chromosome = chromosome, nperm = 0)
	ha = HeatmapAnnotation(subtype = SAMPLE[names(hm_list), "type"], col = list(subtype = SAMPLE_COLOR))
	# ht = Heatmap(res$foldChange, top_annotation = ha, column_title = qq("@{name}, fold change, stat: jaccard coefficient"))
	# draw(ht)
	ht = Heatmap(res$stat, top_annotation = ha, column_title = qq("@{type} enhancers, jaccard coefficient"), 
		cluster_columns = FALSE, ...)
	draw(ht)
}

neg_cr = cr_filtered[cr_filtered$corr < 0]
pos_cr = cr_filtered[cr_filtered$corr > 0]
cr_list = list(neg_cr = neg_cr, pos_cr = pos_cr)
pdf(qq("@{output_dir}/cr_enriched_to_enhancers.pdf"), width = 12, height = 4)
CR_corr_to_enhancer(type = "active", cluster_rows = FALSE)
CR_corr_to_enhancer(type = "poised", cluster_rows = FALSE)
dev.off()


label = ifelse(cr_filtered$corr < 0, "neg_cr", "pos_cr")
label2 = paste(cr_filtered$IDH_ss, cr_filtered$MES_ss, cr_filtered$PDGFRA_ss, cr_filtered$RTK_II_ss, sep = "|")
label_merged = paste(label, label2, sep = "_")
cr_list = split(cr_filtered, label_merged)
pdf(qq("@{output_dir}/cr_subtype_enriched_to_enhancers_subtype_specific_cr.pdf"), width = 12, height = 16)
CR_corr_to_enhancer(type = "active", split = ifelse(grepl("neg", names(cr_list)), "neg", "pos"), cluster_rows = FALSE)
CR_corr_to_enhancer(type = "poised", split = ifelse(grepl("neg", names(cr_list)), "neg", "pos"), cluster_rows = FALSE)
dev.off()

