library(GetoptLong)

GetoptLong(c("chr=s", ""))

source("~/project/development/cotools/script/load_all.R")
source("~/project/development/cotools/test/head.R")

co_opt$wd = "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/"

##### initialize the folder 


########## part 1: generate cr and cr_filtered ##############
qid = NULL
for(chr in chromosome) {

    qid = c(qid, qsub(name = qq("cr_@{chr}"), 
                      options = qq("-l walltime=10:00:00,mem=10G"), 
                      private = "chr", 
                      share = TRUE, 
                      {
                        # !!!!! expr should only contain protein_coding genes!!!!!
                        gr = correlated_regions(SAMPLE$id, expr, txdb, chr = chr, factor = SAMPLE$type, col = SAMPLE_COLOR)
                        saveRDS(gr, file = qq("@{co_opt$wd}/rds/@{chr}_cr.rds"))
                      })
            )

    qid = c(qid, qsub(name = qq("cr_raw_meth_@{chr}"), 
                      options = qq("-l walltime=10:00:00,mem=10G"), 
                      private = "chr", 
                      share = TRUE, 
                      {
                        gr = correlated_regions(SAMPLE$id, expr, txdb, raw_meth = TRUE, chr = chr, factor = SAMPLE$type, col = SAMPLE_COLOR, cov_cutoff = 5, min_dp = 8)
                        saveRDS(gr, file = qq("@{co_opt$wd}/rds/@{chr}_cr_raw_meth.rds"))
                      })
            )
}


##### compare smoothed and raw methylation

qid = qsub(name = "cr_filter",
           options = "-l walltime=10:00:00,mem=10G", 
           dependency = qid, 
           {
            cr_filtered = filter_correlated_regions(template = qq("@{co_opt$wd}/rds/@{chr}_cr.rds"))
            saveRDS(cr_filtered, file = qq("@{co_opt$wd}/rds/cr_filtered_fdr_0.05.rds"))
            cr_filtered_raw_meth = filter_correlated_regions(template = qq("@{co_opt$wd}/rds/@{chr}_cr_raw_meth.rds"))
            saveRDS(cr_filtered_raw_meth, file = qq("@{co_opt$wd}/rds/cr_filtered_cr_filtered_raw_meth_fdr_0.05.rds"))
           }
)

cv = rowIQRs(expr)/rowMedians(expr)


cr_filtered = add_subtype_specificity(cr_filtered)
cr_filtered$expr_cv = cv[cr_filtered$gene_id]

cr_filtered_raw_meth = add_subtype_specificity(cr_filtered_raw_meth)
cr_filtered_raw_meth$expr_cv = cv[cr_filtered_raw_meth$gene_id]



############## part 2: downstream analysis ################

cr_filtered = readRDS(qq("@{co_opt$wd}/rds/cr_filtered_fdr_0.05.rds"))

### some per chromosome barplots
neg_cr = cr_filtered[cr_filtered$corr < 0]
pos_cr = cr_filtered[cr_filtered$corr > 0]

n1 = length(neg_cr) * 5
n2 = length(pos_cr) * 5
w1 = sum(width(neg_cr))
w2 = sum(width(pos_cr))

pdf(qq("@{co_opt$wd}/image/cr_number.pdf"), width = 4, height = 6)
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
cr = readRDS(qq("@{co_opt$wd}/rds/chr12_cr.rds"))
x = cr[cr$gene_id == "ENSG00000060982.10"]

cr = readRDS(qq("@{co_opt$wd}/rds/chr12_cr_raw_meth.rds"))
x2 = cr[cr$gene_id == "ENSG00000060982.10"]

methylation_hooks$set("chr12")
s = start(x)
e = end(x)
cov = sapply(seq_along(s), function(i) {
        ind = extract_sites(s[i], e[i], methylation_hooks$site())
        mean(methylation_hooks$coverage(row_index = ind))
    })

pdf(qq("@{co_opt$wd}/image/BCAT1_correlation_compare.pdf"))
plot(x$corr, x2$corr, col = ifelse(cov > 5, "red", "black"), pch = ifelse(x$meth_anova < 0.01, 16, 1),
xlab = "corr calculated from smoothed meth", 
    ylab = "corr calculated from raw meth", 
    main = "Compare correlation (to expression) for BCAT1\n(50kb extended, every 5CpG window)")
legend("topleft", pch = c(1, 1, 16), col = c("red", "black", "black"), legend = c("cov > 5", "cov <= 5", "meth_anova < 0.01"))
text(0.2, -0.6, qq("corr_all = @{sprintf('%.2f', cor(x$corr, x2$corr))}\ncorr_diff_meth = @{l = x$meth_anova < 0.01; sprintf('%.2f', cor(x$corr[l], x2$corr[l]))}"), adj = c(0, 0.5))
dev.off()

pdf(qq("@{co_opt$wd}/image/BCAT1_meth_plot.pdf"), width = 10, height = 8)
compare_meth(cr_filtered, "chr12", 25050000, 25060000, x, x2)
dev.off()


# #cpg and sum_width changing with cutoff
qsub(name = "cr_qc", options = "-l walltime=3:00:00,mem=5G", 
    share = TRUE, {
        pdf(qq("@{co_opt$wd}/image/cr_qc.pdf"), width = 8, height = 16)
        foo = cr_qc(template = paste0(co_opt$wd, "/rds/@{chr}_cr.rds"))
        dev.off()
})


qsub(name = "hilbert_sig", options = "-l walltime=3:00:00,mem=5G", 
    share = TRUE, {
        pdf(qq("@{co_opt$wd}/image/hilbert_sig.pdf"), width = 12, height = 12)
        cr_hilbert(cr = cr_filtered, txdb = txdb, merge_chr = TRUE)
        dev.off()
})

qsub(name = "hilbert_all", options = "-l walltime=3:00:00,mem=5G", 
    share = TRUE, {
        pdf(qq("@{co_opt$wd}/image/hilbert_all.pdf"), width = 18, height = 12)
        cr_hilbert(template = paste0(co_opt$wd, "/rds/@{chr}_cr.rds"), txdb = txdb)
        dev.off()
})

qsub(name = "cr_annotate", options = "-l walltime=3:00:00,mem=5G", 
    share = TRUE, {
        pdf(qq("@{co_opt$wd}/image/cr_annotate.pdf"), width = 8, height = 5)
        cr_correlated_to_genomic_features(cr_filtered, GENOMIC_FEATURE_LIST)
        dev.off()
})


## how cr enriched at tss
qsub(name = "cr_tss", options = "-l walltime=3:00:00,mem=5G", 
    share = TRUE, {
        pdf(qq("@{co_opt$wd}/image/cr_tss.pdf"), width = 10, height = 4)
        cr_enriched_at_tss(cr_filtered, txdb)
        dev.off()
})


## test a cgi
pdf("LRRC3_meth_plot.pdf", width = 10, height = 8)
compare_meth("chr21", 45874000, 45878000)
dev.off()


### scatterplot
qsub(name = "cr_scatterplot_me", options = "-l walltime=3:00:00,mem=5G", 
    share = TRUE, {
        pdf(qq("@{co_opt$wd}/image/cr_scatter.pdf"), width = 8, height = 8)
        cr_scatterplot_me(cr_filtered, expr, gi = "ENSG00000060982.10", text_column = c("corr_p", "meth_anova", "meth_diameter", "gene_tss_dist"))
        dev.off()
})



######## part 3: complex heatmaps ###############


get_hm = function(sid, histone_mark) {
    df = read.table(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/@{sid}_@{histone_mark}/sicer/@{sid}_@{histone_mark}.mkdup-W200-G200-islands-summary-FDR0.01"))
    gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), density = df[[7]])
}


histone_mark = "H3K4me3"
by = "gene"
on = "tss"
which = "neg"

# load histone marks

sample = dir("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/", pattern = qq("AK.*@{histone_mark}"))
sample = gsub(qq("_@{histone_mark}"), "", sample)

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
pdf(qq("@{co_opt$wd}/image/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf"), width = 18, height = 10)
enriched_heatmap_list_on_gene(cr, GENOMIC_FEATURE_LIST$cgi, txdb, log2(expr+1), hm_list, 
    hm_name = histone_mark, on = on, by = by)
dev.off()

system(qq("convert @{co_opt$wd}/image/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf @{co_opt$wd}/image/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.png"))

pdf(qq("@{co_opt$wd}/image/heatmap_cgi_@{by}_@{on}_@{which}_@{histone_mark}.pdf"), width = 12, height = 10)
enriched_heatmap_list_on_tss_cgi(cr, GENOMIC_FEATURE_LIST$cgi, txdb, log2(expr+1), hm_list, 
    hm_name = histone_mark, by = "gene")
dev.off()

system(qq("convert @{co_opt$wd}/image/heatmap_cgi_@{by}_@{on}_@{which}_@{histone_mark}.pdf @{co_opt$wd}/image/heatmap_cgi_@{by}_@{on}_@{which}_@{histone_mark}.png"))

### gviz
tx_list = transcriptsBy(txdb, by = "gene")
gi = "ENSG00000060982.10"
for(chr in unique(cr_filtered$gi)) {
    qsub(name = qq("cr_gviz_@{chr}"), options = "-l walltime=3:00:00,mem=5G", 
        share = TRUE, {
            cr_subset = cr_filtered[seqnames(cr_filtered) == chr]
            for(gi in unique(cr_subset$gi)) {
                pdf(qq("@{co_opt$wd}/image/gviz/gviz_@{chr}_@{gi}.pdf"), width = 16, height = 12)
                gviz_cr(cr_filtered, gi, expr, txdb, tx_list = tx_list[[gi]]$tx_name, gf_list = GENOMIC_FEATURE_LIST[c("cgi", "tfbs")], hm_list = hm_list)
                dev.off()
            }
        })
}

