
source("~/project/development/cotools/script/load_all.R")
source("~/project/development/cotools/test/head.R")

co_opt$wd = "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/"

cr_filtered = filter_correlated_regions(template = "@{co_opt$wd}/rds/@{chr}_cr.rds", cutoff = 0.005)
#cr_filtered_raw_meth = filter_correlated_regions(template = "@{co_opt$wd}/rds/@{chr}_cr_raw_meth.rds")

cv = rowIQRs(expr)/rowMedians(expr)

cr_filtered = add_subtype_specificity(cr_filtered)
cr_filtered$expr_cv = cv[cr_filtered$gene_id]

#cr_filtered_raw_meth = add_subtype_specificity(cr_filtered_raw_meth)
#cr_filtered_raw_meth$expr_cv = cv[cr_filtered_raw_meth$gene_id]

saveRDS(cr_filtered, file = qq("@{co_opt$wd}/rds/cr_filtered_fdr_0.005.rds"))
#saveRDS(cr_filtered_raw_meth, file = qq("@{co_opt$wd}/rds/cr_filtered_cr_filtered_raw_meth_fdr_0.05.rds"))


neg_cr = cr_filtered[cr_filtered$corr < 0]
neg_cr_reduced = reduce_cr(neg_cr, expr, txdb, gap = 1)

pos_cr = cr_filtered[cr_filtered$corr > 0]
pos_cr_reduced = reduce_cr(pos_cr, expr, txdb, gap = 2)


