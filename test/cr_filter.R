library(GetoptLong)

GetoptLong(c("chr=s", ""))

setwd("~/project/development/cotools/R")
source("~/project/development/cotools/test/head.R")
for(f in dir(pattern = "\\.R$")) source(f)

co_opt$wd = "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/"

cr_filtered = filter_correlated_regions(template = qq("@{co_opt$wd}/rds/@{chr}_cr.rds"))
saveRDS(cr_filtered, file = qq("@{co_opt$wd}/rds/cr_filtered_fdr_0.05.rds"))
cr_filtered_raw_meth = filter_correlated_regions(template = qq("@{co_opt$wd}/rds/@{chr}_cr_raw_meth.rds"))
saveRDS(cr_filtered_raw_meth, file = qq("@{co_opt$wd}/rds/cr_filtered_cr_filtered_raw_meth_fdr_0.05.rds"))

cv = rowIQRs(expr)/rowMedians(expr)

cr_filtered = add_subtype_specificity(cr_filtered)
cr_filtered$expr_cv = cv[cr_filtered$gene_id]

cr_filtered_raw_meth = add_subtype_specificity(cr_filtered_raw_meth)
cr_filtered_raw_meth$expr_cv = cv[cr_filtered_raw_meth$gene_id]

