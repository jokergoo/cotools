library(GetoptLong)

GetoptLong(c("chr=s", ""))

source("~/project/development/cotools/script/load_all.R")
source("~/project/development/cotools/test/head.R")

co_opt$wd = "/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/"

gr = correlated_regions(SAMPLE$id, expr, txdb, chr = chr, factor = SAMPLE$type, col = SAMPLE_COLOR)
saveRDS(gr, file = qq("@{co_opt$wd}/rds/@{chr}_cr.rds"))

gr = correlated_regions(SAMPLE$id, expr, txdb, raw_meth = TRUE, chr = chr, factor = SAMPLE$type, col = SAMPLE_COLOR, cov_cutoff = 5, min_dp = 8)
saveRDS(gr, file = qq("@{co_opt$wd}/rds/@{chr}_cr_raw_meth.rds"))
