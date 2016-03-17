library(GetoptLong)


source("~/project/development/cotools/script/load_all.R")
source("~/project/development/cotools/pipeline/head/head.R")

GetoptLong(c("chr=s", ""))

gr = correlated_regions(SAMPLE$id, expr, txdb, chr = chr, factor = SAMPLE$type, col = SAMPLE_COLOR)
saveRDS(gr, file = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/rds/@{chr}_cr.rds"))

gr = correlated_regions(SAMPLE$id, expr, txdb, raw_meth = TRUE, chr = chr, factor = SAMPLE$type, col = SAMPLE_COLOR, cov_cutoff = 5, min_dp = 8)
saveRDS(gr, file = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/rds/@{chr}_cr_raw_meth.rds"))
