library(GetoptLong)

GetoptLong(c("chr=s", ""))

setwd("~/project/development/cotools/R")
source("~/project/development/cotools/test/head.R")
for(f in dir(pattern = "\\.R$")) source(f)

gr = correlated_regions(SAMPLE$id, expr[l, ], txdb, chr = chr, factor = SAMPLE$type, col = SAMPLE_COLOR, raw_meth = TRUE, cov_cutoff = 5, min_dp = 8)

saveRDS(gr, file = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/results/correlated_region/@{chr}_cr_raw_meth.rds"))
