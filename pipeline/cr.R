
library(GetoptLong)

head = "~/project/development/cotools/pipeline/head/head.R"
GetoptLong(c("chrxx=s", "which chromosome",
	         "head=s", "head R script"))

source("~/project/development/cotools/script/load_all.R")
source(head)



gr = correlated_regions(SAMPLE$id, expr, txdb, chr = chrxx, factor = SAMPLE$type, col = SAMPLE_COLOR)
saveRDS(gr, file = qq("@{RDS_FOLDER}/@{chrxx}_cr.rds"))

# gr = correlated_regions(SAMPLE$id, expr, txdb, raw_meth = TRUE, chr = chrxx, factor = SAMPLE$type, col = SAMPLE_COLOR, cov_cutoff = 5, min_dp = 8)
# saveRDS(gr, file = qq("@{RDS_FOLDER}/@{chrxx}_cr_raw_meth.rds"))
