library(GetoptLong)

head = "~/project/development/cotools/pipeline/head/head.R"
GetoptLong("head=s", "head R script")

source("~/project/development/cotools/script/load_all.R")
source(head)

cr_filtered = filter_correlated_regions(chromosome = chromosome, template = qq("`RDS_FOLDER`/@{chr}_cr.rds", code.pattern = "`CODE`"), cutoff = 0.01)
#cr_filtered_raw_meth = filter_correlated_regions(template = "@{co_opt$wd}/rds/@{chr}_cr_raw_meth.rds")

cv = rowIQRs(expr)/rowMedians(expr)
names(cv) = rownames(expr)

cr_filtered = add_subtype_specificity(cr_filtered)
cr_filtered$expr_cv = cv[cr_filtered$gene_id]

#cr_filtered_raw_meth = add_subtype_specificity(cr_filtered_raw_meth)
#cr_filtered_raw_meth$expr_cv = cv[cr_filtered_raw_meth$gene_id]

saveRDS(cr_filtered, file = qq("@{RDS_FOLDER}/cr_filtered_fdr_0.01.rds"))
#saveRDS(cr_filtered_raw_meth, file = qq("@{co_opt$wd}/rds/cr_filtered_raw_meth_fdr_0.01.rds"))

cr_filtered = readRDS(qq("@{RDS_FOLDER}/cr_filtered_fdr_0.01.rds"))

reduce_cr_gap_test(cr_filtered)

neg_cr = cr_filtered[cr_filtered$corr < 0]
pos_cr = cr_filtered[cr_filtered$corr > 0]

neg_cr_reduced = reduce_cr(neg_cr, expr, txdb, gap = bp(1000), mc.cores = 4)
pos_cr_reduced = reduce_cr(pos_cr, expr, txdb, gap = bp(1000), mc.cores = 4)

cr_reduced = c(neg_cr_reduced, pos_cr_reduced)
cr_reduced = cr_reduced[order(as.vector(seqnames(cr_reduced)), cr_reduced$gene_id, start(cr_reduced))]

attr(cr_reduced, "factor") = attr(cr_filtered, "factor")
attr(cr_reduced, "sample_id") = attr(cr_filtered, "sample_id")
cr_reduced = add_subtype_specificity(cr_reduced)

saveRDS(cr_reduced, file = qq("@{RDS_FOLDER}/cr_reduced_fdr_0.01_gap_1kb.rds"))
