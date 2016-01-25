
mat_list = list()
load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/transition-IDH-MES.RData")
mat_list[["IDH-MES"]] = res$mat
load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/transition-IDH-RTK_I.RData")
mat_list[["IDH-RTK_I"]] = res$mat
load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/transition-IDH-RTK_II.RData")
mat_list[["IDH-RTK_II"]] = res$mat
load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/transition-MES-RTK_I.RData")
mat_list[["MES-RTK_I"]] = res$mat
load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/transition-MES-RTK_II.RData")
mat_list[["MES-RTK_II"]] = res$mat
load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/transition-RTK_I-RTK_II.RData")
mat_list[["RTK_I-RTK_II"]] = res$mat


od = c("E11", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E12", "E1", "E2", "E3")

# lapply(names(mat_list), function(cmp) {
# 	s = strsplit(cmp, "-")[[1]]
# 	pdf(qq("~/project/analysis/hipo16_new/figure_prepare/chromHMM_transition_@{cmp}_all.pdf"), width = 8, height = 8)
# 	chord_diagram(mat_list[[cmp]], subtype1 = s[1], subtype2 = s[2], max_mat = mat_list[[i]])
# 	dev.off()
# })

pdf(qq("~/project/analysis/hipo16_new/figure_prepare/chromHMM_transitions.pdf"), width = 15, height = 10)
par(mfrow = c(2, 3))
mat_list = lapply(mat_list, function(m) {
	m["E11", "E11"] = 0
	m[od, od]
})

i = which.max(sapply(mat_list, sum))

for(cmp in names(mat_list)) {
	s = strsplit(cmp, "-")[[1]]
	chord_diagram(mat_list[[cmp]], subtype1 = s[1], subtype2 = s[2], max_mat = mat_list[[i]])
}


par(mfrow = c(2, 3))
mat_list = lapply(mat_list, function(m) {
	diag(m) = 0
	m[od, od]
})

i = which.max(sapply(mat_list, sum))

for(cmp in names(mat_list)) {
	s = strsplit(cmp, "-")[[1]]
	chord_diagram(mat_list[[cmp]], subtype1 = s[1], subtype2 = s[2], max_mat = mat_list[[i]])
}

dev.off()
