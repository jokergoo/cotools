
# files = dir("chromHMM_segmentation", pattern = ".bed$")
# for(f in files) {
# 	cat(f, "\n")
# 	cmd = qq("bedtools intersect -wao -a /icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/hg19.200bp_windows_autosomes.bed -b /icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/chromHMM_segmentation/@{f} | cut -f 7 > /icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/chromHMM_segmentation/@{f}.states")
# 	system(qq("perl ~/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=1:00:00,mem=4G -N bed_intersect_@{f}' '@{cmd}'"))
# }

library(GetoptLong)
GetoptLong(c("s1=s", "",
	     "s2=s", ""))
library(reshape2)
library(circlize)

sample_text = 
"id   type   age
AK015    IDH  28
AK041    IDH  47
AK066    IDH  34
AK068    IDH  59
AK076    IDH  47
AK085    IDH  38
AK102    IDH  49
AK103    IDH  40
AK124    IDH  32
AK199    IDH  33
AK213    IDH  48
AK231    IDH  56
AK030    MES  66
AK055    MES  55
AK071    MES  49
AK091    MES  58
AK139    MES  46
AK153    MES  63
AK185    MES  36
AK188    MES  55
AK195    MES  58
AK227    MES  55
AK236    MES  44
AK003    RTK_I  44
AK017    RTK_I  19
AK049    RTK_I  73
AK051    RTK_I  55
AK142    RTK_I  49
AK149    RTK_I  73
AK156    RTK_I  63
AK165    RTK_I  36
AK173    RTK_I  59
AK183    RTK_I  49
AK203    RTK_I  20
AK217    RTK_I  61
AK053    RTK_II  62
AK072    RTK_II  51
AK074    RTK_II  40
AK088    RTK_II  63
AK089    RTK_II  62
AK098    RTK_II  36
AK100    RTK_II  38
AK132    RTK_II  60
AK158    RTK_II  59
AK167    RTK_II  69
AK178    RTK_II  57
AK216    RTK_II  61
#normal_occipital     normal
#normal_parietal     normal
#normal_temporal     normal
#normal_frontal     normal
"

SAMPLE = read.table(textConnection(sample_text), header = TRUE, stringsAsFactors = FALSE)
rownames(SAMPLE) = SAMPLE$id


setwd("/home/guz/project/analysis/hipo16_new/figure_prepare/")

library(GetoptLong)

files = dir("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/chromHMM_segmentation", pattern = "states$")
names(files) = gsub("(AK\\d\\d\\d).*$", "\\1", files)

states = lapply(files, function(f) {
	scan(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/chromHMM_segmentation/@{f}"), what = "character")
})
names(states) = gsub("(AK\\d\\d\\d).*$", "\\1", files)
states = as.data.frame(states)

make_matrix = function(subtype1, subtype2, min = 3) {

	s1 = intersect(SAMPLE$id[SAMPLE$type == subtype1], colnames(states))
	s2 = intersect(SAMPLE$id[SAMPLE$type == subtype2], colnames(states))

	m1 = states[, s1]
	m2 = states[, s2]

	states1 = apply(m1, 1, function(x) {
		t = table(x)
		if(max(t) >= min) {
			names(t[which.max(t)])[1]
		} else {
			NA
		}
	})

	states2 = apply(m2, 1, function(x) {
		t = table(x)
		if(max(t) >= min) {
			names(t[which.max(t)])[1]
		} else {
			NA
		}
	})

	l = !is.na(states1) & !is.na(states2)
	index = which(l)
	x = paste(states1[l], states2[l])
	tb = table(x)
	gr_index_list = lapply(names(tb), function(t) {
			index[x == t]
		})
	names(gr_index_list) = names(tb)

	x = gsub("^(E\\d+).*$", "\\1", names(tb))
	y = gsub("^E\\d+ (E\\d+).*$", "\\1", names(tb))
	df = data.frame(x = x, y = y, v = tb)
	mat = dcast(x~y, data = df)
	rownames(mat) = mat[, 1]
	mat = as.matrix(mat[, -1])
	mat[is.na(mat)] = 0

	return(list(mat = mat, gr_index_list = gr_index_list))
}


res = make_matrix(s1, s2); save(res, file = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/transition-@{s1}-@{s2}.RData"))

chord_diagram = function(mat, max_mat = mat, subtype1, subtype2) {
	grid.col = c("black", "red", "red", "red", "orange", "orange", "purple", "purple", "purple", "blue", "blue", "blue")

	par(xpd = NA)

	tb = mat
	rownames(tb) = 0:11
	colnames(tb) = 0:11
	rownames(tb) = paste0("R", 0:11)
	colnames(tb) = paste0("C", 0:11)
	colmat = rep(grid.col, 12); dim(colmat) = dim(tb); colmat = rgb(t(col2rgb(colmat)), max = 255)
	qati = quantile(tb, 0.7)
	colmat[tb > qati] = paste0(colmat[tb > qati], "A0")
	colmat[tb <= qati] = paste0(colmat[tb <= qati], "20")
	dim(colmat) = dim(tb);

	de = 360 - (360 - 20 - 30) * sum(mat)/sum(max_mat) - 30
	circos.par(start.degree = -de/4, gap.degree = c(rep(1, 11), de/2, rep(1, 11), de/2))
	gcd = rep(grid.col, 2); names(gcd) = c(rownames(tb), colnames(tb))
	chordDiagram(tb, grid.col = gcd, col = colmat,
		directional = TRUE, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1))

	text(-1, 1, qq("state\nin @{subtype2}"), adj = c(0, 1), cex = 1.5)
	text(1, -1, qq("state\nin @{subtype1}"), adj = c(1, 0), cex = 1.5)

	for(sn in get.all.sector.index()) {
		if(abs(get.cell.meta.data("cell.start.degree", sector.index = sn) - get.cell.meta.data("cell.end.degree", sector.index = sn)) > 3) {
			circos.text(get.cell.meta.data("xcenter", sector.index = sn, track.index = 2), get.cell.meta.data("ycenter", sector.index = sn, track.index = 2), 
				gsub("C|R", "", sn), col = "white", font = 2, sector.index = sn, track.index = 2, adj = c(0.5, 0.5), niceFacing = TRUE)
			circos.axis(sector.index = sn, track.index = 2, major.tick.percentage = 0.2, labels.away.percentage = 0.2, labels.cex = 0.5)
		}
	}

	circos.clear()
}

mat = res$mat
diag(mat) = 0
od = c("E11", "E4", "E5", "E6", "E7", "E8", "E9", "E10", "E12", "E1", "E2", "E3")
pdf(qq("~/project/analysis/hipo16_new/figure_prepare/chromHMM_transition_@{s1}_@{s2}.pdf"), width = 8, height = 8)
chord_diagram(mat[od, od], subtype1 = s1, subtype2 = s2)
dev.off()

if(!interactive()) q(save = "no")

#################### read all matrix


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
