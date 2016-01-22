
# files = dir("chromHMM_segmentation", pattern = ".bed$")
# for(f in files) {
# 	cat(f, "\n")
# 	cmd = qq("bedtools intersect -wao -a /icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/hg19.200bp_windows_autosomes.bed -b /icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/chromHMM_segmentation/@{f} | cut -f 7 > /icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/chipseq_transitions/chromHMM_segmentation/@{f}.states")
# 	system(qq("perl ~/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=1:00:00,mem=4G -N bed_intersect_@{f}' '@{cmd}'"))
# }


library(GetoptLong)
GetoptLong(c("s1=s", "",
	     "s2=s", ""))


# 'normal' is a hard coded type
sample_text = 
"id   type
AK015    IDH
AK041    IDH
AK066    IDH
AK068    IDH
AK076    IDH
AK085    IDH
AK102    IDH
AK103    IDH
AK124    IDH
AK199    IDH
AK213    IDH
AK231    IDH
AK030    MES
AK055    MES
AK071    MES
AK072    MES
AK091    MES
AK139    MES
AK153    MES
AK185    MES
AK195    MES
AK227    MES
AK236    MES
AK003    RTK_I
AK017    RTK_I
AK049    RTK_I
AK051    RTK_I
AK142    RTK_I
AK149    RTK_I
AK156    RTK_I
AK165    RTK_I
AK173    RTK_I
AK183    RTK_I
AK203    RTK_I
AK217    RTK_I
AK053    RTK_II
AK074    RTK_II
AK088    RTK_II
AK089    RTK_II
AK098    RTK_II
AK100    RTK_II
AK132    RTK_II
AK158    RTK_II
AK167    RTK_II
AK178    RTK_II
AK188    RTK_II
AK216    RTK_II
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

	by = sum(max_mat)/360
	by = 
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
pdf(qq("chromHMM_transition_@{s1}_@{s2}.pdf"), width = 8, height = 8)
chord_diagram(mat[od, od], subtype1 = s1, subtype2 = s2)
dev.off()
