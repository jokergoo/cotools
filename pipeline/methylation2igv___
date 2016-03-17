
############################################
#
# save methylation as .tdf files for IGV
#
############################################

library(GetoptLong)
NCORES = 1

GetoptLong(c("NCORES=i", "number of cores"))

setwd(dirname(get_scriptname()))

source("../head.R")


# save raw methylation into text files
t = unique(sample$type)
con = vector("list", length = length(t))
names(con) = t

for(type in t) {
	con[[type]] = file(qq("igv/raw_methylation_@{type}.igv"), "w")
	sample_id = sample$id[sample$type == type]
	writeLines("#track graphType=heatmap midRange=0.5:0.5 midColor=255,255,255 color=255,0,0 altColor=0,0,255 viewLimits=0:1", con = con[[type]])
	writeLines(paste("Chromosome", "Start", "End", "Feature", paste(sample_id, collapse = "\t"), sep = "\t", collapse = "\t"), con = con[[type]])
}

for(chr in chromosome) {
	
	bsseq = get_bsmooth(chr)
	gc(verbose = FALSE)

	# ensure the column order is the same
	m = get_meth(bsseq, NULL, sample$id)
	n = nrow(m)
	site = get_cpg_site(bsseq)

	for(type in unique(sample$type)) {
		qqcat("save raw methylation in @{type}\n")
		m2 = cbind(rep(chr, n), site, site+1, rep(".", n), as.data.frame(m[, sample$id[sample$type == type]]))
		write.table(m2, file = con[[type]], row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
	}

}

for(type in t) {
	close(con[[type]])
	
	system(qq("igvtools toTDF igv/raw_methylation_@{type}.igv igv/raw_methylation_@{type}.tdf hg19"))
	system(qq("rm igv/raw_methylation_@{type}.igv"))
}
