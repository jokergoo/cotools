#####################################################################
#
# plot average coverage per sample
#
#####################################################################

library(GetoptLong)
NCORES = 1

GetoptLong(c("NCORES=i", "number of cores"))

setwd(dirname(get_scriptname()))

source("hipo16_head.R")


sum_cov = numeric(length(SAMPLE$id[SAMPLE$WGBS]))
n = sum_cov
coverage_by_chr = list()
for(chr in CHROMOSOME) {
	bsseq = get_bsmooth(chr)
	gc(verbose = FALSE)

	cov = getCoverage(bsseq)[, SAMPLE$id[SAMPLE$WGBS], drop = FALSE]

	qqcat("@{dim(cov)[1]} CpG sites, @{dim(cov)[2]} samples\n")

	t1 = colSums(cov)
	sum_cov = sum_cov + t1
	t2 = apply(cov, 2, function(x) sum(x != 0))
	n = n + t2

	coverage_by_chr[[chr]] = t1/t2
	names(coverage_by_chr[[chr]]) = paste(SAMPLE$type[SAMPLE$WGBS], SAMPLE$id[SAMPLE$WGBS], sep = "-")
}

coverage = sum_cov / n
names(coverage) = paste(SAMPLE$type[SAMPLE$WGBS], SAMPLE$id[SAMPLE$WGBS], sep = "-")


plot_coverage = function(coverage, main = "CpG coverage per sample") {
	par(mar = c(12, 4, 4, 1))
	plot(c(0, length(coverage)), c(0, max(coverage)), axes = FALSE, ann = FALSE, type = "n")
	abline(v = 0, lty = 3, col = "grey")
	for(i in seq_along(coverage)) {
		abline(v = i, lty = 3, col = "grey")
		lines(c(i-1, i), c(coverage[i], coverage[i]), lwd = 2)
	}
	par(las = 3)
	axis(side = 1, at = seq_along(coverage) - 0.5, labels = names(coverage))
	axis(side = 2)
	box()
	par(las = 0)
	title(main, ylab = "Mean CpG coverage")
}

pdf(qq("@{DIR_IMAGE}/average_CpG_coverage.pdf"), width=12, height = 6)
plot_coverage(coverage)
dev.off()

