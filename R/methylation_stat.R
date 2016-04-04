### functions related to QC
## the functions make related plots and output some tables

### mat_density_plot for coverage and methylation, add mean and median lines on the heatmap


# == title
# Basic qc plot for bisulfite sequencing data
#
# == param
# -sample_id a vector of sample ids
# -chromosome a vector of chromosomes
#
# == detail
# For each sample id, it will produce five plots:
#
# 1. mean/median CpG coverage per chromosome
# 2. histogram of CpG coverage
# 3. methylation per chromosome 
# 4. histogram of methylation
# 5. mean Methylation for each CpG coverage 
#
wgbs_qcplot = function(sample_id, chromosome = paste0("chr", 1:22)) {

	# coverage and methylation per chromosome
	data = rep(list(list(cov = NULL, meth = NULL, strand = NULL)), length(sample_id))
	names(data) = sample_id
	
	for(chr in chromosome) {
			
		methylation_hooks$set(chr)
		for(sid in sample_id) {

			cv = methylation_hooks$coverage(col_index = sid)
			l = which(cv != 0)
			if(length(l) > 10000) {
				l = sort(sample(l, 10000))
			}
			cv = cv[l]

			mh = as.vector(methylation_hooks$meth(row_index = l, col_index = sid))
			
			strd = rep("*", length(mh))
			
			data[[sid]]$cov[[chr]] = cv
			data[[sid]]$meth[[chr]] = mh
			data[[sid]]$strand[[chr]] = strd
		}
	}

	for(sid in sample_id) {

		cov = data[[sid]]$cov
		meth = data[[sid]]$meth
		strand = data[[sid]]$strand

		par(mfrow = c(2, 3))
		
		# mean coverage per chromosome
		cpg_coverage_mean = sapply(cov, mean)
		cpg_coverage_median = sapply(cov, median)
		plot(c(0, length(cpg_coverage_mean)), c(0, max(c(cpg_coverage_mean, cpg_coverage_median))), axes = FALSE, ann = FALSE, type="n")
		for(i in seq_along(cpg_coverage_mean)) {
			abline(v = i, lty = 3, col = "grey")
			lines(c(i-1, i), c(cpg_coverage_mean[i], cpg_coverage_mean[i]), lwd = 2)
			lines(c(i-1, i), c(cpg_coverage_median[i], cpg_coverage_median[i]), lwd = 2, col = "red")
		}
		abline(v = 0, lty = 3, col = "grey")
		par(las = 3)
		axis(side = 1, at = seq_along(cpg_coverage_mean)-0.5, labels = names(cpg_coverage_mean))
		axis(side = 2)
		box()
		par(las = 0)
		title(main = qq("Coverage per chromosome (@{sid})"), ylab = "mean and median CpG coverage")
		legend("bottomleft", lty=1, col = c("black", "red"), legend = c("mean", "median"))

		
		# coverage distribution
		if(all(unique(unlist(strand)) %in% c("+", "-"))) {
			x = unlist(cov)
			y = unlist(strand)
			x1 = x[y == "+"]
			x2 = x[y == "-"]
			ta = table(x)
			ta1 = table(x1)
			ta2 = table(x2)
			xlim = range(c(as.numeric(names(ta)), as.numeric(names(ta1)), as.numeric(names(ta2))))
			ylim = range(ta, ta1, ta2)
			plot(as.numeric(names(ta)), ta, xlim = xlim, ylim = ylim, main = qq("histogram of CpG coverage (@{sid})"), log = "x", axes = FALSE, type = "h", ylab = "", xlab="CpG coverage")
			axis(side = 1)
			breaks = seq(0, max(ta)/sum(ta), by = 0.02)
			axis(side = 2, at = breaks*sum(ta), labels = breaks)
			box()
			par(new = TRUE)
			plot(as.numeric(names(ta1))+0.2, ta1, xlim = xlim, ylim = ylim, type = "h", col = "red", log = "x", axes = FALSE, ann = FALSE)
			par(new = TRUE)
			plot(as.numeric(names(ta2))+0.4, ta2, xlim = xlim, ylim = ylim, type = "h", col = "blue", log = "x", axes = FALSE, ann = FALSE)
			#axis(side = 2)
			legend("topright", lty = 1, col = c("black", "red", "blue"), legend = c("strand *", "strand +", "strand -"))
			par(new = FALSE)
		} else {
			ta = table(unlist(cov))
			plot(as.numeric(names(ta)), ta, main = qq("histogram of CpG coverage (@{sid})"), log = "x", axes = FALSE, type = "h", ylab = "", xlab="CpG coverage")
			axis(side = 1)
			breaks = seq(0, max(ta)/sum(ta), by = 0.02)
			axis(side = 2, at = breaks*sum(ta), labels = breaks)
			box()
		}
		
		# mean methylation per chromosome
		cpg_methyrate_mean = sapply(meth, mean)
		cpg_methyrate_median = sapply(meth, median)
		plot(c(0, length(cpg_methyrate_mean)), c(0, 1), axes = FALSE, ann = FALSE, type = "n")
		for(i in seq_along(cpg_methyrate_mean)) {
			abline(v = i, lty = 3, col = "grey")
			lines(c(i-1, i), c(cpg_methyrate_mean[i], cpg_methyrate_mean[i]), lwd = 2)
			lines(c(i-1, i), c(cpg_methyrate_median[i], cpg_methyrate_median[i]), lwd = 2, col = "red")
		}
		abline(v = 0, lty = 3, col = "grey")
		par(las = 3)
		axis(side = 1, at = seq_along(cpg_methyrate_mean) - 0.5, labels = names(cpg_methyrate_mean))
		axis(side = 2)
		box()
		par(las = 0)
		title(main = qq("methylation per chromosome (@{sid})"), ylab = "mean and median methylation")
		legend("bottomleft", lty=1, col = c("black", "red"), legend = c("mean", "median"))

		# distribution of methylation on all chromosomes
		hist(unlist(meth), main = qq("histogram of methylation (@{sid})"), xlab = "methylation")

		
		# methylation to coverage
		coverage2methyrate = tapply(unlist(meth), unlist(cov), mean)
		plot(as.numeric(names(coverage2methyrate)), coverage2methyrate, ylim = c(0, 1), pch=16, log = "x", cex = 0.8, xlab = "CpG coverage", ylab = "mean methylation", main = qq("Mean Methylation for each CpG coverage (@{sid})"))
		coverage2methyrate = tapply(unlist(meth), unlist(cov), median)
		points(as.numeric(names(coverage2methyrate)), coverage2methyrate, pch=16, cex = 0.8, col = "red")
		legend("bottomleft", pch = 16, col = c("black", "red"), legend = c("mean", "median"))

		par(mfrow = c(1, 1))
	}
	
	return2(data, invisible = TRUE)
}

# == title
# coverage and methylation for one sample
#
# == param
# -sid a single sample id
# -chromosome chromosome
# -species species
# -nw number of windows
# -... pass to `gtrellis::initialize_layout`
#
# == details
# The whole genome is segented by ``nw`` windows and mean methylation and mean CpG coverage
# are visualized as two tracks.
#
plot_coverage_and_methylation_on_genome = function(sid, chromosome = paste0("chr", 1:22), 
	species = "hg19", nw = 10000, ...) {

	w = round(read.chromInfo(species = species)$chr.len["chr1"]/nw)
	flag = 0
	col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"), transparency = 0.8)
	for(chr in chromosome) {
		methylation_hooks$set(chr)

		meth = methylation_hooks$meth(col_index = sid)[,1]
		cov = methylation_hooks$coverage(col_index = sid)[,1]; cov = log10(cov+1)
		site = methylation_hooks$site()
		gr = methylation_hooks$GRanges()
		chr_len = read.chromInfo(species = species)$chr.len[chr]
		chr_gr = GRanges(seqname = chr, ranges = IRanges(1, chr_len))
		chr_window = makeWindows(chr_gr, w = w)
		mtch = as.matrix(findOverlaps(chr_window, gr))

		gr2 = chr_window[unique(mtch[,1])]
		meth = tapply(mtch[,2], mtch[,1], function(i) mean(meth[i]))
		cov = tapply(mtch[,2], mtch[,1], function(i) mean(cov[i]))
		
		if(flag == 0) {
			gtrellis_layout(category = chromosome, species = species, 
				n_track = 2, track_ylim = c(0, quantile(cov, 0.95), 0, 1), 
				track_ylab = c("log10(coverage)", "methylation"),
				add_name_track = TRUE, add_ideogram_track = TRUE, ...)
			flag = 1
		}

		add_track(gr2, track = 2, cate = chr, panel.fun = function(gr) {
			x = (start(gr) + end(gr))/2
			y = cov
			grid.points(x, y, pch = ".", gp = gpar(col = "#FF000010"))
		})
		add_track(gr2, track = 3, cate = chr, panel.fun = function(gr) {
			x = (start(gr) + end(gr))/2
			y = meth
			grid.points(x, y, pch = ".", gp = gpar(col = col_fun(y)))
		})
	}
}




# == title
# methylation for more than one samples
#
# == param
# -sample_id a vector of sample ids
# -annotation annotation of samples (e.g. subtypes)
# -chromosome chromosome
# -species species
# -nw number of windows
# -... pass to `gtrellis::initialize_layout`
#
# == details
# The whole genome is segented by ``nw`` windows
#
plot_multiple_samples_methylation_on_genome = function(sample_id, annotation, 
	chromosome = paste0("chr", 1:22), species = "hg19", nw = 1000, ...) {
	
	col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"), transparency = 0.8)
	type = unique(annotation)
		
	n = table(annotation)[type]
	ty = numeric(2*length(n))
	ty[seq_along(n)*2-1] = 0.5
	ty[seq_along(n)*2] = n + 0.5

	gtrellis_layout(category = chromosome, species = species, 
		n_track = length(type), track_ylab = type, track_ylim = ty, track_height = n,
		add_name_track = TRUE, add_ideogram_track = TRUE, ...)

	w = round(read.chromInfo(species = species)$chr.len["chr1"]/nw)
	for(chr in chromosome) {
		methylation_hooks$set(chr)

		meth = methylation_hooks$meth(col_index = sample_id)
		site = methylation_hooks$site()
		gr = methylation_hooks$GRanges()
		chr_len = read.chromInfo(species = species)$chr.len[chr]
		chr_gr = GRanges(seqname = chr, ranges = IRanges(1, chr_len))
		chr_window = makeWindows(chr_gr, w = w)
		mtch = as.matrix(findOverlaps(chr_window, gr))

		gr2 = chr_window[unique(mtch[,1])]
		meth = tapply(mtch[,2], mtch[,1], function(i) colMeans(meth[i, , drop = FALSE]))

		meth = do.call("rbind", meth)

		for(i in seq_along(type)) {
			sid = sample_id[annotation == type[i]]
			m = meth[, sid, drop = FALSE]
			add_track(gr2, track = i+1, cate = chr, panel.fun = function(gr) {
				x = (start(gr2) + end(gr2))/2
				for(i in seq_along(sid)) {
					y = rep(i, length(x)) + (runif(length(x))-0.5)*0.8
					grid.points(x, y, pch = ".", gp = gpar(col = col_fun(m[, i])))
				}
			})
		}
	}
}

.mat_dist = function(mat, anno, col, ha = NULL, title = NULL, ...) {

	od = order(factor(anno, levels = unique(anno), ordered = TRUE), colMedians(mat, na.rm = TRUE))
	# mat = mat[, od]

	if(is.null(ha)) ha = HeatmapAnnotation(df = data.frame(type = anno), col = list(type = col))

	densityHeatmap(mat, anno = ha, title = title, column_order = od, ...)
	for(an in sapply(ha@anno_list, function(x) x@name)) {
		decorate_annotation(an, {
			grid.text(an, x = unit(-2, "mm"), just = "right")
		})
	}
	## MDS plot
	loc = cmdscale(dist2(t(mat), pairwise_fun = function(x, y) {l = is.na(x) | is.na(y); x = x[!l]; y = y[!l]; sqrt(sum((x-y)^2))}))
	plot(loc[, 1], loc[, 2], pch = 16, cex = 1, col = col[anno], main = qq("MDS:@{title}"), xlab = "dimension 1", ylab = "dimension 2")
	legend("bottomleft", pch = 16, legend = names(col), col = col)
}


# == title
# Global methylation distribution
# 
# == param
# -sample_id a vector of sample ids
# -annotation subtype information
# -annotation_color color for subtypes
# -ha additional annotation can be specified as a `ComplexHeatmap::HeatmapAnnotation` object
# -chromosome chromosomes
# -by_chr whether make the plot by chromosome
# -max_cov maximum coverage (used to get rid of extremely high coverage which affect visualization of CpG coverage distribution)
# -background background to look into
# -p probability to randomly sample CpG sites
#
# == details
# It visualize distribution of methylation valus and CpG coverages through heatmaps.
#
global_methylation_distribution = function(sample_id, annotation, 
	annotation_color = structure(seq_along(unique(annotation)), names = unique(annotation)),
	ha = NULL, chromosome = paste0("chr", 1:22), by_chr = FALSE, max_cov = 100,
	background = NULL, p = 0.001) {

	annotation_color = annotation_color[intersect(names(annotation_color), unique(annotation))]
	
	###############################################
	# distribution of global methylation
	meth_mat = NULL
	cov_mat = NULL
	for(chr in chromosome) {

		methylation_hooks$set(chr)

		meth_gr = methylation_hooks$GRanges()
		if(!is.null(background)) {
			mtch = as.matrix(findOverlaps(meth_gr, background))
			ind = unique(mtch[, 1])
		} else {
			ind = seq_len(length(meth_gr))
		}
		
		nr = length(ind)
		ind = ind[sample(c(FALSE, TRUE), nr, replace = TRUE, p = c(1-p, p))]
		qqcat("random sampled @{length(ind)} sites from @{nr} sites on @{chr} (with p = @{p})\n")
		mm = methylation_hooks$meth(row_index = ind, col_index = sample_id)
		meth_mat = rbind(meth_mat, mm)
		
		cm = methylation_hooks$coverage(row_index = ind, col_index = sample_id)
		cov_mat = rbind(cov_mat, cm)

		if(by_chr) {
			cm[cm == 0] = NA
			cm[cm > max_cov] = NA
			.mat_dist(mm, anno = annotation, col = annotation_color, ha = ha, title = qq("methylation:@{chr}"), range = c(0, 1), ylab = "methylation")
			.mat_dist(log10(cm), anno = annotation, col = annotation_color, ha = ha, title = qq("coverage:@{chr}"), range = c(0, Inf), ylab = qq("log10(CpG coverage, 1~@{max_cov})"))
		}
	}

	if(!by_chr) {
		nr = nrow(meth_mat)
		if(nr > 100000) {
			meth_mat = meth_mat[sample(nr, 100000), ]
			cov_mat = cov_mat[sample(nr, 100000), ]
		}
		cov_mat[cov_mat == 0] = NA
		cov_mat[cov_mat > max_cov] = NA

		.mat_dist(meth_mat, anno = annotation, col = annotation_color, ha = ha, title = "methylation", range = c(0, 1), ylab = "methylation")
		.mat_dist(log10(cov_mat), anno = annotation, col = annotation_color, ha = ha, title = "coverage", range = c(0, Inf), ylab = qq("log10(CpG coverage, 1~@{max_cov})"))
		
		od = order(factor(annotation, levels = unique(annotation), ordered = TRUE), colMedians(meth_mat, na.rm = TRUE))
		return(od)
	}

}
