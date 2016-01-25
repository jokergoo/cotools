 
chromatin_states_transition_chord_diagram = function(mat, max_mat = mat, cate1, cate2, ...) {

	par(xpd = NA)

	tb = mat
	n = nrow(mat)
	rownames(tb) = 1:n
	colnames(tb) = 1:n
	rownames(tb) = paste0("R", 1:n)
	colnames(tb) = paste0("C", 1:n)
	colmat = rep(grid.col, n); dim(colmat) = dim(tb); colmat = rgb(t(col2rgb(colmat)), max = 255)
	qati = quantile(tb, 0.7)
	colmat[tb > qati] = paste0(colmat[tb > qati], "A0")
	colmat[tb <= qati] = paste0(colmat[tb <= qati], "20")
	dim(colmat) = dim(tb);

	de = 360 - (360 - 20 - 30) * sum(mat)/sum(max_mat) - 30
	circos.par(start.degree = -de/4, gap.degree = c(rep(1, n-1), de/2, rep(1, n-1), de/2))
	gcd = rep(grid.col, 2); names(gcd) = c(rownames(tb), colnames(tb))
	chordDiagram(tb, col = colmat,
		directional = TRUE, annotationTrack = "grid", preAllocateTracks = list(track.height = 0.1), ...)

	text(-1, 1, qq("states\nin @{cate1}"), adj = c(0, 1), cex = 1.5)
	text(1, -1, qq("states\nin @{cate2}"), adj = c(1, 0), cex = 1.5)

	for(sn in get.all.sector.index()) {
		if(abs(get.cell.meta.data("cell.start.degree", sector.index = sn) - get.cell.meta.data("cell.end.degree", sector.index = sn)) > 3) {
			circos.text(get.cell.meta.data("xcenter", sector.index = sn, track.index = 2), get.cell.meta.data("ycenter", sector.index = sn, track.index = 2), 
				gsub("C|R", "", sn), col = "white", font = 2, sector.index = sn, track.index = 2, adj = c(0.5, 0.5), niceFacing = TRUE)
			circos.axis(sector.index = sn, track.index = 2, major.tick.percentage = 0.2, labels.away.percentage = 0.2, labels.cex = 0.5)
		}
	}

	circos.clear()
}
