
# == title
# Scatter plot between methylation and expression
#
# == param
# -cr correlated regions from `correlated_regions`
# -expr expression matrix which was used to detect ``cr``
# -gi gene id
# -text_column which column in ``cr`` should be put as annotation text in the plot
# -xlab xlab in the plot
# -ylab ylab in the plot
#
# == details
# Scatterplot for all CRs corresponding to the gene will be made.
#
cr_scatterplot_me = function(cr, expr, gi = NULL, text_column = NULL,
	xlab = "Methylation", ylab = "Expression") {
	
	annotation = attr(cr, "factor")
	annotation_color = attr(cr, "col")
	sample_id = attr(cr, "sample_id")

	if(!is.null(gi)) {
		cr = cr[cr$gene_id %in% gi]
	}

	mmat = mcols(cr)
	mmat = as.matrix(mmat[, paste0("mean_meth_", sample_id)])
	colnames(mmat) = sample_id

	for(i in seq_len(length(cr))) {
		cr_now = cr[i]
		gi = cr_now$gene_id
		chr = as.vector(seqnames(cr_now))

		v = mmat[i, ]
		e = expr[gi, colnames(mmat)]

		if(is.null(annotation)) annotation = rep("unknown", length(e))
		if(is.null(annotation_color)) annotation_color = structure(seq_along(unique(annotation)), names = unique(annotation))

		scatterplot_with_boxplot(v, e, annotation, annotation_color, 
			main = qq("@{gi} cor = @{sprintf('%.2f', cr_now$corr)}\nn_CpG = @{cr_now$n}"),
			xlab = xlab, ylab = ylab,
			text_list = unlist(mcols(cr_now)[, text_column, drop = FALSE]))
		if(i %% 50 == 0) {
			qqcat("@{i}/@{nrow(cr)} CRs finished\n")
		}
	}
}

# == title
# scatterplot with boxplots on both sides
#
# == param
# -x x
# -y y
# -annotation annotations
# -annotation_color colors for annotation
# -main title for the plot
# -xlab xlab
# -ylab ylab
# -xlim xlim
# -ylim ylim
# -text_list additional text
#
scatterplot_with_boxplot = function(x, y, annotation = rep("unknown", length(x)), 
	annotation_color = structure(seq_along(levels(annotation)), names = levels(annotation)),
	main = NULL, xlab = NULL, ylab = NULL, xlim = range(x), ylim = range(y), text_list = NULL) {

	if(!is.factor(annotation)) {
		annotation = factor(annotation)
	}
	xrange = xlim
	yrange = ylim

	layout(rbind(1:2, 3:4), width=c(1, 2), height=c(2, 1))
    par(mar = c(0, 5, 5, 0))
    boxplot(y ~ annotation, ylim = yrange, axes = FALSE, ann = FALSE, col = annotation_color[levels(annotation)])
    box()
    axis(side = 2, cex.axis = 1.5)
    title(ylab = ylab, cex.lab = 1.5)

    par(mar = c(0, 0, 5, 5))
    plot(x, y, xlim = xrange, ylim = yrange, axes = FALSE, ann = FALSE, cex = 1.5, pch = 16, col = annotation_color[annotation])
    box()
    title(main, cex.main = 1)

   # lm.res = lm(y ~ x)
   # abline(lm.res, col = "#E41A1C")

    par(mar = c(5, 5, 0, 0), xpd = NA)
    plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, ann = FALSE)
    legend("bottomleft", legend = levels(annotation), pch = 16, col = annotation_color[levels(annotation)])
    if(length(text_list)) {
	    text_list_name = names(text_list)
	    text_list = as.character(text_list)
	    for(i in seq_along(text_list)) {
	    	if(!is.na(suppressWarnings(as.numeric(text_list[i])))) {
	    		if(nchar(text_list[i]) > 5) {
	    			text_list[i] = sprintf("%.2e", as.numeric(text_list[i]))
	    		}
	    	}
	    }
	    text(0, 1, qq("@{text_list_name} = @{text_list}\n"), adj = c(0, 1))	
	}

    par(mar = c(5, 0, 0, 5))
    boxplot(x ~ annotation, ylim = xrange, horizontal = TRUE, axes = FALSE, ann = FALSE, col = annotation_color[levels(annotation)])
    box()
    axis(side = 1, cex.axis = 1.5)
    #axis(side = 4, at = seq_along(levels(d$cate)), labels=levels(d$cate))
    title(xlab = xlab, cex.lab = 1.5)
}
