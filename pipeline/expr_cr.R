cv = tapply(cr_filtered$expr_cv, cr_filtered$gene_id, function(x) x[1])

cv = sort(cv)

pdf("/home/guz/project/analysis/hipo16_new/figure_prepare/cr_select_cv_cutoff.pdf", width = 8, height = 8)

hist(log(cv), main = "histogram of IQR/meidan of expression for each gene")
plot(log(cv), sapply(cv, function(x) sum(cv <= x))/length(cv), type = "l", ylab = "p(X >= x)", main = "CDF")
for(q in 1:99) {
	x = cv[round(length(cv)*q/100)]
	gi = names(x)
	x = sprintf("%.2f", x)
	cr_temp = cr_filtered[cr_filtered$gene_id == gi][1]
	cr_scatterplot_me(cr_temp, expr, gi = gi, ylab = qq("expression (@{q}% quantile), IQR/median = @{x}"),
		text_column = c("corr_p", "meth_anova", "meth_diameter", "gene_tss_dist"))
	cat(q, "/100\n")
}

dev.off()
