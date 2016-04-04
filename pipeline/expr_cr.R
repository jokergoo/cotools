 
suppressPackageStartupMessages(library(GetoptLong))

head = "~/project/development/cotools/pipeline/head/head.R"
GetoptLong(c("head=s", "head R script"))

source("~/project/development/cotools/script/load_all.R")
source(head)

cr_filtered = readRDS(qq("@{RDS_FOLDER}/cr_filtered_fdr_0.01.rds"))


pdf(qq("@{output_dir}/cr_meth_vs_expr_0.005.pdf"), width = 8, height = 8)


expr = log2(expr + 1)
cr_gi = unique(cr_filtered$gene_id)
base_mean = rowMeans(expr[cr_gi, SAMPLE$id])

q = quantile(base_mean, seq(0, 1, by = 0.1))
plot(density(base_mean), main = "dist of base_mean")
for(x in q) abline(v = x, col = "#AAAAAA")

expr_cv = tapply(cr_filtered$expr_cv, cr_filtered$gene_id, function(x) x[1])

x = base_mean
y = log(expr_cv[cr_gi])

layout(matrix(c(1, 3, 0, 2), 2), width = c(2, 1), height = c(1, 2))
par(mar = c(1, 4, 4, 1))
d = density(x)
plot(d$x, d$y, xlim = range(x), type = "l", axes = FALSE, ann = FALSE)
box()
axis(side = 2)
q1 = quantile(x)
abline(v = q1, col = "#AAAAAA")


par(mar = c(4, 1, 1, 4))
d = density(y)
plot(d$y, d$x, ylim = range(y[y < Inf]), type = "l", axes = FALSE, ann = FALSE)
box()
axis(side = 1)
q2 = quantile(y[y < Inf])
abline(h = q2, col = "#AAAAAA")

par(mar = c(4, 4, 1, 1))
plot(x, y, pch = 16, cex = 0.5, col = "#00000080", xlab = "base_mean", ylab = "log(expr_cv)")
abline(v = q1, col = "#AAAAAA")
abline(h = q2, col = "#AAAAAA")

layout(1)
######################################

cv = tapply(cr_filtered$expr_cv, cr_filtered$gene_id, function(x) x[1])

cv = sort(cv)

op = par(no.readonly = TRUE)

hist(log(cv), main = "histogram of IQR/meidan of expression for each gene")
plot(log(cv), sapply(cv, function(x) sum(cv <= x))/length(cv), type = "l", ylab = "p(X >= x)", main = "CDF")
corr = NULL
for(q in seq(0.5, 99, by = 0.5)) {
	x = cv[round(length(cv)*q/100)]
	gi = names(x)
	x = sprintf("%.2f", x)
	cr_temp = cr_filtered[cr_filtered$gene_id == gi]
	cr_temp = cr_temp[which.max(abs(cr_temp$corr))[1]]
	cr_scatterplot_me(cr_temp, expr, gi = gi, ylab = qq("expression (@{q}% quantile), IQR/median = @{x}"),
		text_column = c("corr_p", "meth_anova", "meth_diameter", "gene_tss_dist"))
	corr = c(corr, cr_temp$corr)
	cat(q, "/100\n")
}

par(op)
par(mfrow = c(2, 1), mar = c(5, 4, 1, 1))
x = cv[round(length(cv)*(seq(0.5, 99, by = 0.5))/100)]
plot(x, abs(corr), pch = 16, col = ifelse(corr > 0, "red", "green"),
	xlab = "IQR/median", ylab = "best correlation for the corresponding CR gene")
plot(x, base_mean[names(x)], pch = 16, xlab = "IQR/median",
	ylab = "base mean of expression")

dev.off()

