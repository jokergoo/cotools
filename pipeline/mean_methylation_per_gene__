
library(GetoptLong)
NCORES = 1

GetoptLong(c("NCORES=i", "number of cores"))

setwd(dirname(get_scriptname()))

source("hipo16_head.R")

load(qq("@{DIR_RDATA}/gencode_broad_merged.RData"))

bed = data.frame(chr = sapply(gencode, function(x) x$chr),
	             start = sapply(gencode, function(x) x$start),
	             end = sapply(gencode, function(x) x$end),
	             name = names(gencode),
	             stringsAsFactors = FALSE)

mat_mean = matrix(nrow = nrow(bed), ncol = nrow(SAMPLE))
rownames(mat_mean) = names(gencode)
colnames(mat_mean) = SAMPLE$id
mat_sd = mat_mean
mat_mean_upstream = mat_mean
mat_sd_upstream = mat_mean
mat_mean_downstream = mat_mean
mat_sd_downstream = mat_mean


for(chr in CHROMOSOME) {
	bsseq = get_bsmooth(chr)
	gc(verbose = FALSE)
	
	qqcat("building index for CpG site in @{chr}...\n") 
	site = get_cpg_site(bsseq)

	site.index = build.position.index(site, by = 100000)

	bed_subset = bed[bed[[1]] == chr, ]
	for(i in seq_len(nrow(bed_subset))) {
		CpG_site_index = get_CpG_site_index(bed_subset[i, 2], bed_subset[i, 3], site, site.index)
		if(length(CpG_site_index)) {
			m = get_meth(bsseq, CpG_site_index, SAMPLE$id)
			mat_mean[bed_subset[i, 4], ] = colMeans(m)
			mat_sd[bed_subset[i, 4], ] = colSds(m)
		}

		CpG_site_index = get_CpG_site_index(bed_subset[i, 2] - 10000, bed_subset[i, 2], site, site.index)
		if(length(CpG_site_index)) {
			m = get_meth(bsseq, CpG_site_index, SAMPLE$id)
			mat_mean_upstream[bed_subset[i, 4], ] = colMeans(m)
			mat_sd_upstream[bed_subset[i, 4], ] = colSds(m)
		}

		CpG_site_index = get_CpG_site_index(bed_subset[i, 3], bed_subset[i, 3] + 10000, site, site.index)
		if(length(CpG_site_index)) {
			m = get_meth(bsseq, CpG_site_index, SAMPLE$id)
			mat_mean_downstream[bed_subset[i, 4], ] = colMeans(m)
			mat_sd_downstream[bed_subset[i, 4], ] = colSds(m)
		}
	}
	rm(bsseq)
	gc(verbose = FALSE)
}

gene_methylation = list(
	mat_mean = mat_mean,
	mat_sd = mat_sd,
	mat_mean_upstream = mat_mean_upstream,
	mat_sd_upstream = mat_sd_upstream,
	mat_mean_downstream = mat_mean_downstream,
	mat_sd_downstream = mat_sd_downstream
)

save(gene_methylation, file = "RData/methylation_for_all_genes.RData")

