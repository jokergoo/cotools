 
suppressPackageStartupMessages(library(GetoptLong))

head = "~/project/development/cotools/pipeline/head/head.R"
whichxx = "neg"
GetoptLong(c("head=s", "head R script",
	         "histone_mark=s", "histone mark",
	         "whichxx=s", "pos|neg"))

source("~/project/development/cotools/script/load_all.R")
source(head)

cr_filtered = readRDS(qq("@{RDS_FOLDER}/cr_filtered_fdr_0.01.rds"))

make_enriched_heatmap = function(histone_mark, by, on, which, type = 1:3) {
	
	hm_list = get_hm_list(histone_mark)

	if(which == "neg") {
	    cr = cr_filtered[cr_filtered$corr < 0]
	} else {
	    cr = cr_filtered[cr_filtered$corr > 0]
	}

	n_factor = length(unique(attr(cr, "factor")[attr(cr, "sample_id") %in% names(hm_list)]))
	if(n_factor == 0) n_factor = 1

	ht_global_opt(heatmap_column_title_gp = gpar(fontsize = 10))
	if(any(type %in% c(1, 2, 3))) {
		qqcat("heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}\n")
		pdf(qq("@{output_dir}/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf"), width = 8 + n_factor*2, height = 10)
		enriched_heatmap_list_on_gene(cr, GENOMIC_FEATURE_LIST$cgi, txdb, log2(expr+1), hm_list, 
		    hm_name = histone_mark, on = on, by = by)
		dev.off()

		system(qq("convert @{output_dir}/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.pdf @{output_dir}/heatmap_simple_@{by}_@{on}_@{which}_@{histone_mark}.png"))
	}

	if(any(type %in% c(3, 2))) {
		qqcat("heatmap_cgi_@{by}_@{which}_@{histone_mark}\n")
		pdf(qq("@{output_dir}/heatmap_cgi_@{by}_@{which}_@{histone_mark}.pdf"), width = 6 + n_factor*2, height = 10)
		enriched_heatmap_list_on_tss_cgi(cr, GENOMIC_FEATURE_LIST$cgi, txdb, log2(expr+1), hm_list, 
		    hm_name = histone_mark, by = by)
		dev.off()

		system(qq("convert @{output_dir}/heatmap_cgi_@{by}_@{which}_@{histone_mark}.pdf @{output_dir}/heatmap_cgi_@{by}_@{which}_@{histone_mark}.png"))
	}

	if(any(type %in% 3)) {
		qqcat("heatmap_encode_tfbs_@{which}_@{histone_mark}")
		pdf(qq("@{output_dir}/heatmap_encode_tfbs_@{which}_@{histone_mark}.pdf"), width = 6 + n_factor*2, height = 10)
		enriched_heatmap_list_on_genomic_features(cr, GENOMIC_FEATURE_LIST$tfbs, hm_list, 
		    hm_name = histone_mark)
		dev.off()

		system(qq("convert @{output_dir}/heatmap_encode_tfbs_@{which}_@{histone_mark}.pdf heatmap_encode_tfbs_@{which}_@{histone_mark}.png"))

		qqcat("heatmap_encode_strong_enhancer_@{which}_@{histone_mark}")
		pdf(qq("@{output_dir}/heatmap_encode_strong_enhancer_@{which}_@{histone_mark}.pdf"), width = 6 + n_factor*2, height = 10)
		enriched_heatmap_list_on_genomic_features(cr, GENOMIC_FEATURE_LIST$enhancer, hm_list, 
		    hm_name = histone_mark)
		dev.off()

		system(qq("convert @{output_dir}/heatmap_encode_strong_enhancer_@{which}_@{histone_mark}.pdf heatmap_encode_strong_enhancer_@{which}_@{histone_mark}.png"))
	}
}

make_enriched_heatmap(histone_mark, by = "gene", on = "tss", which = whichxx, type = 1:3)
make_enriched_heatmap(histone_mark, by = "tx", on = "tss", which = whichxx, type = 1:2)
make_enriched_heatmap(histone_mark, by = "gene", on = "body", which = whichxx, type = 1)


# for(hm in c("H3K4me3", "H3K27Ac", "H3K27me3", "H3K36me3", "H3K4me1", "H3K9me3")) {
# 	cmd = qq("Rscript-3.1.2 ~/project/development/cotools/pipeline/cr_enriched.R --histone_mark @{hm} --whichxx neg")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=18G -N gbm_enriched_heatmap_@{hm}_neg' '@{cmd}'")
# 	system(cmd)
#   cmd = qq("Rscript-3.1.2 ~/project/development/cotools/pipeline/cr_enriched.R --histone_mark @{hm} --whichxx pos")
# 	cmd = qq("perl /home/guz/project/development/ngspipeline2/qsub_single_line.pl '-l walltime=40:00:00,mem=18G -N gbm_enriched_heatmap_@{hm}_pos' '@{cmd}'")
# 	system(cmd)
# }
