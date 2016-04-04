library(GetoptLong)

head = "~/project/development/cotools/pipeline/head/head.R"
chrxx = paste0("chr", 1:22)
GetoptLong("head=s", "head R script",
	       "chrxx=s{1,}", "chromosomes")

source("~/project/development/cotools/script/load_all.R")
source(head)

cr_filtered = readRDS(qq("@{RDS_FOLDER}/cr_filtered_fdr_0.01.rds"))

attr(cr_filtered, "col") = SAMPLE_COLOR

######################## gviz ########################
hm_list = get_hm_list("H3K4me3")

gn = extract_field_from_gencode(gencode_gtf_file, level = "gene", primary_key = "gene_id", field = "gene_name")
tx_list = transcriptsBy(txdb, by = "gene")
gi = "ENSG00000060982.10"
for(chr in intersect(chrxx, unique(as.character(seqnames(cr_filtered))))) {
	cat(chr, "...\n")
    cr_subset = cr_filtered[seqnames(cr_filtered) == chr]
    for(gi in unique(cr_subset$gene_id)) {
        pdf(qq("@{output_dir}/gviz/gviz_@{chr}_@{gi}_@{gn[gi]}.pdf"), width = 16, height = 12)
        cr_gviz(cr_filtered, gi, expr, txdb, tx_list = tx_list[[gi]]$tx_name, gf_list = GENOMIC_FEATURE_LIST[c("cgi", "tfbs")], 
        	hm_list = hm_list, symbol = gn[gi])
        dev.off()
    }
}
##############################################
