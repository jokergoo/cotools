
# universal
# 'normal' is a hard coded type
sample_text = 
"id   type   age
AK015    IDH  28
AK041    IDH  47
AK066    IDH  34
AK068    IDH  59
AK076    IDH  47
AK085    IDH  38
AK102    IDH  49
AK103    IDH  40
AK124    IDH  32
AK199    IDH  33
AK213    IDH  48
AK231    IDH  56
AK030    MES  66
AK055    MES  55
AK071    MES  49
AK091    MES  58
AK139    MES  46
AK153    MES  63
AK185    MES  36
AK188    MES  55
AK195    MES  58
AK227    MES  55
AK236    MES  44
AK003    RTK_I  44
AK017    RTK_I  19
AK049    RTK_I  73
AK051    RTK_I  55
AK142    RTK_I  49
AK149    RTK_I  73
AK156    RTK_I  63
AK165    RTK_I  36
AK173    RTK_I  59
AK183    RTK_I  49
AK203    RTK_I  20
AK217    RTK_I  61
AK053    RTK_II  62
AK072    RTK_II  51
AK074    RTK_II  40
AK088    RTK_II  63
AK089    RTK_II  62
AK098    RTK_II  36
AK100    RTK_II  38
AK132    RTK_II  60
AK158    RTK_II  59
AK167    RTK_II  69
AK178    RTK_II  57
AK216    RTK_II  61
#normal_occipital     normal
#normal_parietal     normal
#normal_temporal     normal
#normal_frontal     normal
"

SAMPLE = read.table(textConnection(sample_text), header = TRUE, stringsAsFactors = FALSE)
rownames(SAMPLE) = SAMPLE$id

#SAMPLE_COLOR = c("IDH" = "red", "MES" = "green", "PDGFRA" = "blue", "RTK_II" = "orange", "normal" = "grey")
#SAMPLE_COLOR = c("IDH" = "red", "MES" = "green", "RTK_I" = "blue", "RTK_II" = "orange")
SAMPLE_COLOR = brewer.pal(4, "Set1")
names(SAMPLE_COLOR) = c("IDH", "MES", "RTK_I", "RTK_II")

AGE_COL_FUN = colorRamp2(c(20, 70), c("black", "white"))

############################################################
#
methylation_hooks$set = function(chr) {

    if(!is.null(methylation_hooks$obj)) {
        if(attr(methylation_hooks$obj, "chr") == chr) {
            qqcat("[@{chr}] @{chr} is already set.\n")
            return(invisible(NULL))
        }
    }

    qqcat("[@{chr}] loading /icgc/dkfzlsdf/analysis/hipo/hipo_016/wgbs_gbm_smoothed/smoothed_object_list_@{chr}.RData\n")
    
    load(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/wgbs_gbm_smoothed/smoothed_object_list_@{chr}.RData"))
    
    class(bs.fit) = "bs_fit"
    attr(bs.fit, "chr") = chr

    methylation_hooks$obj = bs.fit

    return(invisible(NULL))
}


methylation_hooks$meth = function(bs_fit = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

   if(is.null(row_index) && is.null(col_index)) {
        bs_fit$meth[, , drop = FALSE]
    } else if(is.null(row_index)) {
        bs_fit$meth[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        bs_fit$meth[row_index, , drop = FALSE]
    } else {
        bs_fit$meth[row_index, col_index, drop = FALSE]
    }

}

methylation_hooks$raw = function(bs_fit = methylation_hooks$obj, row_index = NULL, col_index = NULL) {

   if(is.null(row_index) && is.null(col_index)) {
        bs_fit$raw[, , drop = FALSE]
    } else if(is.null(row_index)) {
        bs_fit$raw[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        bs_fit$raw[row_index, , drop = FALSE]
    } else {
        bs_fit$raw[row_index, col_index, drop = FALSE]
    }
}


methylation_hooks$site = function(bs_fit = methylation_hooks$obj, index = NULL) {
    if(is.null(index))
        bs_fit$site
    else bs_fit$site[index]
}

methylation_hooks$GRanges = function(bs_fit = methylation_hooks$obj) {
    chr = attr(bs_fit, "chr")
    site = bs_fit$site
    GRanges(seqnames = rep(chr, length(site)), ranges = IRanges(start = site, end = site))
}

methylation_hooks$coverage = function(bs_fit = methylation_hooks$obj, 
    row_index = NULL, col_index = NULL) {
    
    if(is.null(row_index) && is.null(col_index)) {
        bs_fit$coverage[, , drop = FALSE]
    } else if(is.null(row_index)) {
        bs_fit$coverage[, col_index, drop = FALSE]
    } else if(is.null(col_index)) {
        bs_fit$coverage[row_index, , drop = FALSE]
    } else {
        bs_fit$coverage[row_index, col_index, drop = FALSE]
    }
}


##############
gencode_gtf_file = "/icgc/dkfzlsdf/analysis/B080/guz/gencode/gencode.v19.annotation.gtf"

if(!exists("expression", envir = .GlobalEnv, inherits = FALSE)) {
    cat("load expression...\n")
    load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/data/hipo_016_rnaseq_gencode19expression.RData")
}
if(!exists("txdb", envir = .GlobalEnv, inherits = FALSE)) {
    cat("load txdb...\n")
    txdb = loadDb("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/data/gencode.v19.sqlite")
}

gt = extract_field_from_gencode(gencode_gtf_file, level = "gene", primary_key = "gene_id", field = "gene_type")
gt = gt[gt == "protein_coding"]

expr = expression$rpkm[names(gt), ]
count = expression$count[names(gt), ]
l = apply(count, 1, function(x) sum(x > 0) > length(x)/2)

SAMPLE = SAMPLE[SAMPLE$type != "normal", , drop = FALSE]

expr = expr[l, SAMPLE$id, drop = FALSE]

genes = genes(txdb)
chr = as.vector(seqnames(genes))
names(chr) = names(genes)

chromosome = paste0("chr", c(1:22))
expr = expr[chr[rownames(expr)] %in% chromosome, ]

genome = "hg19"

output_dir = "/home/guz/project/analysis/hipo16_new/figure_prepare/"

dir.create(output_dir, showWarnings = FALSE)
dir.create(paste0(output_dir, "/gviz"), showWarnings = FALSE)

if(!exists("GENOMIC_FEATURE_LIST", envir = .GlobalEnv)) {
    # annotate to other regions
    GENOMIC_FEATURE_LIST = c(
        gene                   = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/gene.bed"),
        exon                   = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/exon.bed"),
        intron                 = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/intron.bed"),
        tss_2k                 = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/tss_2k.bed"),
        intergenic             = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/intergenic.bed"),
        cgi                    = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/cpgIslandExt.bed"),
        cgi_shore              = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/cgi_shore_2k.bed"),
        dnase                  = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/dnase.bed"),
        enhancer               = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/enhancer.bed"),
        repeats_LINE           = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/repeats_LINE.bed"),
        repeats_SINE           = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/repeats_SINE.bed"),
        tfbs                   = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/encode_uniform_tfbs_merged_1kb.bed")        # too large for memory
    )
    cat("loading genomic features...\n")
    GENOMIC_FEATURE_LIST = lapply(GENOMIC_FEATURE_LIST, function(x) {
        qqcat("reading @{x} as GRanges object...\n")
        df = read.table(x, sep = "\t", stringsAsFactors = FALSE)
        df = df[df[[1]] %in% chromosome, , drop = FALSE]
        makeGRangesFromDataFrame(df[1:3], 
                    seqnames.field = "V1",
                    start.field = "V2",
                    end.field = "V3")
    })
}

get_hm_sample = function(hm) {
    hm_sample = dir("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/", pattern = qq("AK.*@{hm}"))
    sample_id = gsub(qq("_@{hm}"), "", hm_sample)
    intersect(SAMPLE$id, sample_id)
}

get_hm = function(sid, histone_mark) {
    df = read.table(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/@{sid}_@{histone_mark}/sicer/@{sid}_@{histone_mark}.mkdup-W200-G200-islands-summary-FDR0.01"), stringsAsFactors = FALSE)
    df = df[df[[1]] %in% chromosome]
    gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), density = df[[7]])
}

# cp = initialize_project(
#     sample_id = SAMPLE$id,
#     type = SAMPLE$type,
#     col = SAMPLE_COLOR[SAMPLE$type],
#     chromosome = chromosome,
#     species = "hg19",
#     methylation_hooks = methylation_hooks,
#     expr = expr,
#     gencode_gtf_file = gencode_gtf_file,
#     txdb = txdb,
#     genomic_features = GENOMIC_FEATURE_LIST,
#     chipseq = list(H3K4me3 = get_hm_sample("H3K4me3"), 
#                    H3K4me1 = get_hm_sample("H3K4me1"),
#                    H3K27Ac = get_hm_sample("H3K27Ac"),
#                    H3K27me3 = get_hm_sample("H3K27me3"),
#                    H3K36me3 = get_hm_sample("H3K36me3"),
#                    H3K9me3 = get_hm_sample("H3K9me3")),
#     get_hm = get_hm
# )
