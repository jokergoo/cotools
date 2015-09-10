
# universal
# 'normal' is a hard coded type
sample_text = 
"id   type     WGBS  RNASEQ   CHIPSEQ    col
AK015    IDH    TRUE    TRUE    FALSE    #FF0000
AK041    IDH    TRUE    TRUE    FALSE    #FF0000
AK066    IDH    TRUE    TRUE    TRUE    #FF0000
AK068    IDH    TRUE    TRUE    FALSE    #FF0000
AK076    IDH    TRUE    TRUE    TRUE    #FF0000
AK085    IDH    TRUE    TRUE    FALSE    #FF0000
AK102    IDH    TRUE    TRUE    FALSE    #FF0000
AK103    IDH    TRUE    TRUE    FALSE    #FF0000
AK124    IDH    TRUE    TRUE    TRUE    #FF0000
AK199    IDH    TRUE    TRUE    TRUE    #FF0000
AK213    IDH    TRUE    TRUE    TRUE    #FF0000
AK231    IDH    TRUE    TRUE    TRUE    #FF0000
AK030    MES    TRUE    TRUE    FALSE    #00FF00
AK055    MES    TRUE    TRUE    FALSE    #00FF00
AK071    MES    TRUE    TRUE    TRUE    #00FF00
AK072    MES    TRUE    TRUE    FALSE    #00FF00
AK091    MES    TRUE    TRUE    TRUE    #00FF00
AK139    MES    TRUE    TRUE    TRUE    #00FF00
AK153    MES    TRUE    TRUE    TRUE    #00FF00
AK168    MES    TRUE    TRUE    TRUE    #00FF00
AK185    MES    TRUE    TRUE    FALSE    #00FF00
AK195    MES    TRUE    TRUE    FALSE    #00FF00
AK227    MES    TRUE    TRUE    FALSE    #00FF00
AK236    MES    TRUE    TRUE    FALSE    #00FF00
AK003    PDGFRA    TRUE    TRUE    FALSE    #0000FF
AK017    PDGFRA    TRUE    TRUE    FALSE    #0000FF
AK049    PDGFRA    TRUE    TRUE    FALSE    #0000FF
AK051    PDGFRA    TRUE    TRUE    FALSE    #0000FF
AK142    PDGFRA    TRUE    TRUE    FALSE    #0000FF
AK149    PDGFRA    TRUE    TRUE    TRUE    #0000FF
AK156    PDGFRA    TRUE    TRUE    TRUE    #0000FF
AK165    PDGFRA    TRUE    TRUE    FALSE    #0000FF
AK173    PDGFRA    TRUE    TRUE    TRUE    #0000FF
AK183    PDGFRA    TRUE    TRUE    TRUE    #0000FF
AK203    PDGFRA    TRUE    TRUE    FALSE    #0000FF
AK217    PDGFRA    TRUE    TRUE    FALSE    #0000FF
AK053    RTK_II    TRUE    TRUE    FALSE    #00FFFF
AK074    RTK_II    TRUE    TRUE    FALSE    #00FFFF
AK088    RTK_II    TRUE    TRUE    FALSE    #00FFFF
AK089    RTK_II    TRUE    TRUE    TRUE    #00FFFF
AK098    RTK_II    TRUE    TRUE    FALSE    #00FFFF
AK100    RTK_II    TRUE    TRUE    TRUE    #00FFFF
AK132    RTK_II    TRUE    TRUE    FALSE    #00FFFF
AK158    RTK_II    TRUE    TRUE    TRUE    #00FFFF
AK167    RTK_II    TRUE    TRUE    FALSE    #00FFFF
AK178    RTK_II    TRUE    TRUE    TRUE    #00FFFF
AK188    RTK_II    TRUE    TRUE    FALSE    #00FFFF
AK216    RTK_II    TRUE    TRUE    TRUE    #00FFFF
normal_occipital     normal  TRUE    TRUE  FALSE   grey
normal_parietal     normal  TRUE    TRUE  FALSE   grey
normal_temporal     normal  TRUE    TRUE  FALSE   grey
normal_frontal     normal  TRUE    TRUE   FALSE  grey
"

SAMPLE = read.table(textConnection(sample_text), header = TRUE, stringsAsFactors = FALSE, comment.char = "")
rownames(SAMPLE) = SAMPLE$id

SAMPLE_COLOR = unique(SAMPLE$col)
names(SAMPLE_COLOR) = unique(SAMPLE$type)



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


###########################################################################

# wgbs_qcplot("AK100")
# plot_coverage_and_methylation_on_genome("AK100", chromosome = c("chr21", "chr22"))

# global_methylation_distribution(SAMPLE$id)
# global_methylation_distribution(SAMPLE$id, by_chr = TRUE)


##############
if(!exists("gene_annotation", envir = .GlobalEnv, inherits = FALSE)) {
    cat("Loading gencode...\n")
    load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/data/gencode_v19_transcript_merged.RData")
}
if(!exists("expression", envir = .GlobalEnv, inherits = FALSE)) {
    cat("load expression...\n")
    load("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/data/hipo_016_rnaseq_gencode19expression.RData")
}
if(!exists("txdb", envir = .GlobalEnv, inherits = FALSE)) {
    cat("load txdb...\n")
    txdb = loadDb("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/data/gencode.v19.sqlite")
}

gt = sapply(gene_annotation$gtf, function(x) x$type)
gt = gt[gt == "protein_coding"]

expr = expression$rpkm[names(gt), ]
count = expression$count[names(gt), ]
l = apply(count, 1, function(x) sum(x > 0) > length(x)/2)

SAMPLE = SAMPLE[SAMPLE$type != "normal", , drop = FALSE]

expr = expr[l, SAMPLE$id, drop = FALSE]


chromosome = paste0("chr", c(1:22, "X"))


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
        tfbs                   = qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/analysis/WGBS_final/bed/tfbs_cluster.bed")        # too large for memory
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
    gsub(qq("_@{hm}"), "", hm_sample)
}

get_hm = function(hm, sid) {
    df = read.table(qq("/icgc/dkfzlsdf/analysis/hipo/hipo_016/chipseq_gbm/@{sid}_@{hm}/sicer/@{sid}_@{hm}.mkdup-W200-G200-islands-summary-FDR0.01"))
    gr = GRanges(seqnames = df[[1]], ranges = IRanges(df[[2]], df[[3]]), density = df[[7]])
}

cp = initialize_project(
    sample_id = SAMPLE$id,
    type = SAMPLE$type,
    col = SAMPLE$col,
    chromosome = chromosome,
    species = "hg19",
    methylation_hooks = methylation_hooks,
    expr = expr,
    txdb = txdb,
    genomic_features = GENOMIC_FEATURE_LIST,
    chipseq = list(H3K4me3 = get_hm_sample("H3K4me3"), 
                   H3K4me1 = get_hm_sample("H3K4me1"),
                   H3K27Ac = get_hm_sample("H3K27Ac"),
                   H3K27me3 = get_hm_sample("H3K27me3"),
                   H3K36me3 = get_hm_sample("H3K36me3"),
                   H3K9me3 = get_hm_sample("H3K9me3")),
    get_hm = get_hm
)
