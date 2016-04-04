
# # old version of bsseq which is used in the pipeline
# library(parallel)
# library(MASS)
# library(grid)
# library(rjson)  # for reading gencode data
# library(DESeq2)  # differential expression analysis
# library(gplots)  # venn plot

# library(RColorBrewer)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(png)
# library(Gviz)
# library(matrixStats)

# library(GTF) # /home/guz/project/development/GTF_0.0.1.tar.gz
# library(GetoptLong) 
# library(circlize)
# library(ComplexHeatmap)  # devtools::install_github("jokergoo/ComplexHeamtap")
# library(EnrichedHeatmap) # devtools::install_github("jokergoo/EnrichedHeamtap")
# library(HilbertCurve) # devtools::install_github("jokergoo/HilbertCurve")
# library(GlobalOptions)
# library(gtrellis)

# library(ConsensusClusterPlus)
# for(f in dir("~/project/development/cotools/R")) {
# 	source(paste0("~/project/development/cotools/R/", f))
# }


load_package = function(path) {
	pkg = basename(path)
	if(paste0("package:", pkg) %in% search()) {
		detach(paste0("package:", pkg), unload = TRUE, character.only = TRUE)
	}
	tmp_dir = tempfile(pattern = "Rtmp_dir_", tmpdir = "/icgc/dkfzlsdf/analysis/B080/guz/temp", fileext = "")
	dir.create(tmp_dir)
	install.packages(path, repos = NULL, lib = tmp_dir)
	library(pkg, character.only = TRUE, lib.loc = tmp_dir)
}

load_package("~/project/development/cotools")
