
# old version of bsseq which is used in the pipeline
library(parallel)
library(MASS)
library(grid)
library(rjson)  # for reading gencode data
library(DESeq2)  # differential expression analysis
library(gplots)  # venn plot

library(RColorBrewer)
library(GenomicRanges)
library(GenomicFeatures)
library(png)
library(Gviz)
library(matrixStats)

library(GTF) # /home/guz/project/development/GTF_0.0.1.tar.gz
library(GetoptLong) 
library(circlize)
library(ComplexHeatmap)  # devtools::install_github("jokergoo/ComplexHeamtap")
library(EnrichedHeatmap) # devtools::install_github("jokergoo/EnrichedHeamtap")
library(HilbertCurve) # devtools::install_github("jokergoo/HilbertCurve")
library(GlobalOptions)
library(gtrellis)

for(f in dir("~/project/development/cotools/R")) {
	source(paste0("~/project/development/cotools/R/", f))
}

