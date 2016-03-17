\name{annotate_to_gene_models}
\alias{annotate_to_gene_models}
\title{
Annotate to gene models
}
\description{
Annotate to gene models
}
\usage{
annotate_to_gene_models(gr, txdb, gene_model =c("tx", "gene"),
    species = "hg19", promoters_upstream = 2000, promoters_downstream = 200,
    annotation_type = c("percent", "number"),
    annotation_prefix = "overlap_to_")
}
\arguments{

  \item{gr}{a \code{\link[GenomicRanges]{GRanges}} object}
  \item{txdb}{a \code{TranscritpDb} object. At least the object should contain 'gene_id', 'tx_name'. ('gene_id' and 'tx_name' columns are used to identify genes and transcripts)}
  \item{gene_model}{type of gene model. By transcripts or by genes}
  \item{species}{We need this information to find out proper intergenic regions}
  \item{promoters_upstream}{length of upstream promoter from TSS, pass to \code{\link[GenomicRanges]{promoters}}}
  \item{promoters_downstream}{length of downstream promoter from TSS, pass o \code{\link[GenomicRanges]{promoters}}}
  \item{annotation_type}{Pass to \code{\link{annotate_to_genomic_features}}}
  \item{annotation_prefix}{Pass to \code{\link{annotate_to_genomic_features}}}

}
\value{
Following columns are attached to \code{gr}:

\describe{
  \item{nearest_tss}{the nearest tss (depending on \code{gene_model})}
  \item{dist_to_tss}{distance to the closest tss (depending on \code{gene_model})}
  \item{nearest_gm}{the closest gene model (depending on \code{gene_model})}
  \item{dist_to_gm}{distance to teh closest gene model (depending on \code{gene_model})}
  \item{prefix_to_exon}{percent of the region which is covered by exon}
  \item{prefix_to_intron}{percent of the region which is covered by intron}
  \item{prefix_to_promoter}{percent of the region which is covered by promoter}
  \item{prefix_to_intergenic}{percent of the region which is covered by intergenic regions}
  \item{prefix_to_fiveUTR}{percent of the region which is covered by 5'UTR}
  \item{prefix_to_threeUTR}{percent of the region which is covered by 3'UTR}
}
}
\examples{
# There is no example
NULL
}
