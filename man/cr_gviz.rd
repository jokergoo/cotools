\name{cr_gviz}
\alias{cr_gviz}
\title{
Gviz plot for a gene model with tracks
}
\description{
Gviz plot for a gene model with tracks
}
\usage{
cr_gviz(cr, gi, expr, txdb, gene_start = NULL, gene_end = NULL, tx_list = NULL,
    species = "hg19", gf_list = NULL, hm_list = NULL, symbol = NULL)
}
\arguments{

  \item{gi}{gene id}
  \item{expe}{expression matrix}
  \item{txdb}{transcriptDb object}
  \item{gene_start}{start of gene}
  \item{gene_end}{end of the gene}
  \item{tx_list}{a list of tx}
  \item{species}{species}
  \item{gf_list}{a list of gf}
  \item{hm_list}{a list of hm}
  \item{symbol}{gene symbol}

}
\examples{
# There is no example
NULL
}
