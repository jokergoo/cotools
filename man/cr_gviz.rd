\name{cr_gviz}
\alias{cr_gviz}
\title{
Customized Gviz plot for a gene model
}
\description{
Customized Gviz plot for a gene model
}
\usage{
cr_gviz(cr, gi, expr, txdb, gene_start = NULL, gene_end = NULL, tx_list = NULL,
    species = "hg19", gf_list = NULL, hm_list = NULL, symbol = NULL)
}
\arguments{

  \item{cr}{correlated regions}
  \item{gi}{gene id}
  \item{expr}{expression matrix}
  \item{txdb}{txDb object}
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
