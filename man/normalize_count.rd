\name{normalize_count}
\alias{normalize_count}
\title{
normalize raw count
}
\description{
normalize raw count
}
\usage{
normalize_count(count, txdb = NULL, method = "rpkm",
    param = NULL, gene_length_type = c("exon", "gene"))
}
\arguments{

  \item{count}{count matrix, rownames are used to map genes in \code{txdb}}
  \item{txdb}{transcriptome, used if normalized by gene length}
  \item{method}{rpkm, voom, tpm, tc, md, deseq2, tmm}
  \item{param}{additional parameters for different normalization methods}
  \item{gene_length_type}{if normalize to gene length, type}

}
\details{
useful links:
\url{http://bib.oxfordjournals.org/content/14/6/671.full.pdf+html}
\url{https://www.biostars.org/p/56919/}
\url{http://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/}
}
\value{
normalized expression matrix
}
\examples{
# There is no example
NULL
}
