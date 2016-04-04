\name{extract_field_from_gencode}
\alias{extract_field_from_gencode}
\title{
Extract field from gencode GTF file
}
\description{
Extract field from gencode GTF file
}
\usage{
extract_field_from_gencode(file, level = "gene", primary_key = "gene_id", field = "gene_name")
}
\arguments{

  \item{file}{the input gtf file}
  \item{level}{level of the annotation (e.g. gene, transcript, exon, ...)}
  \item{primary_key}{primary field}
  \item{field}{field to be retrieved}

}
\details{
Although gtf file can be imported by \code{\link[GenomicFeatures]{makeTranscriptDbFromGFF}}, some information
in the original gtf file will not be imported. This function aims to extract additionally information
from gtf file.

The function calls external perl script, so you need to have perl installed.
}
\value{
A vector in which 'primary_key' corresponds to the name and 'field' corresponds to the value
}
\examples{
# There is no example
NULL

}
