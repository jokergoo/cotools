\name{wgbs_qcplot}
\alias{wgbs_qcplot}
\title{
Basic qc plot for bisulfite sequencing data
}
\description{
Basic qc plot for bisulfite sequencing data
}
\usage{
wgbs_qcplot(sample_id, chromosome = paste0("chr", 1:22))
}
\arguments{

  \item{sample_id}{a vector of sample ids}
  \item{chromosome}{a vector of chromosomes}

}
\details{
For each sample id, it will produce five plots:

1. mean/median CpG coverage per chromosome
2. histogram of CpG coverage
3. methylation per chromosome 
4. histogram of methylation
5. mean Methylation for each CpG coverage
}
\examples{
# There is no example
NULL

}
