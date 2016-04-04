\name{enriched_heatmap_list_on_gene}
\alias{enriched_heatmap_list_on_gene}
\title{
Enriched heatmap
}
\description{
Enriched heatmap
}
\usage{
enriched_heatmap_list_on_gene(cr, cgi, txdb, expr, hm_list = NULL, hm_name = NULL, on = "tss", by = "gene",
    hm_cor_p_cutoff = 0.05, show_expr = TRUE, ...)
}
\arguments{

  \item{cr}{cr}
  \item{cgi}{CpG Island}
  \item{txdb}{txdb}
  \item{expr}{expr}
  \item{hm_list}{a list of histome marks}
  \item{hm_name}{names for hm}
  \item{on}{tss or body}
  \item{by}{gene or tx}
  \item{hm_cor_p_cutoff}{cutoff for the correlation between hm and expression}
  \item{show_expr}{whether show heatmap of expression}
  \item{...}{pass to}

}
\examples{
# There is no example
NULL

}
