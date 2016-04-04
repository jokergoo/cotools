\name{scatterplot_with_boxplot}
\alias{scatterplot_with_boxplot}
\title{
scatterplot with boxplots on both sides
}
\description{
scatterplot with boxplots on both sides
}
\usage{
scatterplot_with_boxplot(x, y, annotation = rep("unknown", length(x)),
    annotation_color = structure(seq_along(levels(annotation)), names = levels(annotation)),
    main = NULL, xlab = NULL, ylab = NULL, xlim = range(x), ylim = range(y), text_list = NULL)
}
\arguments{

  \item{x}{x}
  \item{y}{y}
  \item{annotation}{annotations}
  \item{annotation_color}{colors for annotation}
  \item{main}{title for the plot}
  \item{xlab}{xlab}
  \item{ylab}{ylab}
  \item{xlim}{xlim}
  \item{ylim}{ylim}
  \item{text_list}{additional text}

}
\examples{
# There is no example
NULL

}
