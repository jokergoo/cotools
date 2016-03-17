\name{chromatin_states_transition_chord_diagram}
\alias{chromatin_states_transition_chord_diagram}
\title{
chord diagram for chromatin states transistion
}
\description{
chord diagram for chromatin states transistion
}
\usage{
chromatin_states_transition_chord_diagram(mat, max_mat = mat, cate1, cate2, ...)
}
\arguments{

  \item{mat}{the transition matrix}
  \item{max_mat}{if there are several matrix, set it to the matrix with maximum sum}
  \item{cate1}{name of row states}
  \item{cate2}{name fo column states}
  \item{...}{pass to \code{\link[circlize]{chordDiagram}}}

}
\examples{
# There is no example
NULL
}
