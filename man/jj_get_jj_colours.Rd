% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jj_get_jj_colours.R
\name{jj_get_jj_colours}
\alias{jj_get_jj_colours}
\title{get a set of predefined colours for a vector}
\usage{
jj_get_jj_colours(vec, symbols = FALSE)
}
\arguments{
\item{vec}{vector containing the names, for which colours should be assinged}

\item{symbols}{if TRUE, assign numbers instead of colors (useful if groups should be encoded by symbols)}
}
\description{
Returns a vector of colors with names according to input vector levels.
}
\examples{
jj_get_jj_colours(c("a","d","c","a"))
jj_get_jj_colours(c("a","d","c","a"), symbols=TRUE)
}
\keyword{colours}
