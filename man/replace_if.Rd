% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utilities.R
\name{replace_if}
\alias{replace_if}
\title{replace values that meet a condition}
\usage{
replace_if(vec, count_below = NULL, replacement)
}
\arguments{
\item{vec}{vector with values to replace}

\item{count_below}{replace values that occur less often than this threshold}

\item{replacement}{value that is used for replacement}
}
\description{
Currently only replacing values below a specified abundance is implemented
}
\examples{
example_vec = sample(letters[1:4], 1000,replace=T, prob=c(0.5,0.3,0.1,0.1))
table(example_vec)
replaced_vec = replace_if(example_vec, count_below=200, replacement = 'rare letter')
head(replaced_vec)
table(replaced_vec)
}
