% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jj_list_to_df.R
\name{jj_list_to_df}
\alias{jj_list_to_df}
\title{Convert list to data.frame in long/wide format}
\usage{
jj_list_to_df(flist, long = TRUE, fill_with = NA, maintain_classes = TRUE)
}
\arguments{
\item{flist}{Full path to the file to read or write}

\item{long}{If TRUE (default), return two column data.frame with names in first and values in second column. Else, return data.frame with one column per list element.}

\item{fill_with}{If long=FALSE, value to fill cells in case of unequal lenghts of list elements. Default=NA}

\item{maintain_classes}{keep the class of each entry in the list, otherwise coerces to character columns}
}
\description{
Convert list to data.frame in long/wide format
}
\examples{
example_list = list(signature1=c('CD19A','CD79B','SDC1'), signature2=c('CD8A','CD8B1'))
jj_list_to_df(flist=example_list)
jj_list_to_df(flist=example_list, long = F)
}
