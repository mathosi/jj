#' Convert list to data.frame in long/wide format
#'
#' @name jj_list_to_df
#' @param flist Full path to the file to read or write
#' @param long If TRUE (default), return two column data.frame with names in first and values in second column. Else, return data.frame with one column per list element.
#' @param fill_with If long=FALSE, value to fill cells in case of unequal lenghts of list elements. Default=NA
#' @param maintain_classes keep the class of each entry in the list, otherwise coerces to character columns
#' @export
#' @examples
#' example_list = list(signature1=c('CD19A','CD79B','SDC1'), signature2=c('CD8A','CD8B1'))
#' jj_list_to_df(flist=example_list)
#' jj_list_to_df(flist=example_list, long = F)

jj_list_to_df = function(flist, long=TRUE, fill_with=NA, maintain_classes=TRUE){
  sets_length = sapply(flist, length)
  if(long){
    df = data.frame(name = rep(names(flist), times = sets_length), feature=unlist(flist))
  }else{
    eq_sets_length_list <- lapply(flist, function(x) c(x, rep(fill_with, max(sets_length) - length(x))))
    df <- as.data.frame(do.call(cbind, eq_sets_length_list))
    if(maintain_classes){
      classes = sapply(flist, class)
      for(i in 1:ncol(df)){
        class(df[, i]) = classes[i]
      }
    }
  }
  rownames(df) = NULL
  return(df)
}
