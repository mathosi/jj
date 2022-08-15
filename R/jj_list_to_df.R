#' Convert list to data.frame in long/wide format
#'
#' load excel workbook as list of data.frames or save list of data.frames as excel
#'
#' @name jj_list_to_df
#' @param flist Full path to the file to read or write
#' @param long If TRUE (default), return two column data.frame with names in first and values in second column. Else, return data.frame with one column per list element.
#' @param fill_with If long=FALSE, value to fill cells in case of unequal lenghts of list elements. Default=NA
#' @param maintain_classes keep the class of each entry in the list, otherwise coerces to character columns
#' @keywords list,data.frame
#' @export
#' @examples
#' 

jj_list_to_df = function(flist, long=TRUE, fill_with=NA, maintain_classes=TRUE){
  #fill_with is used to make columns equally long when long=FALSE
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
