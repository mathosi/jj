#' describe all the columns in a data.frame
#'
#' Gives an overview of a data.frame (dimensions, class of each column, mode/min/mean/max, distinct values, na_count, table, head)
#'
#' @param df data.frame
#' @param convert if convert = T, convert character columns that only contain digits to numeric columns before summarising
#' @param return_df if return_df = T, return the summary tibble instead of printing it out
#' @export
#' @examples
#' jj_describe_df(pbmc_small[[]])

jj_describe_df = function(df, convert=FALSE, return_df=FALSE){
  
  cat(sprintf('nrows=%i\nncols=%i\nnentries=%i', nrow(df), ncol(df), nrow(df)*ncol(df)))
  #convert factors to character
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  
  res = purrr::imap_dfr(df, .report, convert=convert)
  if(return_df){
    return(res)
  }else{
    print(res, n=100)
  } 
}


.report <- function(x, name, convert=F) { 
  ro = function(v){ifelse(is.numeric(v),round(v, 2),v)}
  # getmode <- function(v) {
  #   uniqv <- unique(na.omit(v))
  #   uniqv[which.max(tabulate(match(v, uniqv)))]
  # }
  gettable = function(v) {
    res = table(v) %>% sort(., decreasing = T) %>% head
    paste(mapply(function(x, y) paste(x,y, sep=': '), names(res), res, USE.NAMES = F), 
          collapse=', ')
  }
  mode_x = mode(x)
  
  if(convert){
    #Tests, without issuing warnings, whether all elements of a character vector are 
    #legal numeric values, or optionally converts the vector to a numeric vector
    x = Hmisc::all.is.numeric(x, what='vector')
  }
  if(Hmisc::all.is.numeric(x)){
    tibble(
      name  = name, 
      type=mode_x,
      `mode/min/mean/max` = paste(c(ro(get_mode(x)),ro(min(x)),
                                    ro(suppressWarnings(mean(x, na.rm = T))),
                                    ro(max(x))),collapse = '/'), 
      n_distinct = dplyr::n_distinct(x),
      na_count = sum(is.na(x)),
      table = gettable(x),
      head=paste(head(x), collapse = ',')
    )
  }else{
    tibble(
      name  = name, 
      type=mode_x,
      `mode/min/mean/max`  = paste(c( get_mode(x),'','',''),collapse = '/'), 
      n_distinct = dplyr::n_distinct(x),
      na_count = sum(is.na(x)),
      table = gettable(x),
      head=paste(head(x), collapse = ',')
    )
  }
}
