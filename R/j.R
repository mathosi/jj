#' print small subset of R object as overview
#'
#' Works for all R objects in contrast to head().
#'
#' @param obj Object, for which an overview is desired
#' @param headtail return specified number of items from the start and the end of the object (not lists). Can result in duplicated items when headtail is bigger than the object length.
#' @param print_nested also print overview for lists within lists. Default=FALSE, Beta feature!
#' @keywords overview, head, glimpse
#' @export
#' @examples
#'

j <- function(obj, headtail=NULL, print_nested=FALSE, show_max=5){
  if(!is.null(headtail)){
    headtail = as.integer(headtail)
    stopifnot(!is.na(headtail))
  }
  
  check_vector = function(obj){
    message(paste0("Vector of length ", length(obj)))
    if(!is.null(headtail)){
      stopifnot(is.integer(headtail))
      return(c(head(obj, headtail), tail(obj, headtail)))
    }
    return(head(obj, show_max))
  }  
  
  check_matrix_data.frame = function(obj){
    nrows = nrow(obj)
    ncols = ncol(obj)
    message(sprintf("%s with %i rows and %i columns", class(obj)[1], nrows, ncols))
    if(!is.null(headtail)){
      return(rbind(head(obj, headtail), tail(obj, headtail)))
    }
    lookRows <- min(show_max, nrows)
    lookCols <- min(show_max, ncols)
    return(obj[1:lookRows, 1:lookCols])
  }
  
  check_DataFrame = function(obj){
    obj <- as.data.frame(obj)
    nrows = nrow(obj)
    ncols = ncol(obj)
    message(sprintf("DataFrame with %i rows and %i columns", nrows, ncols))
    if(!is.null(headtail)){
      return(rbind(head(obj, headtail), tail(obj, headtail)))
    }
    lookRows <- min(show_max, nrows)
    lookCols <- min(show_max, ncols)
    return(obj[1:lookRows,1:lookCols])
  }
  
  check_list = function(obj, print_nested){
    lookLists <- min(show_max, length(obj))
    message(sprintf("Length of list: %i, first %i elements:",length(obj), lookLists))
    lapply(obj[1:lookLists], function(x) {
      if(is.matrix(x) | is.data.frame(x) | is(x, 'sparseMatrix')){
        check_matrix_data.frame(x)
      }else if (class(x)=="DataFrame"){
        check_DataFrame(x)
      }else if(is.list(x)){
        message('List of length ', length(x))
        if(print_nested){
          check_list(x)
        }
      }else if(is.vector(x)){
        check_vector(x)
      }else{
        return(dplyr::glimpse(x))
      }
    })
  }
  
  if(is.matrix(obj) | is.data.frame(obj) | is(obj, 'sparseMatrix')){
    check_matrix_data.frame(obj, headtail=headtail)
  }else if(class(obj)=="DataFrame"){
    check_DataFrame(obj, headtail=headtail)
  }else if(isS4(obj)){
    if(!is.null(headtail)){
      warning('Dont know how to return head and tail for this')
    }
    return(obj)
  }else if(is.list(obj)){
    check_list(obj, print_nested=print_nested)
  }else if(is.vector(obj) | is.factor(obj)){
    #careful: lists are also vectors -> put after is.list
    check_vector(obj, headtail=headtail)
  }else{
    if(!is.null(headtail)){
      warning('Dont know how to return1 head and tail for this')
    }
    return(dplyr::glimpse(obj))
  }
}
