#' print small subset of R object as overview
#'
#' Works for all R objects in contrast to head().
#'
#' @param obj Object, for which an overview is desired
#' @keywords overview, head, glimpse
#' @export
#' @examples
#'

jj <- function(obj, headtail=NULL){
  if(!is.null(headtail)){
    headtail <- as.integer(headtail)
    stopifnot(is.integer(headtail))
  }
  if(is.matrix(obj) | is.data.frame(obj)){
    message(paste0("Dimensions of ", class(obj), ": ", paste(dim(obj),collapse = ", ")))
    if(!is.null(headtail)){
      return(rbind(head(obj, headtail), tail(obj, headtail)))
    }
    lookRows <- ifelse(dim(obj)[1]>4, 5, dim(obj)[1])
    lookCols <- ifelse(dim(obj)[2]>4, 5, dim(obj)[2])
    return(obj[1:lookRows,1:lookCols])
  }else if(is.list(obj)){
    message(paste0("Length of list: ",length(obj)))
    lookLists <- ifelse(length(obj)>4, 5, length(obj))
    lapply(obj[1:lookLists], class)
    lapply(obj[1:lookLists], function(x) {
      if(is.matrix(x) | is.data.frame(x)){
        lookRows <- ifelse(dim(x)[1]>4, 5, dim(x)[1])
        lookCols <- ifelse(dim(x)[2]>4, 5, dim(x)[2])
        if(!is.null(headtail)){
          return(rbind(head(x, headtail), tail(x, headtail)))
        }
        message(paste0("   Dimensions of ", class(x), ": ", paste(dim(x),collapse = ", ")))
        return(x[1:lookRows,1:lookCols])
        
      }else if (class(x)=="DataFrame"){
        x <- as.data.frame(x)
        lookRows <- ifelse(dim(x)[1]>4, 5, dim(x)[1])
        lookCols <- ifelse(dim(x)[2]>4, 5, dim(x)[2])
        if(!is.null(headtail)){
          return(rbind(head(x, headtail), tail(x, headtail)))
        }
        message(paste0("   Dimensions of DataFrame: ", paste(dim(x),collapse = ", ")))
        return(x[1:lookRows,1:lookCols])
      }else if(is.vector(x)){
        if(!is.null(headtail)){
          return(c(head(x, headtail), tail(x, headtail)))
        }
        message(paste0("   Length of vector: ", length(x)))
        return(head(x))
      }else if(class(x) %in% c('dgCMatrix')){
        if(!is.null(headtail)){
          return(rbind(head(x, headtail), tail(x, headtail)))
        }
        message(paste0("Dimensions of ", class(x), ": ", paste(dim(x),collapse = ", ")))
        return(print(x[1:4, 1:4]))
      }else{
        library(dplyr)
        if(!is.null(headtail)){
          warning('Dont know how to return head and tial for this')
        }
        return(dplyr::glimpse(x))
      }
    })
  }else if(is.vector(obj)){
    message(paste0("Length of vector: ", length(obj)))
    if(!is.null(headtail)){
      stopifnot(is.integer(headtail))
      return(c(head(obj, headtail), tail(obj, headtail)))
    }
    return(head(obj))
  }else if(class(obj)=="DataFrame"){
    obj <- as.data.frame(obj)
    message(paste0("Dimensions of DataFrame: ", paste(dim(obj),collapse = ", ")))
    if(!is.null(headtail)){
      return(rbind(head(obj, headtail), tail(obj, headtail)))
    }
    lookRows <- ifelse(dim(obj)[1]>4, 5, dim(obj)[1])
    lookCols <- ifelse(dim(obj)[2]>4, 5, dim(obj)[2])
    return(obj[1:lookRows,1:lookCols])
  }else if(class(obj) %in% c('dgCMatrix')){
    if(!is.null(headtail)){
      return(rbind(head(obj, headtail), tail(obj, headtail)))
    }
    message(paste0("Dimensions of ", class(obj), ": ", paste(dim(obj),collapse = ", ")))
    return(print(obj[1:4, 1:4]))
  }else if(isS4(obj)){
    if(!is.null(headtail)){
      warning('Dont know how to return head and tial for this')
    }
    return(obj)
  }else{
    library(dplyr)
    if(!is.null(headtail)){
      warning('Dont know how to return head and tial for this')
    }
    return(dplyr::glimpse(obj))
  }
}
