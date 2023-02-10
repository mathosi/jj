#' Cap values at the top/bottom of a vector by a threshold
#'
#' Useful for plotting when there are few outlier values stretching the colour scale
#'
#' @param vec numeric vector
#' @param cap_top cap values above this threshold, either a single number or a quantile as string, eg. 'q95' for quantile 95
#' @param cap_bottom same as cap_top, but cap values below this threshold
#' @keywords cap
#' @export
#' @examples
#' jj_cap_vals(seq(0,100), cap_top='q95', cap_bottom=10)


jj_cap_vals <- function(vec, cap_top='q99', cap_bottom=NULL){
  if(!is.numeric(vec)){
    message('Converting vector to numeric before capping.')
    vec <- as.numeric(as.character(vec))
  }
  if(!is.numeric(cap_top) & !is.null(cap_top)){
    if(grepl('top', cap_top)){
      cap_top <- as.numeric(gsub('top', '', cap_top))
      cap_top <- sort(vec, decreasing = TRUE)[cap_top]
    }else if(grepl('q', cap_top)){
      cap_top <- quantile(vec, as.numeric(gsub('q([0-9]{1,2})([0-9]*)', '\\1.\\2', cap_top))/100, na.rm=TRUE)
    }else if(cap_top=='auto'){
      q_use = 0.95
      cap_top <- quantile(vec, q_use, , na.rm=TRUE)
      vec2 = vec
      vec2[vec2 > cap_top] <- cap_top
      if(n_distinct(vec2) < 5){
        q_use = 0.99
        cap_top <- quantile(vec, q_use, na.rm=TRUE)
        vec2 = vec
        vec2[vec2 > cap_top] <- cap_top
        if(n_distinct(vec2) < 5) cap_top = NULL; q_use='not capped'
      }
      message('cap_top: ', q_use)
      #if(cap_top==0) cap_top = NULL
    }
  }
  
  if(!is.numeric(cap_bottom) & !is.null(cap_bottom)){
    if(grepl('top', cap_bottom)){
      cap_bottom <- as.numeric(gsub('top', '', cap_bottom))
      cap_bottom <- sort(vec, decreasing = FALSE)[cap_bottom]
    }else{
      cap_bottom <- quantile(vec, as.numeric(gsub('q', '', cap_bottom))/100, na.rm=TRUE)
    }
  }
  
  if(!is.null(cap_top)){
    vec[vec > cap_top] <- cap_top
  }
  if(!is.null(cap_bottom)){
    vec[vec < cap_bottom] <- cap_bottom
  }
  return(vec)
}

