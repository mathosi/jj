#' summarise values as mean, mode, or sum by a group vector 
#'
#' jj_summarize_sparse_mat calculates the sum or mean of values for each group in a sparse matrix.
#' jj_summarize_vector summarizes a vector based on a group vector
#' jj_summarize_dataframe summarizes each column in a data.frame based on a group vector
#'
#' @name summarise_vals
#' @param summarize_obj sparse matrix, vector, or data.frame
#' @param summarize_by_vec vector with group annotation, has to have length equal to ncol(sparse_mat)/nrow(data.frame)/length(vector)
#' @param method method used to summarize. Only mode available for non-numeric vectors, sum/mean for sparse matrix and mode/mean/median for numeric vectors
#' @param return_matrix return as matrix instead of sparse matrix
#' @param order order the summarized vector by the groups
#' @param return_vec return only summarized vector instead of data.frame with summarized vector and groups
#' @keywords summarize
#' @export
#' @examples
#' df = data.frame(a=seq(1,5), b=c('a','b','b','b','a'), d=c(3,0,0,1,5))
#' jj_summarize_vector(df$d, df$b, method='mean')
#' jj_summarize_dataframe(df, df$b)

#' @rdname summarise_vals
#' @export
jj_summarize_sparse_mat <- function(summarize_obj, summarize_by_vec, method='mean', return_matrix=TRUE){ 
  #also works for normal matrices
  stopifnot(is.vector(summarize_by_vec) | is.factor(summarize_by_vec))
  #stopifnot(!is.factor(summarize_by_vec))
  if(!identical(ncol(summarize_obj), length(summarize_by_vec))){
    stop('Number of columns in the assay and length of group vector must be identical')
  }
  method = match.arg(method, choices=c('mean', 'sum', 'sd'))
  if(method == 'sd'){
    sd_mat = as(jj_initialize_df(ncol = length(unique(summarize_by_vec)), nrow = nrow(summarize_obj), init = 0, return_matrix = T,
                                 col.names = unique(summarize_by_vec), row.names= rownames(summarize_obj)), 'dgCMatrix')
    for(i in unique(summarize_by_vec)){
      sd_mat[, i] = proxyC::rowSds(summarize_obj[, summarize_by_vec == i]) #also considers all . (=0)
    }
    summarize_obj = sd_mat
  }else{
    #be careful: function 'mean' ignores . completely
    ##summarize_obj <- Matrix.utils::aggregate.Matrix(summarize_obj, groupings = summarize_by_vec, fun = 'sum')
    summarize_obj = Matrix::tcrossprod(Matrix::fac2sparse(summarize_by_vec), summarize_obj)
    groups_in_vec <- rownames(summarize_obj)
    #divide by n cells per group if mean should be returned
    if(method=='mean'){
      groups_freq <- table(summarize_by_vec)
      summarize_obj <- summarize_obj /  as.vector(groups_freq)[match( groups_in_vec, names(groups_freq))]
    }
    summarize_obj <- Matrix::t(summarize_obj)
  }
  summarize_obj = summarize_obj[, gtools::mixedorder(colnames(summarize_obj))]
  
  if(return_matrix){
    return(as.matrix(summarize_obj))
  }
  return(summarize_obj)
}

#' @rdname summarise_vals
#' @export
jj_summarize_vector = function(summarize_obj, summarize_by_vec, method='mean', order=T, return_vec=F){
  stopifnot((is.vector(summarize_obj) | is.factor(summarize_obj)) & (is.vector(summarize_by_vec) | is.factor(summarize_by_vec)))
  options(dplyr.summarise.inform = FALSE)
  method = match.arg(method, choices = c('mode','mean','median','count'))
  df = data.frame(a=summarize_obj, b=summarize_by_vec)
  if(method=='mode'){
    df <- df %>% 
      dplyr::group_by_all() %>% 
      dplyr::summarise(nr_cells=n()) %>% 
      dplyr::ungroup() %>%  
      dplyr::group_by(b) %>% 
      dplyr::filter(nr_cells == max(nr_cells)) %>% 
      .[!duplicated(.$b), ] %>% 
      dplyr::select(-nr_cells)
  }else if(method=='mean'){
    df <- df %>% 
      dplyr::group_by(b) %>% 
      dplyr::summarise(a=mean(a, na.rm=T)) %>% 
      dplyr::select(a, b)
  }else if(method == 'median'){
    df <- df %>% 
      dplyr::group_by(b) %>% 
      dplyr::summarise(a=median(a, na.rm=T)) %>% 
      dplyr::select(a, b)
  }else if(method == 'count'){
    df <- df %>% 
      dplyr::group_by_all() %>% 
      dplyr::summarise(count=n()) 
    return(df)
  }
  if(order){
    df = df %>% dplyr::arrange(b)
  }
  if(return_vec){
    a = df$a
    names(a) = df$b
    return(a)
  }
  return(df)
}

#' @rdname summarise_vals
#' @export
jj_summarize_dataframe = function(summarize_obj, summarize_by_vec, method = 'mean'){
  stopifnot(is.vector(summarize_by_vec) | is.factor(summarize_by_vec))
  summarize_df = as.data.frame(summarize_obj)
  stopifnot(identical(nrow(summarize_df), length(summarize_by_vec)))
  sdf = jj_initialize_df(ncol = ncol(summarize_df), 
                         nrow = length(unique(summarize_by_vec)),
                         init = NA, 
                         col.names = colnames(summarize_df),
                         row.names =  dplyr::pull(dplyr::arrange(data.frame(b= unique(summarize_by_vec)), b), b)) 
  
  for(i in 1:ncol(summarize_df)){
    if(class(summarize_df[, i]) %in% c('numeric','integer')){
      message(sprintf('%i/%i: Summarize %s by %s', i, ncol(summarize_df), colnames(summarize_df[i]), method))
      svec = jj_summarize_vector(summarize_df[, i], summarize_by_vec, method=method, order = T, return_vec = T)
      stopifnot(identical(names(svec), rownames(sdf)))
      sdf[, i]  = svec
    }else{
      message(sprintf('%i/%i: Summarize %s by mode', i, ncol(summarize_df), colnames(summarize_df[i])))
      svec = jj_summarize_vector(summarize_df[, i], summarize_by_vec, method='mode', order = T, return_vec = T)
      stopifnot(identical(names(svec), rownames(sdf)))
      sdf[, i] = svec
    }
  }
  return(sdf)
}
