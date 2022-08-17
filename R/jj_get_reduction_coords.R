#' get cell embeddings together with meta data 
#'
#' Extract cell embeddings together with meta data from a Seurat, cisTopic or ArchR-object
#'
#' @param obj Seurat, cisTopic or ArchR object
#' @param redname name of the reduction, for which cell embeddings should be obtained
#' @param exact_match if TRUE, `redname` must exactly match the name of the reduction in the object
#' @keywords embeddings
#' @export
#' @examples
#' head(jj_get_reduction_coords(pbmc_small, 'tsne'))

jj_get_reduction_coords <- function(obj, redname=NULL, exact_match=TRUE){
  if(class(obj)=='Seurat'){
    dr_df_list <- lapply(obj@reductions, function(x) x@cell.embeddings[, 1:min(12, ncol(x@cell.embeddings))])
    if(!is.null(redname)){
      if(exact_match){
        dr_df_list <- dr_df_list[grepl(paste0('^',redname,'$'), names(dr_df_list))]
      }else{
        dr_df_list <- dr_df_list[grepl(redname, names(dr_df_list))]
      }
      if(length(dr_df_list)==0) stop('No reduction matches string provided in redname')
    } 
    dr_df <- do.call('cbind', dr_df_list)
    dr_df <- cbind(dr_df, obj@meta.data)
  }else if(class(obj)=='ArchRProject'){  
    dr_df_list = lapply(obj@embeddings, function(x) x@listData$df[, 1:min(12, ncol(x@listData$df))])
    if(!is.null(redname)){
      if(exact_match){
        dr_df_list <- dr_df_list[grepl(paste0('^',redname,'$'), names(dr_df_list))]
      }else{
        dr_df_list <- dr_df_list[grepl(redname, names(dr_df_list))]
      }
      if(length(dr_df_list)==0) stop('No reduction matches string provided in redname')
    } 
    dr_df <- do.call('cbind', dr_df_list)
    dr_df <- cbind(dr_df, obj@cellColData)
  }else if(class(obj)=='cisTopic'){
    dr_df_list <- lapply(obj@dr$cell, function(x) x[, 1:min(12, ncol(x))])
    if(exact_match){
      dr_df_list <- dr_df_list[grepl(paste0('^',redname,'$'), names(dr_df_list))]
    }else{
      dr_df_list <- dr_df_list[grepl(redname, names(dr_df_list))]
    }
    if(length(dr_df_list)==0) stop('No reduction matches string provided in redname')
    dr_df <- do.call('cbind', dr_df_list)
    dr_df <- cbind(dr_df, obj@cell.data)
  }
  return(dr_df)
}
