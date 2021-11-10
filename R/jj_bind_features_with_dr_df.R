#' bind features together with meta.data
#'
#' add values from a matrix/seurat object to data.frame containing metadata (dr_df)
#'
#' @param obj Seurat object or matrix of features*samples
#' @param assay assay in Seurat object to use
#' @param slot slot in assay to use
#' @param features Features to extract from obj
#' @param dr_df optional data.frame, which is joined together with the feature matrix
#' @param cap_top value passed to jj_cap_vals
#' @param cap_bottom value passed to jj_cap_vals
#' @param log10Transform if TRUE, return log10 transformed feature values, default: F
#' @keywords cap
#' @export
#' @examples
#' jj_bind_features_with_dr_df(seurat_rna, assay='RNA', slot='data', features=c('CD4','CD8A'), cap_top='q95')


jj_bind_features_with_dr_df <- function(obj, assay='RNA', slot='counts', 
                                     features, dr_df=NULL, cap_top=NULL, 
                                     cap_bottom=NULL, log10Transform=FALSE,...){
  #
  #if dr_df=NULL, just return the matrix from the seurat assay
  if(is.character(dr_df)){
    dr_df <- jj_get_reduction_coords(obj, dr_df)
  }
  if(class(obj) == 'Seurat'){
    obj = GetAssayData(obj, assay=assay, slot=slot)
  }
  goi_df <- t(jj_get_feature_mat(obj,
                             features = features))
  
  if(!all(features %in% colnames(goi_df))){
    stop(sprintf('%s not found in the assay.', paste(features[!features %in% colnames(goi_df)], collapse = ', ')))
  }
  if(log10Transform){
    message('Performing log transformation:  log10(x+0.1) + 1')
    goi_df <- apply(goi_df, 2, function(x) log10(x+0.1) + 1)
  }
  if(!is.null(cap_top) | !is.null(cap_bottom)){
    message('Capping values')
    goi_df <- apply(goi_df, 2, jj_cap_vals, cap_top=cap_top,  cap_bottom=cap_bottom)
  }
  colnames(goi_df) <- goi <- gsub('-', '_', colnames(goi_df))
  if(is.null(dr_df)){
    return(goi_df)
  }
  dr_df <- cbind(dr_df, goi_df)
  return(dr_df)
}