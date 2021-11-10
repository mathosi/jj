#' extract expression values for features of interest from feature matrix
#'
#' extracts gene expression values for a given vector of features of interest
#'
#' @param feature_mat matrix with n features (rows) and m observations (columns)
#' @param features features for which the values should be extracted from the matrix
#' @keywords get expression, expression, extract expression, features
#' @export
#' @examples
#'

jj_get_feature_mat <- function(feature_mat, features=NULL, verbose=TRUE, ...){
  library(dplyr)
  if(is.null(features)){
    features <- unique(rownames(feature_mat))
  }
  geneNames <- rownames(feature_mat)
  featuresMask <- geneNames %in% features
  feature_mat <- feature_mat[featuresMask,, drop=FALSE]
  geneNamesNew <- geneNames[featuresMask]
  if(any(duplicated(geneNamesNew))){
    if(verbose){
      warning("Expression matrix contains duplicated rownames. Summing up counts from rows with same name.")
    }
    feature_mat <- feature_mat %>%
      as.data.frame %>%
      dplyr::group_by(summarizedGenes=!!geneNamesNew) %>%
      dplyr::summarise_all(list(sum)) %>%
      as.data.frame %>%
      tibble::column_to_rownames('summarizedGenes')
  }else{
    feature_mat <- as.data.frame(feature_mat)
    rownames(feature_mat) <- geneNamesNew
  }

  #exclude all features, which cannot be found in the gene expression data
  notFound <- features[!features %in% geneNames]
  if(verbose){
    lengthNotFound <- length(notFound)
    if(lengthNotFound>100){
      notFound <- c(notFound[1:100], "... ")
    }
    if(lengthNotFound>0){
      message(paste0("\n\nNo expression values found for ",
                 lengthNotFound, " out of ", length(features)," genes:\n\n"))
      message(notFound)
    }else{
      message("\n\nExpression values found for all ", length(features), " genes.\n")
    }
  }
  features <- features[features %in% geneNames]
  if(is.null(features)){
    return(NULL)
  }else{
    fmat <- feature_mat[rownames(feature_mat) %in% features,, drop=FALSE]
  }
  fmat <- fmat[match(features, rownames(fmat)), ]
  return(fmat)
}
