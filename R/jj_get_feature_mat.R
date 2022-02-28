#' extract expression values for features of interest from feature matrix
#'
#' extracts gene expression values for a given vector of features of interest
#'
#' @param feature_mat matrix with n features (rows) and m observations (columns)
#' @param features features for which the values should be extracted from the matrix
#' @param fail_on_missing return error if not all features are found
#' @param use_features_order order the rows in the matrix by the order of `features`, otherwise just leave the order as in the supplied `feature_mat`
#' @keywords get expression, expression, extract expression, features
#' @export
#' @examples
#' mat_use = matrix(1:20, nrow=10, byrow=T)
#' rownames(mat_use) = LETTERS[1:10]
#' jj_get_features(mat_use, c('D','J','A'))
#' jj_get_features(mat_use, c('D','J','A', 'Z'),  use_features_order=F)

jj_get_features = function(feature_mat, features, fail_on_missing=FALSE, use_features_order=TRUE, verbose=TRUE){
  is.any.matrix = function(m){ is.matrix(m) | is.data.frame(m) | is(m, 'sparseMatrix') }
  
  stopifnot(is.character(features))
  stopifnot(is.any.matrix(feature_mat))
  
  features = unique(features)
  lf = length(features)
  
  #show which features are found
  feat_avail = rownames(feature_mat)
  features_found = features[features %in% feat_avail]
  lff = length(features_found)
  
  if(lff == 0){
    stop('None of the features is available in the matrix')
  }
  
  if(!identical(lf, lff)){
    features_missing <- features[!features %in% feat_avail]
    if(verbose){
      lfm <- length(features_missing)
      if(lfm>100){
        features_missing <- c(features_missing[1:100], "... ")
      }
      if(lfm>0){
        message(sprintf("\n\nNo expression values found for %i out of %i features:\n", lfm, lf))
        message(features_missing)
      }
    }
    if(fail_on_missing){
      stop('Not all features are available in the matrix')
    }
  }else{
    if(verbose){
      message("\n\nExpression values found for all ", lf, " features\n")
    }
  }
  
  feature_mask = feat_avail %in% features_found
  feature_mat <- feature_mat[feature_mask,, drop=FALSE]
  if(use_features_order){
    feature_mat <- feature_mat[match(features_found, rownames(feature_mat)),, drop=FALSE ]
  }
  return(feature_mat)
}

# jj_get_feature_mat <- function(feature_mat, features=NULL, verbose=TRUE, ...){
#   library(dplyr)
#   if(is.null(features)){
#     features <- unique(rownames(feature_mat))
#   }
#   geneNames <- rownames(feature_mat)
#   featuresMask <- geneNames %in% features
#   feature_mat <- feature_mat[featuresMask,, drop=FALSE]
#   geneNamesNew <- geneNames[featuresMask]
#   if(any(duplicated(geneNamesNew))){
#     if(verbose){
#       warning("Expression matrix contains duplicated rownames. Summing up counts from rows with same name.")
#     }
#     feature_mat <- feature_mat %>%
#       as.data.frame %>%
#       dplyr::group_by(summarizedGenes=!!geneNamesNew) %>%
#       dplyr::summarise_all(list(sum)) %>%
#       as.data.frame %>%
#       tibble::column_to_rownames('summarizedGenes')
#   }else{
#     feature_mat <- as.data.frame(feature_mat)
#     rownames(feature_mat) <- geneNamesNew
#   }
# 
#   #exclude all features, which cannot be found in the gene expression data
#   mf <- features[!features %in% geneNames]
#   if(verbose){
#     lengthmf <- length(mf)
#     if(lengthmf>100){
#       mf <- c(mf[1:100], "... ")
#     }
#     if(lengthmf>0){
#       message(paste0("\n\nNo expression values found for ",
#                  lengthmf, " out of ", length(features)," genes:\n\n"))
#       message(mf)
#     }else{
#       message("\n\nExpression values found for all ", length(features), " genes.\n")
#     }
#   }
#   features <- features[features %in% geneNames]
#   if(is.null(features)){
#     return(NULL)
#   }else{
#     fmat <- feature_mat[rownames(feature_mat) %in% features,, drop=FALSE]
#   }
#   fmat <- fmat[match(features, rownames(fmat)), ]
#   return(fmat)
# }