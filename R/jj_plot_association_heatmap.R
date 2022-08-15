#' find associations of embedding scores with meta data
#'
#' jj_metadata_association_heatmap: Plot heatmap to find metadata associations with a set of scores, e.g. from PCA
#' jj_association_heatmap: Plot mean or scaled mean per group for all columns in scores_mat, good to further understand single metadata associations with scores.
#' Further arguments can be passed to ComplexHeatmap
#' @name jj_association_heatmap
#' @param scores_mat numeric matrix of scores (eg pc scores)
#' @param metadat_mat data.frame/matrix with the metadata for each cell
#' @param group_vec vector of length = nrow(scores_mat) with group annotation
#' @param cor_cutoff Only show correlations with absolute value above this cutoff for better overview
#' @param plotAssocFor Option to specify subset of metadata_mat for which associations should be calulated
#' @param categorical_skip Skip categorical metadata columns with more than this number of distinct values
#' @export
#' @examples
#' #e.g. PC2 has positive correlation with the grouping from RNA_snn_res.1
#' pca_scores = Embeddings(pbmc_small, 'pca')[,1:10]
#' jj_metadata_association_heatmap(score_mat = pca_scores, metadat_mat = pbmc_small[[]],
#'                                 cor_cutoff = 0.2)
#' # Groups 1 and 2 tend to have high scores and group 0 has low scores for PC2
#' jj_association_heatmap(scores_mat = pca_scores, group_vec = pbmc_small$RNA_snn_res.1)
#' # confirm by plotting the data
#' jj_plot_features(pbmc_small, reduction = 'pca', meta_features = 'RNA_snn_res.1', pt.size = 2)

jj_metadata_association_heatmap = function(score_mat, metadat_mat, plotAssocFor=NULL, cor_cutoff=0.2, categorical_skip=30, ...){
  library(ComplexHeatmap)
  if(is.null(plotAssocFor)){
    plotAssocFor = colnames(metadat_mat)
    plotAssocFor = plotAssocFor[apply(metadat_mat,2, function(x) length(unique(x)) > 1 & length(unique(x)) < nrow(metadat_mat))]
  } 
  #numericCols <- colnames(metadat_mat[plotAssocFor])[apply(metadat_mat[,plotAssocFor],2, function(x){all(grepl(x, pattern= "^[0-9]+(\\.[0-9]+)*$"))})]
  categoricalCols <- suppressWarnings(plotAssocFor[apply(metadat_mat[,plotAssocFor],2,function(x) all(is.na(as.numeric(x))))])
  numericCols <- plotAssocFor[!plotAssocFor %in% categoricalCols]
  
  categoricalCols = categoricalCols[!sapply(categoricalCols, function(x) length(unique(metadat_mat[, x])) > categorical_skip)]
  
  assoc_df = calculateAssociation(score_mat, metadat_mat, plotAssocFor = categoricalCols)
  numeric_mat = apply(metadat_mat[, colnames(metadat_mat) %in% numericCols],2, as.numeric)
  assoc_df2 = cor(score_mat, numeric_mat)
  total_df = cbind(assoc_df, assoc_df2)
  total_df2 = round(total_df,2)
  total_df2[!abs(total_df) > cor_cutoff] = ''
  h1 <- Heatmap(total_df, column_title = "Pearson correlation numeric variable / lmfit (categorical) vs observed PC score",
                col =  circlize::colorRamp2(c(-1, -cor_cutoff, 0, cor_cutoff, 1), c("darkblue","white", "white","white", "red")),
                name="Pcor",
                cell_fun = function(j, i, x, y, width, height, fill) {
                  grid.text(sprintf("%s", total_df2[i, j]), x, y, gp = gpar(fontsize = 10))},
                #heatmap_height = unit(12,"cm"),
                na_col = "grey30", ...)
  return(h1)
}

#' @export
jj_association_heatmap <- function(scores_mat, group_vec, do_scale=T, ...){
  stopifnot(identical(nrow(scores_mat), length(group_vec)))
  summary_df = jj_summarize_dataframe(scores_mat, group_vec)
  title_use = 'mean'
  if(do_scale){
    summary_df = scale(summary_df)
    title_use ='scaled mean'
  }
  ComplexHeatmap::Heatmap(summary_df, name=title_use, ...)
}

calculateAssociation <- function(scores_mat, pheno_df, plotAssocFor = NULL){
  #return associations for categorical data with continuous data
  # linear model predicts the pc score for each cell based on its group (=mean pc score for that group is always predicted): if the mean pc score is quite different between the groups, then there is a high correlation
  # correlation between real pc score and predicted pc score (=mean)
  
  if(is.null(plotAssocFor)) plotAssocFor = colnames(pheno_df)
  phenoDf = as.data.frame(cbind(scores_mat, pheno_df))
  largeCoeffDf <- jj_initialize_df(nrow = ncol(scores_mat), ncol = 0, row.names = colnames(scores_mat),
                                  init = NA)
  for(predictor in plotAssocFor){
    #fitted values are just mean LMC props per level in predictor var
    cat(paste0(predictor,"\n"))
    corrVec <- rep(NA,ncol(scores_mat))
    for(pc in 1:ncol(scores_mat)){
      pc_name <- colnames(scores_mat)[pc]
      if(length(unique((phenoDf[,predictor])))>1){
        model <- lm(eval(parse(text=pc_name)) ~ eval(parse(text=predictor)), data= phenoDf, na.action = "na.omit")
        #only get the proportion values for the non-NA positions of the predictor var
        origLMCProp <- phenoDf[names(model$fitted.values), pc_name]
        corrVec[pc] <- cor(origLMCProp, model$fitted.values)
      }
    }
    largeCoeffDf[,predictor] <- corrVec
  }
  return(largeCoeffDf)
}
