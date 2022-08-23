#' find associations of embedding scores with meta data
#'
#' jj_metadata_association_heatmap: Plot heatmap to find metadata associations with a set of scores, e.g. from PCA
#' jj_association_heatmap: Plot mean or scaled mean per group for all columns in scores_mat, good to further understand single metadata associations with scores.
#' Further arguments can be passed to ComplexHeatmap
#' jj_test_associations: Function to test association of metadata with a group of interest. Run logistic regressions with group1/group2 as outcome and a column in the metadata_df as predictor
#' Repeat for all columns and sort by p-value and deviance against a null model (group ~ 1). This can reveal metadata variables enriched for the cluster of interest (group1)
#' @name jj_association_heatmap
#' @param scores_mat numeric matrix of scores (eg pc scores)
#' @param metadata_df data.frame/matrix with the metadata for each cell
#' @param group_vec vector of length = nrow(scores_mat) with group annotation
#' @param cor_cutoff Only show correlations with absolute value above this cutoff for better overview
#' @param plotAssocFor Option to specify subset of metadata_mat for which associations should be calulated
#' @param categorical_skip Skip categorical metadata columns with more than this number of distinct values
#' @param group1 cell names of group, for which jj_test_associations should test associations. Used to select the subset of metadata_df (which needs to contain rownames with the cell names
#' @param group2 group1 is compared against group2. If group2=NULL (default), all cells in the metadata_df not part of group1 are defined as group2
#' @export
#' @examples
#' #e.g. PC2 has positive correlation with the grouping from RNA_snn_res.1
#' pca_scores = Embeddings(pbmc_small, 'pca')[,1:10]
#' jj_metadata_association_heatmap(score_mat = pca_scores, metadata_df = pbmc_small[[]],
#'                                 cor_cutoff = 0.2)
#' # Groups 1 and 2 tend to have high scores and group 0 has low scores for PC2
#' jj_association_heatmap(scores_mat = pca_scores, group_vec = pbmc_small$RNA_snn_res.1)
#' # confirm by plotting the data
#' jj_plot_features(pbmc_small, reduction = 'pca', meta_features = 'RNA_snn_res.1', pt.size = 2)
#' #test association of meta data with one selected group of cells
#' group1 = colnames(pbmc_small)[pbmc_small$RNA_snn_res.1 == '0']
#' jj_test_associations(metadata_df = pbmc_small@meta.data, group1 = group1)

jj_metadata_association_heatmap = function(score_mat, metadata_df, plotAssocFor=NULL, cor_cutoff=0.2, categorical_skip=30, ...){
  library(ComplexHeatmap)
  if(is.null(plotAssocFor)){
    plotAssocFor = colnames(metadata_df)
    plotAssocFor = plotAssocFor[apply(metadata_df,2, function(x) length(unique(x)) > 1 & length(unique(x)) < nrow(metadata_df))]
  } 
  #numericCols <- colnames(metadata_df[plotAssocFor])[apply(metadata_df[,plotAssocFor],2, function(x){all(grepl(x, pattern= "^[0-9]+(\\.[0-9]+)*$"))})]
  categoricalCols <- suppressWarnings(plotAssocFor[apply(metadata_df[,plotAssocFor],2,function(x) all(is.na(as.numeric(x))))])
  numericCols <- plotAssocFor[!plotAssocFor %in% categoricalCols]
  
  categoricalCols = categoricalCols[!sapply(categoricalCols, function(x) length(unique(metadata_df[, x])) > categorical_skip)]
  
  assoc_df = calculateAssociation(score_mat, metadata_df, plotAssocFor = categoricalCols)
  numeric_mat = apply(metadata_df[, colnames(metadata_df) %in% numericCols],2, as.numeric)
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

#' @export
jj_test_associations = function(metadata_df, group1, group2=NULL){

  stopifnot(all(group1 %in% rownames(metadata_df)))
  if(!is.null(group2)){
    stopifnot(all(group2 %in% rownames(metadata_df)))
  }else{
    group2 = rownames(metadata_df)[!rownames(metadata_df) %in% group1]
  }
  metadata_df$association_group = ifelse(rownames(metadata_df) %in% group1, 0, ifelse(rownames(metadata_df) %in% group2, 1, 2))
  metadata_df = metadata_df[metadata_df$association_group %in% c(0, 1), ]
  pval_vec = dev_vec = vector()
  for(i in 1:(ncol(metadata_df)-1)){
    message(i)
    metadata_df_use = metadata_df[!is.na(metadata_df[, i]), ]
    metadata_df_use$ass_var = metadata_df_use[, i]
    lun = length(unique(metadata_df_use$ass_var))
    if(!is.numeric(metadata_df_use$ass_var) & (lun > 100 | lun == 1)){
      message(sprintf('Skipping %s since number of levels is too high: %i', colnames(metadata_df)[i], lun))
      next
    }
    #if(class(metadata_df$ass_var) == 'factor'){
    #  metadata_df$ass_var = as.character(metadata_df$ass_var)
    #}
    #fit logistic regression model
    glm.fit.naive <- glm(association_group ~ 1, data = metadata_df_use, family = binomial)
    #print(summary(glm.fit.naive))
    glm.fit <- glm(association_group ~ ass_var, data = metadata_df_use, family = binomial)
    #print(summary(glm.fit))
    res = anova(glm.fit.naive, glm.fit, test="LRT")
    pval_vec[colnames(metadata_df)[i]] = res$`Pr(>Chi)`[2]
    dev_vec[colnames(metadata_df)[i]] = res$Deviance[2]
  }
  pval_vec = p.adjust(pval_vec)
  res_df = data.frame(variable=names(pval_vec), deviance = dev_vec, p_adj = pval_vec)
  res_df = res_df[with(res_df, order(p_adj, -deviance)), ]
  return(res_df)
}
