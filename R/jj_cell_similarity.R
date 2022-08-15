#' encirle fraction of points that contain x percent of the density mass
#'
#' Option 1: plot correlations between clusters of cells based on values from one matrix (eg. scaled scRNA matrix with cells from 3 different clusters)
#' Option 2: plot correlations between clusters of cells based on common features from two matrices (eg. scaled scRNA matrix with 3 groups and scaled gene activity matrix with 6 groups)
#' Additional arguments can be passed to Heatmap function such as cluster_rows=FALSE, name = 'Spearman_cor', row_title = row_title, column_title = 'mat_b_group'
#' 
#' @param scaled_mat_a scaled gene expression/activity matrices with cells as columns and rownames (genes)
#' @param group_a vector of the same length as columns in scaled_mat_a that defines the groups to summarize by
#' @param scaled_mat_b optional second scaled matrix (for option 2)
#' @param group_b vector of the same length as columns in scaled_mat_b that defines the groups to summarize by
#' @param plot_complex_heatmap if TRUE, plot the correlation matrix as ComplexHeatmap instead of returning it
#' @param label if TRUE, add text of Spearman correlation for each comparison
#' @param plot_triangular if TRUE, plot a triangular heatmap without clustering (only for option 1)
#' @param rname optional row title
#' @return Returns ComplexHeatmap of Spearman correlations between the groups
#' @export
#' @examples
#' 
#' #Option 1: Compare correlation of graph-based clusters based on RNA counts
#' pbmc_small = ScaleData(pbmc_small, features = rownames(pbmc_small))
#' scaled_rna = GetAssayData(pbmc_small, assay='RNA', slot='scale.data')
#' jj_cluster_correlation_heatmap(scaled_mat_a = scaled_rna,
#'                               group_a = as.character(pbmc_small$RNA_snn_res.1))
#' jj_cluster_correlation_heatmap(scaled_mat_a = scaled_rna,
#'                               group_a = as.character(pbmc_small$RNA_snn_res.1),
#'                               plot_triangular=TRUE)
#' #Option2: Dummy example, the idea is to compare clusters from a second dataset                              
#' pbmc_subset = pbmc_small[sample(1:nrow(pbmc_small), 100, replace = F), sample(1:ncol(pbmc_small), 50, replace = F)]
#' pbmc_subset = ScaleData(pbmc_subset, features = rownames(pbmc_subset))
#' scaled_rna2 = GetAssayData(pbmc_subset, assay='RNA', slot='scale.data')
#' jj_cluster_correlation_heatmap(scaled_mat_a = scaled_rna,
#'                               group_a = as.character(pbmc_small$RNA_snn_res.1),
#'                               scaled_mat_b = scaled_rna2,
#'                               group_b = as.character(pbmc_subset$RNA_snn_res.1))

jj_cluster_correlation_heatmap = function(scaled_mat_a, group_a, scaled_mat_b=NULL,
                                          group_b=NULL, plot_complex_heatmap=TRUE,
                                          label=TRUE, plot_triangular=FALSE, ...){
  library(Matrix)
  library(Matrix.utils)
  if(!is.null(scaled_mat_b)){
    if(!identical(rownames(scaled_mat_a), rownames(scaled_mat_b))){
      message('Rownames of matrices are not identical. Trying to match')
      common_genes = intersect(rownames(scaled_mat_a), rownames(scaled_mat_b))
      if(length(common_genes)== 0){stop('No common rownames found.')}
      scaled_mat_a = scaled_mat_a[rownames(scaled_mat_a) %in% common_genes,]
      scaled_mat_b = scaled_mat_b[rownames(scaled_mat_b) %in% common_genes,]
      scaled_mat_b= scaled_mat_b[match(rownames(scaled_mat_a), rownames(scaled_mat_b)), ]
      stopifnot(identical(rownames(scaled_mat_a), rownames(scaled_mat_b)))
    }
    assay_a = jj_summarize_sparse_mat(scaled_mat_a, group_a, method = 'mean')
    assay_b = jj_summarize_sparse_mat(scaled_mat_b, group_b, method = 'mean')
    cor_mat = cor(assay_a, assay_b, method='spearman') #rows are a clusters, columns are b clusters
  }else{
    assay_a = jj_summarize_sparse_mat(scaled_mat_a, group_a, method = 'mean')
    cor_mat = cor(assay_a[complete.cases(assay_a), ], method='spearman') #rows are a clusters, columns are b clusters
  }
  
  if(plot_complex_heatmap){
    library(ComplexHeatmap)
    if(label){
      if(plot_triangular & is.null(scaled_mat_b)){
        ht = triangular_heatmap(cor_mat, ...)
      }else{
        ht = Heatmap(cor_mat, 
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(sprintf("%.2f", cor_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                     }, ...)
      }
    }else{
      ht = Heatmap(cor_mat, ...)
    }
    return(ht)
  }
  return(cor_mat)
}

triangular_heatmap = function(mat, col_fun = NULL, text_format = '%.2f', ...){
  
  if(is.null(col_fun)){
    diag(mat) = 0
    col_fun = circlize::colorRamp2(c(min(mat), 0, max(mat)), c("green", "white", "red"))
  }
  
  ht = Heatmap(mat, col=col_fun,
               cluster_rows = FALSE, cluster_columns = FALSE, 
               rect_gp = gpar(type = "none"), 
               show_column_names = F, show_row_names = F,
               cell_fun = function(j, i, x, y, w, h, col) {
                 if(i > j) {
                   #NULL
                   #grid.rect(x, y, w, h, gp = gpar(fill = 'white', col='white', lwd=2, transparency = 0))
                 } else if(j == i) {
                   #grid.rect(x, y, w, h, gp = gpar(fill = 'white'))
                   if(!j %in% c(1)){
                     grid.text(sprintf("%s", colnames(mat)[j]), x, y+0.4*h, rot=90, just = 'right')
                   }
                 } else {
                   grid.rect(x, y, w, h, gp = gpar(fill = col))
                   grid.text(sprintf(text_format, mat[i, j]), x, y)
                 }
                 #grid.rect(x, y, w, h, gp = gpar(fill = NA, col = "black"))
               }, ...)
  
  rnames = rownames(mat)
  rnames[length(rnames)] = ''
  ht = ht  + rowAnnotation(labels = anno_text(rnames, which = "row"), 
                           width = max(grobWidth(textGrob(rownames(mat))))) # you need to calculate the width by hand
  return(ht)
}


# cell_similarity = function(seurat_obj, group_column, assay, slot, sample_size = 100){
#   #compute cell similarities between cells in different clusters (seurat subset of similar cell types should be used)
#   #TODO: include correction for peak count (only sample cells with similar peak count)
#   ass_dgc = GetAssayData(seurat_obj, assay=assay, slot=slot)
#   group_col = seurat_obj@meta.data[, group_column]
#   #peak_counts = Matrix::colSums(ass_dgc)
#   comp_assay = msInitializeDf(colNr = n_distinct(group_col),
#                               rowNr = nrow(ass_dgc),
#                               initVal = NA, 
#                               col.names = unique(group_col))
#   if(sample_size == 'all'){
#     comp_assay = jj_summarize_sparse_mat(ass_dgc, as.character(group_col), method='sum')
#   }else{
#     for(i in unique(group_col)){
#       cells_use = sample(which(group_col == i), sample_size, replace = T)
#       comp_assay[, i] = Matrix::rowSums(ass_dgc[, cells_use])
#     }
#   }
#   
#   cos_similariy_mat = as.data.frame(wordspace::dist.matrix(t(comp_assay), method='cosine', convert=F))
#   return(cos_similariy_mat)
#   #Heatmap(cell_sim)
# }
#}

