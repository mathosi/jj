#' violin plot of feature by group
#'
#' Plot a gene/numeric/categorical feature by group as violin or boxplot
#'
#' @name plot_feature_by_group 
#' @param df data.frame containing the columns passed in feature_column and group_column
#' @param seurat_obj Seurat object that contains the `assay` with the feature `gene_plot`
#' @param gene_plot Gene from the RNA matrix to plot
#' @param feature_column string of the column with feature that should be quantified
#' @param group_column string of the column with the variable used to group the feature
#' @param group_vec Vector of strings with the grouping information
#' @param type Type of plot, options are 'violin' or 'boxplot'
#' @param plot_mean Plot the mean value per group as horizontal line
#' @param plot_group_size Plot number of cells per group 
#' @param plot_zero_fraction for sparse data, plot the fraction of zero counts per group as pie
#' @param plot_cell_sample if TRUE, plot a sample of cells for each group (equal number)
#' @param order if TRUE, order the groups  by their mean value
#' @param custom_colors named vector of colors to use to fill the violins/boxplots
#' @param theme_use theme to use, default: theme_minimal()
#' @param absolute_numbers if TRUE, plot absolute counts per category instead of relative fractions per group
#' @param return_df if TRUE, instead of plotting, return the data.frame with the data
#' @param x_lab label for the groups
#' @param flip_coordinates flip coordinate system
#' @param text_size size of numbers in plots
#' @export
#' @examples
#' #plot as boxplot with additional mean, number of cells per group and cell sample (requires ggbeeswarm)
#' jj_plot_numeric_by_group(pbmc_small@meta.data, feature_column = 'nFeature_RNA', group_column = 'groups',
#'                          plot_mean = T, plot_cell_sample = T,
#'                          plot_group_size = T, type = 'boxplot')
#' #plot as violin with custom colours
#' jj_plot_numeric_by_group(pbmc_small@meta.data, feature_column = 'nFeature_RNA', group_column = 'groups',
#'                          custom_colors = c(g1='green', g2='blue'),
#'                          plot_group_size = T, type = 'violin')
#' #plot a sparse feature directly from Seurat
#' jj_plot_sparse_by_group_seurat(pbmc_small, 'CD79A', 'groups', assay='RNA', slot='data')
#' #or from a sparse matrix
#' sp_mat = GetAssayData(pbmc_small)
#' jj_plot_sparse_by_group(sp_mat, gene_plot = 'MS4A1', group_vec = pbmc_small$groups, order=T)
#' #barplot of fractions by group, eg fractions of cluster annotations per group
#' jj_plot_categorical_by_group(pbmc_small[[]], feature_column = 'RNA_snn_res.1', group_column =  'groups')
#' #or using absolute counts
#' jj_plot_categorical_by_group(pbmc_small[[]], feature_column = 'RNA_snn_res.1', group_column =  'groups', absolute_numbers = T)

#' @rdname plot_feature_by_group
#' @export
jj_plot_sparse_by_group_seurat = function(seurat_obj, gene_plot, group_column, assay=NULL, slot='data', ...){
  #jj_plot_sparse_by_group_seurat(seurat_rna, gene_plot='CD8A', group_column='cell_state', custom_colors = jj_get_colours(seurat_rna$cell_state, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/mm_scrna/scripts/colour_map.csv'))
  #jj_plot_sparse_by_group_seurat(seurat_rna, 'CD8A', 'decontX_snn_res.0.5') 
  
  if(is.null(assay)){
    assay = DefaultAssay(seurat_obj)
  }
  rna_mat = GetAssayData(seurat_obj, assay=assay, slot=slot)
  gvec = seurat_obj@meta.data[, group_column]
  jj_plot_sparse_by_group(rna_mat = rna_mat, gene_plot = gene_plot, group_vec = gvec, ...)
}

#' @rdname plot_feature_by_group
#' @export
jj_plot_sparse_by_group = function(rna_mat, gene_plot, group_vec, x_lab='Group', theme_use = theme_minimal(), 
                                plot_cell_sample=FALSE, plot_zero_fraction=TRUE, plot_mean = TRUE, plot_group_size = FALSE,
                                type='violin', custom_colors=NULL, order=FALSE, flip_coordinates=FALSE){
  #rna mat n genes (rows) * m samples (columns)
  #eg.
  #rnamat = GetAssayData(seurat_rna)
  #jj_plot_sparse_by_group(rnamat, gene_plot='CD4', group_vec=seurat_rna$sample_name, order=T)
  stopifnot(identical(ncol(rna_mat), length(group_vec)))
  stopifnot(gene_plot %in% rownames(rna_mat))
  type = match.arg(type, choices = c('violin', 'boxplot'))
  
  #get single cell expression data for the feature and add group information
  data_df = as.data.frame(Matrix::t(jj_get_feature_mat(rna_mat, gene_plot, use_features_order = T)))
  data_df$x = group_vec
  
  #mean expression per group
  mean_df = .get_mean_expressing(rna_mat, gene_plot, group_vec) %>% 
    dplyr::rename(x=cluster)
  
  #do not use scaled matrix here!
  #get number of cells per cluster and number of cells with expression > 0
  cell_exp_ct_mat = .get_fraction_expressing(rna_mat, gene_plot, group_vec) %>% 
    dplyr::rename(x=cluster)
  cell_exp_ct_mat$circle_fraction = cell_exp_ct_mat$cell_exp_ct /cell_exp_ct_mat$cell_ct *2*pi
  
  if(order){
    order_use = order(mean_df$count, decreasing = T)
    data_df$x = factor(data_df$x, levels = mean_df$x[order_use])
    cell_exp_ct_mat = cell_exp_ct_mat[order_use, ]
    mean_df = mean_df[order_use, ]
  }
  
  cell_exp_ct_mat$row_nr = 1:nrow(cell_exp_ct_mat)
  mean_df$row_nr = 1:nrow(mean_df)
  #max_y = max(data_df[, genes_plot[i]])
  if(type=='violin'){
    gg = ggplot() +
      geom_violin(data=data_df, mapping = aes_string(x = 'x', y = gene_plot, fill='x'))
  }else if(type=='boxplot'){
    gg = ggplot() +
      geom_boxplot(data=data_df, mapping = aes_string(x = 'x', y = gene_plot, fill='x'), outlier.shape = NA) 
  }
  
  if(plot_cell_sample){
    data_points_df = .get_sample_points(data_df, feature_column = gene_plot, group_column = 'x')
    gg = gg + 
      ggbeeswarm::geom_quasirandom(data = data_points_df, 
                                   mapping = aes_string(x = 'x', y = gene_plot),
                                   method = "tukeyDense", width = 0.5,
                                   alpha = 0.2, color = "black",
                                   size = 1)
  }
  
  if(plot_mean){
    gg = gg + geom_segment(data=mean_df, aes(x = row_nr-0.25, y = count, xend = row_nr+0.25, yend = count, colour = "Mean")) + 
      scale_color_manual(values=(Mean='brown4'))
  }
  
  if(plot_zero_fraction){
    library(ggforce)
    #experimental: scale_pies by group size, problem: does not reflect the true size of the population
    scale_by_size = F
    if(scale_by_size){
      min_size = 0.1
      max_size = 0.2
      max_size = max_size - min_size
      cell_exp_ct_mat$scaled_cluster_size = min_size + max_size *min_max_normalize(cell_exp_ct_mat$cell_ct)
      gg = gg + 
        geom_circle(data = cell_exp_ct_mat, aes(x0 = row_nr,  y0=-0.25, r=scaled_cluster_size), fill = 'grey') + 
        geom_arc_bar(data = cell_exp_ct_mat, aes(x0 = row_nr, y0 = -0.25, start = 0, end = circle_fraction, r0 = 0, r = scaled_cluster_size),
                     fill = 'black') + 
        coord_equal() 
    }else{
      gg = gg + 
        geom_circle(data = cell_exp_ct_mat, aes(x0 = row_nr,  y0=-0.25, r=0.1), fill = 'grey') + 
        geom_arc_bar(data = cell_exp_ct_mat, aes(x0 = row_nr, y0 = -0.25, start = 0, end = circle_fraction, r0 = 0, r = 0.1),
                     fill = 'black') + 
        coord_equal() 
    }
    
  }
  
  if(plot_group_size){
    gg = gg + 
      geom_text(data = cell_exp_ct_mat, mapping = aes(x = row_nr, y = -0.5, label = cell_ct ))
  }
  
  if(flip_coordinates){
    gg = gg + coord_flip()
  }
  
  if(!is.null(custom_colors)){
    if(!is.null(names(custom_colors)) & is.factor(data_df$x)){
      custom_colors = custom_colors[match(levels(data_df$x), names(custom_colors))]
    }
    gg = gg + scale_fill_manual(values=custom_colors)
  }
  
  gg = gg + labs(x=x_lab, fill = x_lab, colour='') + theme_use
  
  return(gg)
}

#' @rdname plot_feature_by_group
#' @export
jj_plot_numeric_by_group = function(df, feature_column, group_column, custom_colors=NULL,
                                 plot_cell_sample=FALSE, plot_mean = TRUE, plot_group_size = FALSE,
                                 theme_use = theme_minimal(), type='violin', order=FALSE, flip_coordinates=FALSE){
  #jj_plot_numeric_by_group(seurat_rna@meta.data, 'dissimilarity.score', 'patient', plot_group_size = T, flip_coordinates=T)
  type = match.arg(type, choices = c('violin', 'boxplot'))
  data_df = as.data.frame(df)[, c(group_column, feature_column)]
  
  #mean expression per group
  mean_df =  data_df %>%
    dplyr::group_by(!!rlang::sym(group_column)) %>% 
    dplyr::summarise(count=mean(!!rlang::sym(feature_column)), n=n()) %>% 
    as.data.frame 
  
  if(order){
    order_use = order(mean_df$count, decreasing = T)
    data_df[, group_column] = factor(data_df[, group_column], levels = mean_df[, group_column][order_use])
    mean_df = mean_df[order_use, ]
  }
  mean_df$row_nr = 1:nrow(mean_df)
  
  if(type=='violin'){
    gg = ggplot() +
      geom_violin(data=data_df, mapping = aes_string(x = group_column, y = feature_column, fill=group_column))
  }else if(type=='boxplot'){
    gg = ggplot() +
      geom_boxplot(data=data_df, mapping = aes_string(x = group_column, y = feature_column, fill=group_column),
                   outlier.shape = NA) 
  }
  
  if(plot_cell_sample){
    data_points_df = .get_sample_points(data_df, feature_column = feature_column, group_column = group_column)
    gg = gg + 
      ggbeeswarm::geom_quasirandom(data = data_points_df, 
                                   mapping = aes_string(x = group_column, y = feature_column),
                                   method = "tukeyDense", width = 0.5,
                                   alpha = 0.2, color = "black",
                                   size = 1)
  }
  
  if(plot_mean){
    gg = gg + geom_segment(data=mean_df, aes(x = row_nr-0.25, y = count, xend = row_nr+0.25, yend = count, colour = "Mean")) + 
      scale_color_manual(values=(Mean='brown4'))
  }
  
  if(plot_group_size){
    gg = gg + 
      geom_text(data = mean_df, mapping = aes(x = row_nr, y = -0.25, label = n ))
  }
  
  if(flip_coordinates){
    gg = gg + coord_flip()
  }
  
  if(!is.null(custom_colors)){
    gg = gg + scale_fill_manual(values=custom_colors)
  }else{
    gg = gg + scale_fill_manual(values=jj_get_jj_colours(data_df[, group_column]))
  }
  
  gg = gg + theme_use
  
  return(gg)
}


#' @rdname plot_feature_by_group
#' @export
jj_plot_categorical_by_group = function(df,
                                        feature_column, 
                                        group_column,
                                        custom_colors=NULL,
                                        absolute_numbers=FALSE, 
                                        return_df = FALSE,
                                        flip_coordinates=FALSE, 
                                        add_text=FALSE,
                                        text_size = 5,
                                        theme_use = theme_minimal()){

  summarise_fractions = function(df, summarise_by){
    summary_df = df %>% 
      dplyr::group_by(!!sym(summarise_by)) %>% 
      dplyr::mutate(total=n()) %>% 
      dplyr::group_by_all() %>% 
      dplyr::summarise(ncells= n()) %>% 
      dplyr::mutate(fraction=ncells/total) %>% 
      dplyr::ungroup() %>% 
      dplyr::arrange(!!sym(summarise_by))
    return(summary_df)
  }
  summary_df = summarise_fractions(df[, c(group_column, feature_column)],
                                   summarise_by = group_column)
  total_cells = sum(summary_df$ncells)
  summary_df$total_frac = summary_df$ncells / total_cells
  summary_df = as.data.frame(summary_df)
  
  if(return_df){
    return(summary_df)
  }
  
  if(!is.null(custom_colors)){
    summary_df[, feature_column] = factor(summary_df[, feature_column], levels = names(custom_colors))
  }
  
  if(!absolute_numbers){
    gg = ggplot(summary_df, aes_string(x=group_column, y = 'fraction', fill = feature_column)) +
      geom_bar(stat='identity', position="fill")
    if(add_text){
      gg = gg + geom_text(aes(label = paste0(round(fraction*100,2),"%")), 
                          position = position_stack(vjust = 0.5), size = text_size)
    }
  }else{
    gg = ggplot(summary_df, aes_string(x=group_column, y = 'ncells', fill = feature_column)) +
      geom_bar(stat='identity')
    if(add_text){
      gg = gg + geom_text(aes(label = ncells), 
                          position = position_stack(vjust = 0.5), size = text_size)
    }
  }
  
  if(flip_coordinates){
    gg = gg + coord_flip()
  }
  
  if(!is.null(custom_colors)){
    gg = gg + scale_fill_manual(values=custom_colors)
  }else{
    gg = gg + scale_fill_manual(values=jj_get_jj_colours(summary_df[, feature_column]))
  }
  
  gg = gg + theme_use
  
  return(gg)
}

.get_sample_points = function(df, feature_column, group_column, sample_size=500){
  df <- df %>%
    dplyr::group_by(!!rlang::sym(group_column)) %>%
    dplyr::add_count(name = "n") %>%
    dplyr::mutate(n = ifelse(n > sample_size, sample_size, n)) %>%
    dplyr::sample_n(n[1]) %>%
    dplyr::select(-n) 
  return(df)
}

.get_mean_expressing =  function(rna_mat, genes_plot, group_vec, scale_data=FALSE){
  #rna_mat genes (rows) * cells (columns)
  rna_mat = rna_mat[rownames(rna_mat) %in% genes_plot,,drop=FALSE ]  %>% Matrix::t()
  
  if(scale_data){
    rna_mat = scale(rna_mat)
  }
  expr_mat <- rna_mat %>%
    as.data.frame %>% 
    dplyr::mutate(cluster=group_vec) %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise_all(list(mean)) %>% 
    as.data.frame %>% column_to_rownames('cluster') %>%
    t %>% 
    as.data.frame %>%
    rownames_to_column('Gene') %>% 
    tidyr::pivot_longer(!Gene, names_to = 'cluster', values_to = 'count') 
  return(expr_mat)
}

.get_fraction_expressing = function(rna_mat, genes_plot, group_vec){
  #rna_mat genes (rows) * cells (columns)
  rna_mat = rna_mat[rownames(rna_mat) %in% genes_plot,,drop=FALSE ]  %>% Matrix::t()
  cell_exp_ct_mat <- rna_mat %>%
    as.data.frame %>% 
    dplyr::mutate(cluster=group_vec) %>% 
    dplyr::group_by(cluster) %>%
    dplyr::summarise_all(list(function(x) sum(x>0))) %>% 
    as.data.frame %>%
    column_to_rownames('cluster') %>%
    t %>% 
    as.data.frame %>%
    rownames_to_column('Gene') %>% 
    tidyr::pivot_longer(!Gene, names_to = 'cluster', values_to = 'cell_exp_ct')
  cell_ct_mat <- data.frame(cluster=group_vec) %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::summarise(cell_ct=n()) %>% 
    as.data.frame 
  cell_exp_ct_mat = cell_exp_ct_mat %>% 
    dplyr::left_join(cell_ct_mat, by = 'cluster')
  return(cell_exp_ct_mat)
}
