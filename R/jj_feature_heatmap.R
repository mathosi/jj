#' feature heatmap
#'
#' 
#' @name jj_feature_heatmap
#' @param obj Seurat object or matrix with features in rows and cells in columns
#' @param features_use Features to plot from the Seurat active assay or the matrix
#' @param group_vec vector with group annotation, has to have length equal to ncol(obj)
#' @param scale_data scale feature-wise (i.e. the columns of the heatmap)
#' @param show_top_annot if TRUE and features_use is a named vector, the names are used as top annotation
#' @param row_annot either TRUE/FALSE to show colour bar by group from group_vec or a HeatmapAnnotation object
#' @param transpose if TRUE, transpose the heatmap and show features as rows
#' @param column_colors list with colour mapping for the column labels
#' @param group_annot_df optional data.frame with group annotations to plot for the columns in jj_plot_dotplot
#' @export
#' @examples
#' jj_plot_heatmap(pbmc_small, features_use = c('CD8A','KLRG1','CD79A','CD79B','CD3E'), group_vec = pbmc_small$groups, scale_data=F)
#' #pass named vector to features_use to show column annotation
#' jj_plot_heatmap(GetAssayData(pbmc_small), features_use = c(T='CD8A',T='KLRG1',B='CD79A',B='CD79B',T='CD3E', myeloid = 'CD14'), group_vec = pbmc_small$groups, row_annot = T)
#' #when transposing the heatmap, the row_annot is shown as top annotation
#' jj_plot_heatmap(GetAssayData(pbmc_small), features_use = c('CD8A','KLRG1','CD79A','CD79B','CD3E'), group_vec = pbmc_small$groups, transpose = T, row_annot = T)
#' jj_plot_heatmap(GetAssayData(pbmc_small), features_use = c('CD8A','KLRG1','CD79A','CD79B','CD3E'), group_vec = pbmc_small$groups, transpose = T, row_annot = HeatmapAnnotation(myannot = c('Group1', 'Group2'), col = list(myannot = c(Group1='green', Group2='cyan'))))
#' jj_plot_dotplot(GetAssayData(pbmc_small), features_use = c('CD8A','KLRG1','CD79A','CD79B','CD3E'), group_vec = pbmc_small$groups)



#' @rdname feature_heatmap
#' @export
jj_plot_heatmap = function(obj, features_use, group_vec, scale_data=TRUE, show_top_annot=TRUE, row_annot = NULL,
                           return_matrix=FALSE, cluster_rows=T, cluster_columns=T, transpose = F, ...){
  #create heatmap of scaled mean normalized expression for `features_use` per group in `group_vec`
  #scaling ensures that color scale is comparable between genes
  library(ComplexHeatmap)
  if('Seurat' %in% class(obj)){
    heatmap_mat = Matrix::t(jj_bind_features_with_dr_df(obj, slot='data', features=features_use))
  }else{
    message('Class of obj is not Seurat, assuming it is a matrix of normalized counts (features in rows, cells in columns)')
    heatmap_mat = obj[rownames(obj) %in% features_use, ]
    stopifnot(nrow(heatmap_mat)>1 & ncol(heatmap_mat) > 1)
    heatmap_mat = heatmap_mat[match(features_use, rownames(heatmap_mat)), ]
  }
  
  heatmap_mat = jj_summarize_sparse_mat(heatmap_mat,summarize_by_vec = group_vec,
                                        method = 'mean')
  heatmap_mat = base::t(heatmap_mat)
  
  top_annot = NULL
  if(transpose){
    if(!is.null(row_annot)){
      if('logical' %in% class(row_annot)){
        if(row_annot){
          library(jj)
          cols_use = jj_get_jj_colours(gtools::mixedsort(rownames(heatmap_mat)))
          cols_use = cols_use[names(cols_use) %in% rownames(heatmap_mat)]
          #print(cols_use)
          top_annot = HeatmapAnnotation(cluster = rownames(heatmap_mat), col=list(cluster=cols_use))
        }
      }else{
        top_annot = row_annot
      }
    }
  }else{
    if(!is.null(row_annot)){
      if('logical' %in% class(row_annot)){
        if(row_annot){
          library(jj)
          cols_use = jj_get_jj_colours(gtools::mixedsort(rownames(heatmap_mat)))
          cols_use = cols_use[names(cols_use) %in% rownames(heatmap_mat)]
          #print(cols_use)
          row_annot = rowAnnotation(cluster = rownames(heatmap_mat), col=list(cluster=cols_use))
        }else{
          row_annot = NULL
        }
      }
    }
    
  }

  if(transpose){
    if(!is.null(names(features_use)) & show_top_annot){
      row_annot =  rowAnnotation(Feature=names(features_use), col = list(Feature=jj_get_jj_colours(names(features_use))))
    }else{
      row_annot = NULL
    }
  }else{
    if(!is.null(names(features_use)) & show_top_annot){
      top_annot =  HeatmapAnnotation(Feature=names(features_use), col = list(Feature=jj_get_jj_colours(names(features_use))))
    }else{
      top_annot = NULL
    }
  }
  
  name_use = 'mean expression'
  if(scale_data){
    heatmap_mat = scale(heatmap_mat)
    name_use = 'scaled mean expression'
  } 
  if(return_matrix) return(heatmap_mat)
  if(transpose){
    heatmap_mat = base::t(heatmap_mat)
  }
  h1 = Heatmap(heatmap_mat, name = name_use, top_annotation = top_annot, right_annotation = row_annot,
               cluster_rows = cluster_rows, cluster_columns = cluster_columns, ...)
  return(h1)
}

#' @rdname feature_heatmap
#' @export
jj_plot_dotplot = function(obj, features_use, group_vec, scale_data = FALSE,
                        cluster_rows = TRUE,  cluster_columns = TRUE,
                        cap_top = NULL, cap_bottom = NULL, transpose = F, max_size = 6){
  
  get_fraction_expressing = function (rna_mat, features_use, group_vec){
    rna_mat = rna_mat[rownames(rna_mat) %in% features_use, ,
                      drop = FALSE] %>% as.matrix %>% base::t()
    cell_exp_ct_mat <- rna_mat %>% as.data.frame %>%
      dplyr::mutate(cluster = group_vec) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise_all(list(function(x) sum(x >  0))) %>%
      as.data.frame %>% column_to_rownames("cluster") %>%
      base::t(.) %>% as.data.frame %>% rownames_to_column("Gene") %>%
      tidyr::pivot_longer(!Gene, names_to = "cluster", values_to =
                            "cell_exp_ct")
    cell_ct_mat <- data.frame(cluster = group_vec) %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarise(cell_ct = n()) %>% as.data.frame
    cell_exp_ct_mat = cell_exp_ct_mat %>% dplyr::left_join(cell_ct_mat,
                                                           by = "cluster")
    return(cell_exp_ct_mat)
  }
  
  get_mean_expressing = function (rna_mat, features_use, group_vec,
                                  scale_data = FALSE){
    rna_mat = rna_mat[rownames(rna_mat) %in% features_use, ,
                      drop = FALSE] %>% as.matrix %>% base::t()
    
    expr_mat <- rna_mat %>% as.data.frame %>% dplyr::mutate(cluster =
                                                              group_vec) %>%
      dplyr::group_by(cluster) %>% dplyr::summarise_all(list(mean)) %>%
      as.data.frame %>% column_to_rownames("cluster")
    
    if (scale_data) {
      expr_mat = scale(expr_mat)
    }
    
    expr_mat = base::t(expr_mat) %>%
      as.data.frame %>%
      rownames_to_column("Gene") %>%
      tidyr::pivot_longer(!Gene, names_to = "cluster", values_to = "count")
    return(expr_mat)
  }
  
  if('Seurat' %in% class(obj)){
    message('Class of obj is Seurat, using the DefaultAssay: ', DefaultAssay(obj))
    obj = GetAssayData(obj)
  }else{
    message('Class of obj is not Seurat, assuming it is a matrix of normalized counts (features in rows, cells in columns)')
  }
  
  stopifnot((is(obj, "Matrix") | class(obj) == "matrix"))
  stopifnot(identical(ncol(obj), length(group_vec)))
  stopifnot(all(features_use %in% rownames(obj)))
  stopifnot(length(features_use) > 1)
  
  features_use = unique(features_use)
  cell_exp_ct_mat = get_fraction_expressing(obj, features_use,
                                            group_vec)
  expr_mat = get_mean_expressing(obj, features_use, group_vec,
                                 scale_data = scale_data)
  expr_mat = expr_mat %>% dplyr::left_join(cell_exp_ct_mat,
                                           by = c("Gene", "cluster"))
  
  markers <- expr_mat$Gene %>% unique()
  mat <- expr_mat %>% dplyr::filter(Gene %in% markers) %>%
    dplyr::select(Gene, cluster, count) %>%
    pivot_wider(names_from = cluster, values_from = count) %>%
    as.data.frame()
  row.names(mat) <- mat$Gene
  mat <- mat[, -1]
  
  if(cluster_rows) {
    clust <- hclust(dist(mat %>% as.matrix()))
  }else {
    clust = data.frame(labels = features_use, order =
                         1:length(features_use))
  }
  
  if(cluster_columns){
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
  }else{
    v_clust = data.frame(labels = levels(as.factor(group_vec)),
                         order = 1:length(unique(as.character(group_vec))))
  }
  
  expr_mat_use = expr_mat %>% dplyr::filter(Gene %in% markers) %>%
    dplyr::mutate(Pct_Expressing = (cell_exp_ct/cell_ct) *
                    100, Gene = factor(Gene, levels =
                                         clust$labels[clust$order]),
                  cluster = factor(cluster, levels =
                                     v_clust$labels[v_clust$order]))
  
  if(!scale_data) {
    expr_mat_use = expr_mat_use %>% dplyr::filter(count > 0, Pct_Expressing > 1)
  }else{
    expr_mat_use = expr_mat_use %>% dplyr::filter(Pct_Expressing >  1)
  }
  expr_mat_use$count = jj::jj_cap_vals(expr_mat_use$count, cap_top =
                                         cap_top, cap_bottom = cap_bottom)
  
  #use theme_cowplot from the cowplot package
  theme_cowplot = function(font_size = 14, font_family = "", line_size = 0.5, 
                           rel_small = 12/14, rel_tiny = 11/14, rel_large = 16/14){
    half_line <- font_size/2
    small_size <- rel_small * font_size
    theme_grey(base_size = font_size, base_family = font_family) %+replace% 
      theme(line = element_line(color = "black", size = line_size, 
                                linetype = 1, lineend = "butt"),
            rect = element_rect(fill = NA,  color = NA, size = line_size, linetype = 1), 
            text = element_text(family = font_family, face = "plain", color = "black",
                                size = font_size, hjust = 0.5, vjust = 0.5,
                                angle = 0, lineheight = 0.9, margin = margin(), debug = FALSE), 
            axis.line = element_line(color = "black", size = line_size, lineend = "square"),
            axis.line.x = NULL, axis.line.y = NULL,
            axis.text = element_text(color = "black", size = small_size), 
            axis.text.x = element_text(margin = margin(t = small_size/4), vjust = 0.5, angle = 90, hjust = 1),
            axis.text.x.top = element_text(margin = margin(b = small_size/4), vjust = 0), 
            axis.text.y = element_text(margin = margin(r = small_size/4), hjust = 1),
            axis.text.y.right = element_text(margin = margin(l = small_size/4), hjust = 0),
            axis.ticks = element_line(color = "black", size = line_size),
            axis.ticks.length = unit(half_line/2,"pt"),
            axis.title.x = element_text(margin = margin(t = half_line/2), vjust = 1),
            axis.title.x.top = element_text(margin = margin(b = half_line/2), vjust = 0), 
            axis.title.y = element_text(angle = 90, margin = margin(r = half_line/2), vjust = 1), 
            axis.title.y.right = element_text(angle = -90, margin = margin(l = half_line/2), vjust = 0), 
            legend.background = element_blank(), 
            legend.spacing = unit(font_size, "pt"), 
            legend.spacing.x = NULL, 
            legend.spacing.y = NULL, legend.margin = margin(0, 0, 0, 0),
            legend.key = element_blank(), legend.key.size = unit(1.1 * font_size, "pt"), 
            legend.key.height = NULL, legend.key.width = NULL, 
            legend.text = element_text(size = rel(rel_small)), 
            legend.text.align = NULL, legend.title = element_text(hjust = 0), 
            legend.title.align = NULL, legend.position = "right", 
            legend.direction = NULL, legend.justification = c("left", "center"),
            legend.box = NULL, legend.box.margin = margin(0,  0, 0, 0), 
            legend.box.background = element_blank(), 
            legend.box.spacing = unit(font_size, "pt"), panel.background = element_blank(), 
            panel.border = element_blank(), panel.grid = element_blank(), 
            panel.grid.major = NULL, panel.grid.minor = NULL, 
            panel.grid.major.x = NULL, panel.grid.major.y = NULL, 
            panel.grid.minor.x = NULL, panel.grid.minor.y = NULL, 
            panel.spacing = unit(half_line, "pt"), panel.spacing.x = NULL, 
            panel.spacing.y = NULL, panel.ontop = FALSE, strip.background = element_rect(fill = "grey80"), 
            strip.text = element_text(size = rel(rel_small), 
                                      margin = margin(half_line/2, half_line/2, half_line/2, 
                                                      half_line/2)), strip.text.x = NULL, 
            strip.text.y = element_text(angle = -90), 
            strip.placement = "inside", strip.placement.x = NULL, 
            strip.placement.y = NULL, strip.switch.pad.grid = unit(half_line/2, "pt"), 
            strip.switch.pad.wrap = unit(half_line/2, "pt"),
            plot.background = element_blank(),
            plot.title = element_text(face = "bold", size = rel(rel_large), hjust = 0, vjust = 1, 
                                      margin = margin(b = half_line)), 
            plot.subtitle = element_text(size = rel(rel_small), hjust = 0, vjust = 1, margin = margin(b = half_line)), 
            plot.caption = element_text(size = rel(rel_tiny), 
                                        hjust = 1, vjust = 1, margin = margin(t = half_line)), 
            plot.tag = element_text(face = "bold", hjust = 0, 
                                    vjust = 0.7), plot.tag.position = c(0, 1),
            plot.margin = margin(half_line, half_line, half_line, half_line), complete = TRUE)
  }
  
  if(transpose){
    x = 'Gene'
    y = 'cluster'
  }else{
    x = 'cluster'
    y = 'Gene'
  }
  
  dotplot <- ggplot(expr_mat_use, aes_string(x = x, y = y,
                                      color = 'count', size = 'Pct_Expressing')) + geom_point() +
    theme_cowplot() + 
    # theme_minimal() +
    # theme(axis.line = element_line(color = "black", size = 0.5, lineend = "square"),
    #       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    #       axis.ticks = element_blank(),
    #       #panel.border = element_blank(),
    #       panel.grid = element_blank(),
    #       axis.ticks = element_line(color = "black", size = 0,5)) +
    viridis::scale_color_viridis() + scale_y_discrete(position =
                                                        "right") +
    scale_radius(range = c(0,max_size))
  if(scale_data){
    dotplot = dotplot + labs(colour='Scaled\nmean count', x='', y='', size = '% detected')
  }else{
    dotplot = dotplot + labs(colour='Mean count', x='', y='', size = '% detected')
  }
  
  dotplot
}

# jj_plot_dotplot = function(obj, features_use, group_vec, group_annot_df = NULL, scale_data=FALSE, 
#                            column_colors=NULL,  cluster_rows=TRUE, cluster_columns=TRUE){
#   library(ggdendro) #install.packages('ggdendro')
#   library(cowplot)
#   #devtools::install_github("YuLab-SMU/ggtree")
#   #install.packages('aplot') #xlim2 ylim2
#   library(aplot)
#   library(ggtree)
#   library(patchwork) 
#   library(tidyverse)
#   #rna mat n genes (rows) * m samples (columns)
#   #group_annot_vec named vector with levels in goup_vec as values and annotations as names
#   stopifnot((is(obj, 'Matrix') | class(obj) == 'matrix'))
#   stopifnot(identical(ncol(obj), length(group_vec)))
#   stopifnot(all(features_use %in% rownames(obj)))
#   if(!is.null(column_colors)){
#     stopifnot(is.list(column_colors))
#   }
#   
#   #do not use scaled matrix here!
#   #get number of cells per cluster and number of cells with expression > 0
#   cell_exp_ct_mat = get_fraction_expressing(obj, features_use, group_vec)
#   
#   #get mean expression per cluster (can use normalized or scaled data here)
#   expr_mat = get_mean_expressing(obj, features_use, group_vec, scale_data = scale_data)
#   
#   expr_mat = expr_mat %>%  
#     dplyr::left_join(cell_exp_ct_mat, by=c('Gene', 'cluster'))
#   
#   if(!is.null(group_annot_df)){
#     #stopifnot(colnames(group_annot_df) %in% c('cluster','Group'))
#     stopifnot('cluster' %in% colnames(group_annot_df))
#     annot_cols = colnames(group_annot_df)[!colnames(group_annot_df) == 'cluster']
#     group_annot_df = group_annot_df[!duplicated(group_annot_df), ]
#     expr_mat = expr_mat %>%
#       dplyr::left_join(group_annot_df, by = 'cluster') 
#   }
#   
#   if(!is.null(names(features_use))){
#     gene_annot_df = data.frame(Gene=features_use, gene_annot = names(features_use))
#     expr_mat = expr_mat %>% dplyr::left_join(gene_annot_df, by = 'Gene')
#   }
#   
#   #should look like
#   # Gene  cluster  count cell_ct cell_exp_ct Group gene_annot
#   # <chr> <chr>    <dbl>   <int>       <int> <chr> <chr>     
#   #   1 Areg  0       0.0301    4579         180 A     tisTreg   
#   #stopifnot(colnames(expr_mat) == c('Gene',  'cluster',  'count', 'cell_ct', 'cell_exp_ct'))
#   
#   markers <- expr_mat$Gene %>% unique()
#   # make data square to calculate euclidean distance
#   mat <- expr_mat %>% 
#     dplyr::filter(Gene %in% markers) %>% 
#     dplyr::select(Gene, cluster, count) %>%  # drop unused columns to faciliate widening
#     pivot_wider(names_from = cluster, values_from = count) %>% 
#     as.data.frame() # make df as tibbles -> matrix annoying
#   row.names(mat) <- mat$Gene  # put gene in `row`
#   mat <- mat[,-1] #drop gene column as now in rows
#   if(cluster_rows){
#     clust <- hclust(dist(mat %>% as.matrix())) # hclust with distance matrix
#     ddgram <- as.dendrogram(clust) # create dendrogram
#     ggtree_plot <- ggtree::ggtree(ddgram)
#   }else{
#     clust= data.frame(labels = features_use, order = 1:length(features_use))
#     ggtree_plot = plot_spacer()
#   }
#   
#   if(cluster_columns){
#     v_clust <- hclust(dist(mat %>% as.matrix() %>% t())) # hclust with distance matrix
#     ddgram_col <- as.dendrogram(v_clust)
#     ggtree_plot_col <- ggtree(ddgram_col) + layout_dendrogram()
#   }else{
#     vclust= data.frame(labels = levels(as.factor(group_vec)), order = 1:length(unique(as.character(group_vec))))
#     ggtree_plot_col = plot_spacer()
#   }
#   
#   
#   expr_mat_use = expr_mat %>% dplyr::filter(Gene %in% markers) %>% 
#     dplyr::mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100,
#                   Gene = factor(Gene, levels = clust$labels[clust$order]),
#                   cluster = factor(cluster, levels = v_clust$labels[v_clust$order]))
#   if(!scale_data){
#     expr_mat_use = expr_mat_use %>% 
#       dplyr::filter(count > 0, `% Expressing` > 1)
#   }else{
#     expr_mat_use = expr_mat_use %>% 
#       dplyr::filter(`% Expressing` > 1)
#   }
#   
#   dotplot <- ggplot(expr_mat_use, aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
#     geom_point() + 
#     cowplot::theme_cowplot() + 
#     theme(axis.line  = element_blank()) +
#     theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#     ylab('') +
#     theme(axis.ticks = element_blank()) + viridis::scale_color_viridis() + 
#     #scale_color_gradientn(colours = viridis::viridis(20), limits = c(0,4), oob = scales::squish, name = 'log2 (count + 1)') +
#     scale_y_discrete(position = "right")
#   #################################################
#   ggtree_plot_col <- suppressMessages(ggtree_plot_col + xlim2(dotplot))
#   ggtree_plot <- suppressMessages(ggtree_plot + ylim2(dotplot))
#   
#   labels = legend_cols = list()
#   if(!is.null(group_annot_df)){
#     for(i in seq_along(annot_cols)){
#       labels[[i]] <- ggplot(expr_mat %>% 
#                               #dplyr::rename()
#                               mutate(#Group = Group,
#                                 cluster = factor(cluster, levels = v_clust$labels[v_clust$order])), 
#                             aes_string(x = 'cluster', y = 1, fill = annot_cols[i])) + 
#         geom_tile() + 
#         theme_nothing() +
#         xlim2(dotplot)
#       #annotate("text", x = 0, y = 1, label = annot_cols[i])  
#       if(!is.null(column_colors)){
#         labels[[i]] = labels[[i]] + scale_fill_manual(values=column_colors[[annot_cols[i]]]) 
#       }else{
#         labels[[i]] = labels[[i]] + scale_fill_manual(values=msPickSampleColors(as.data.frame(expr_mat)[, annot_cols[i]]))
#       }
#       
#       legend_cols[[i]] <- plot_grid(get_legend(labels[[i]] + theme(legend.position="bottom")))
#     }
#   }else{
#     labels[[1]] = legend_cols[[1]] = plot_spacer()
#   }
#   
#   if(!is.null(expr_mat$gene_annot)){
#     labels2 <- ggplot(expr_mat %>% 
#                         mutate(Gene = factor(Gene, levels = clust$labels[clust$order])), 
#                       aes(x = Gene, y = 1, fill = gene_annot)) + 
#       geom_tile() + 
#       scale_fill_brewer(palette = 'Set3') + 
#       theme_nothing() +
#       ylim2(dotplot) + coord_flip()
#     legend2 <- plot_grid(get_legend(labels2 + theme(legend.position="bottom")))
#   }else{
#     labels2 = legend2 = plot_spacer()
#   }
#   
#   gg = plot_spacer() + plot_spacer() + plot_spacer() + ggtree_plot_col
#   
#   for(i in seq_along(labels)){
#     gg = gg +  plot_spacer() + 
#       plot_spacer() + plot_spacer() + labels[[i]] 
#   }
#   
#   gg = gg +
#     plot_spacer() + plot_spacer() +  plot_spacer()+ plot_spacer() +
#     ggtree_plot + labels2 + plot_spacer()+ dotplot 
#   
#   for(i in seq_along(legend_cols)){
#     gg = gg + plot_spacer() + plot_spacer() + plot_spacer() + legend_cols[[i]]
#   }
#   
#   gg = gg  + plot_spacer() + plot_spacer() + plot_spacer()+ legend2
#   
#   gg = gg + plot_layout(ncol = 4, widths = c(0.7, 0.4, -0.6, 4), heights = c(0.9, rep(0.2, length(labels)), -0.5, 4, rep(0.7, length(legend_cols)), 1))
#   
#   return(gg)
# }
# 
# 
# get_fraction_expressing = function(rna_mat, features_use, group_vec){
#   #rna_mat genes (rows) * cells (columns)
#   rna_mat = rna_mat[rownames(rna_mat) %in% features_use,,drop=FALSE ]  %>% base::t()
#   cell_exp_ct_mat <- rna_mat %>%
#     as.data.frame %>% 
#     dplyr::mutate(cluster=group_vec) %>% 
#     dplyr::group_by(cluster) %>%
#     dplyr::summarise_all(list(function(x) sum(x>0))) %>% 
#     as.data.frame %>%
#     column_to_rownames('cluster') %>%
#     base::t %>% 
#     as.data.frame %>%
#     rownames_to_column('Gene') %>% 
#     tidyr::pivot_longer(!Gene, names_to = 'cluster', values_to = 'cell_exp_ct')
#   cell_ct_mat <- data.frame(cluster=group_vec) %>% 
#     dplyr::group_by(cluster) %>% 
#     dplyr::summarise(cell_ct=n()) %>% 
#     as.data.frame 
#   cell_exp_ct_mat = cell_exp_ct_mat %>% 
#     dplyr::left_join(cell_ct_mat, by = 'cluster')
#   return(cell_exp_ct_mat)
# }
# 
# get_mean_expressing =  function(rna_mat, features_use, group_vec, scale_data=FALSE){
#   #rna_mat genes (rows) * cells (columns)
#   rna_mat = rna_mat[rownames(rna_mat) %in% features_use,,drop=FALSE ]  %>% base::t()
#   
#   if(scale_data){
#     rna_mat = scale(rna_mat)
#   }
#   expr_mat <- rna_mat %>%
#     as.data.frame %>% 
#     dplyr::mutate(cluster=group_vec) %>% 
#     dplyr::group_by(cluster) %>% 
#     dplyr::summarise_all(list(mean)) %>% 
#     as.data.frame %>% column_to_rownames('cluster') %>%
#     base::t %>% 
#     as.data.frame %>%
#     rownames_to_column('Gene') %>% 
#     tidyr::pivot_longer(!Gene, names_to = 'cluster', values_to = 'count') 
#   return(expr_mat)
# }
# 
# get_sample_points = function(df, feature_column, group_column, sample_size=500){
#   df <- df %>%
#     dplyr::group_by(!!rlang::sym(group_column)) %>%
#     dplyr::add_count(name = "n") %>%
#     dplyr::mutate(n = ifelse(n > sample_size, sample_size, n)) %>%
#     dplyr::sample_n(n[1]) %>%
#     dplyr::select(-n) 
#   return(df)
# }
# 
# 
# get_mean_pct_per_group = function(rna_mat, features_use, group_vec, return_summarized_genes=TRUE){
#   stopifnot(identical(ncol(rna_mat), length(group_vec)))
#   stopifnot(all(features_use %in% rownames(rna_mat)))
#   
#   rna_mat = rna_mat[rownames(rna_mat) %in% features_use, ]  %>% base::t()
#   
#   expr_mat <- rna_mat %>%
#     as.data.frame %>% dplyr::mutate(cluster=group_vec) %>% 
#     dplyr::group_by(cluster) %>% dplyr::summarise_all(list(mean)) %>% 
#     as.data.frame %>% column_to_rownames('cluster') %>% base::t %>%  as.data.frame %>%
#     rownames_to_column('Gene') %>% 
#     tidyr::pivot_longer(!Gene, names_to = 'cluster', values_to = 'count')
#   
#   scaled_expr_mat<- scale(rna_mat) %>%
#     as.data.frame %>% dplyr::mutate(cluster=group_vec) %>% 
#     dplyr::group_by(cluster) %>% dplyr::summarise_all(list(mean)) %>% 
#     as.data.frame %>% column_to_rownames('cluster') %>% base::t %>%  as.data.frame %>%
#     rownames_to_column('Gene') %>% 
#     tidyr::pivot_longer(!Gene, names_to = 'cluster', values_to = 'scaled_count')
#   expr_mat = expr_mat %>% dplyr::left_join(scaled_expr_mat, by = c('Gene', 'cluster'))
#   
#   cell_exp_ct_mat <- rna_mat %>%
#     as.data.frame %>% dplyr::mutate(cluster=group_vec) %>% 
#     dplyr::group_by(cluster) %>% dplyr::summarise_all(list(function(x) sum(x>0))) %>% 
#     as.data.frame %>% column_to_rownames('cluster') %>% base::t %>%  as.data.frame %>%
#     rownames_to_column('Gene') %>% 
#     tidyr::pivot_longer(!Gene, names_to = 'cluster', values_to = 'cell_exp_ct')
#   
#   cell_ct_mat <- data.frame(cluster=group_vec) %>% 
#     dplyr::group_by(cluster) %>% dplyr::summarise(cell_ct=n()) %>% 
#     as.data.frame 
#   
#   expr_mat = expr_mat %>%
#     dplyr::left_join(cell_ct_mat, by = 'cluster') %>% 
#     dplyr::left_join(cell_exp_ct_mat, by=c('Gene', 'cluster')) %>% 
#     dplyr::mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100)
#   
#   if(return_summarized_genes){
#     summarized_expr_mat = expr_mat %>% 
#       dplyr::group_by(cluster) %>% 
#       dplyr::summarise(mean_pct_expressing = mean(`% Expressing`),
#                        mean_expr_count = mean(count),
#                        mean_scaled_expr_count = mean(scaled_count),
#                        cluster_cell_count = unique(cell_ct))
#     return(summarized_expr_mat)
#   }else{
#     return(expr_mat)
#   }
# }
