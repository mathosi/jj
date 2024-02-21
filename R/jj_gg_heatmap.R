#' feature heatmap based on ggplot
#'
#' @name jj_gg_heatmap
#' @param df data.frame in long format
#' @param x_group Group plotted in rows
#' @param y_group Group plotted in columns
#' @param value Continuous feature used to colour the tiles
#' @param y_annot Optional group used to annotate the columns (be careful not to mix up with rows)
#' @param x_annot Optional group used to annotate the rows (be careful not to mix up with columns)
#' @param cluster_rows logical, should rows be clustered
#' @param cluster_columns logical, should rows be clustered
#' @param order_by_annot instead of clustering columns by values, cluster columns by mixedorder of y annotation
#' @param show_numbers Plot the values as number in each tile
#' @param col_scale Colorscale to use, one of rwb0 (red <0, white = 0, blue >0), viridis
#' @param size_x size of x axis labels
#' @param size_y size of y axis labels
#' @export
#' @examples
#' data('mtcars')
#' mtcars2 = mtcars %>% dplyr::group_by(cyl, gear) %>% dplyr::summarise(mpg = mean(mpg))
#' mtcars2$cyl_annot = from_to(mtcars2$cyl, c('4' = 'few', '6' = 'many', '8' = 'many'))
#' jj_gg_heatmap(mtcars2, x_group = 'cyl', y_group = 'gear', value='mpg', show_numbers = T, col_scale = 'viridis', x_annot = 'cyl_annot')

jj_gg_heatmap = function(df, x_group, y_group, value, y_annot = NULL, x_annot = NULL, order_by_annot = FALSE, cluster_columns = TRUE, cluster_rows = T, show_numbers = FALSE, size_x = 8, size_y = 10, col_scale = 'rwb0'){
  library(cowplot)
  library(tibble)
  library(tidyr)
  library(dplyr)
  stopifnot(all(c(x_group, y_group, value) %in% colnames(df)))
  if(!is.null(y_annot)) stopifnot(y_annot %in% colnames(df))
  if(!is.null(x_annot)) stopifnot(x_annot %in% colnames(df))
  df = as.data.frame(df)
  stopifnot(!anyDuplicated(df[, c(x_group, y_group, value)]))
  df$AnnotationString = 'dummy'
  
  if(is.character(col_scale)){
    col_scale = match.arg(col_scale, choices = c('rwb0', 'viridis'))
    if(col_scale == 'rwb0'){
      col_scale = scale_fill_gradient2(low = "indianred",  mid = "whitesmoke", high = "darkblue", midpoint = 0)
    }else if(col_scale == 'viridis'){
      col_scale = viridis::scale_fill_viridis()
    }
  }
  corder = unique(df[, x_group])
  rorder = unique(df[, y_group])
  if(cluster_columns | cluster_rows){
    df_wide = df %>% pivot_wider(id_cols = x_group, names_from = y_group, values_from = value) %>% column_to_rownames(x_group)
    if(anyNA(df_wide)){
      warning('Some values are NA. Using imputated values (for clustering only)')
      df_wide = impute::impute.knn(as.matrix(df_wide))$data
    }
    if(cluster_columns){
      corder = rownames(df_wide)[hclust(dist(df_wide))$order]
      df[, x_group] = factor(df[, x_group], levels = corder)
    }
    if(cluster_rows){
      rorder = colnames(df_wide)[hclust(dist(t(df_wide)))$order]
      df[, y_group] = factor(df[, y_group], levels = rorder)
    }
  }
  
  if(!is.null(y_annot)){
    if(length(unique(df[, y_annot])) > length(rorder)) message('Did you accidentally set y_annot instead of x_annot?')
    if(order_by_annot){
      if(cluster_rows) warning('Both cluster_rows and order_by_annot are set to TRUE. Ordering by annotation.')
      rorder = df[, c(y_group, y_annot)] %>% #dplyr::mutate(y_annot = mixsort_factor(y_annot)) %>%  
        dplyr::arrange_(y_annot) %>% pull(y_group) %>% unique
    }
    df[, y_group] = factor(df[, y_group], levels = rorder)
    annot_y_gg = ggplot(data = df, mapping = aes_string(x='AnnotationString', y = y_group , fill = y_annot)) + geom_tile() + theme_void()
  }
  
  if(!is.null(x_annot)){
    if(length(unique(df[, x_annot])) > length(corder)) message('Did you accidentally set x_annot instead of y_annot?')
    if(order_by_annot){
      if(cluster_columns) warning('Both cluster_rows and order_by_annot are set to TRUE. Ordering by annotation.')
      corder = df[, c(x_group, x_annot)] %>% #dplyr::mutate(y_annot = mixsort_factor(y_annot)) %>%  
        dplyr::arrange_(x_annot) %>% pull(x_group) %>% unique
    }
    df[, x_group] = factor(df[, x_group], levels = corder)
    annot_x_gg = ggplot(data = df, mapping = aes_string(y='AnnotationString', x = x_group , fill = x_annot)) + geom_tile() + theme_void()
  }
  
  gg = ggplot(data = df, mapping = aes_string(x=x_group, y = y_group, fill = value)) + geom_tile() + 
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1, size = size_x), axis.text.y = element_text(size = size_y)) + 
    col_scale
  
  if(show_numbers){
    gg = gg + geom_text(aes_string(label = value, colour=value)) + theme(legend.title = element_text(size = 12)) +  
      binned_scale(aesthetics = "color", scale_name = "stepsn", guide = 'none', palette = function(x) c("white", "black"),
                   breaks =  max(df[, value])/3, limits = c(min(df[, value]),   max(df[, value]) ))
  }
  
  if(!is.null(y_annot) & !is.null(x_annot)){
    lgd = plot_grid(get_legend(gg), get_legend(annot_x_gg), get_legend(annot_y_gg), align = 'v', ncol = 1)
    annot_x_gg = annot_x_gg + theme(legend.position =  'none')
    annot_y_gg = annot_y_gg + theme(legend.position =  'none')
    gg = gg + theme(legend.position =  'none')
    gg2 = get_panel(gg)
    gg = plot_grid(annot_x_gg, NULL, gg, annot_y_gg,  nrow = 2,ncol = 2, align =  'hv', rel_heights = c(0.1, 0.8), rel_widths = c(0.8, 0.1))
    #annot_y_gg = plot_grid(annot_y_gg, gg2, align = 'h', ncol = 2, rel_widths = c(1,0))
    #annot_y_gg = plot_grid(NULL, annot_y_gg, align = 'v', nrow = 2, rel_heights = c(0.05, 0.8)) 
    plot_grid(gg, NULL, lgd, align = 'h', ncol = 3, rel_widths = c(0.8, 0.05, 0.1)) 
    #plot_grid(gg, annot_y_gg, NULL, lgd, align = 'h', ncol = 4, rel_widths = c(0.8, 0.05, 0.05, 0.1)) 
  }else if(!is.null(y_annot)){
    lgd = plot_grid(get_legend(gg), get_legend(annot_y_gg), align = 'v', ncol = 1)
    annot_y_gg = annot_y_gg + theme(legend.position =  'none')
    gg = gg + theme(legend.position =  'none')
    plot_grid(gg, annot_y_gg, NULL, lgd, align = 'h', nrow = 1, rel_widths = c(0.8, 0.05, 0.05, 0.1)) 
  }else if(!is.null(x_annot)){
    lgd = plot_grid(get_legend(gg), get_legend(annot_x_gg), align = 'h', nrow = 1)
    annot_x_gg = annot_x_gg + theme(legend.position =  'none')
    gg = gg + theme(legend.position =  'none')
    gg = plot_grid(annot_x_gg, gg, align = 'v', nrow = 2, rel_heights = c(0.05, 0.8)) 
    plot_grid(gg, NULL, lgd, align = 'h', nrow = 1, rel_widths = c(0.75, 0.05, 0.1)) 
  }else{
    return(gg)
  }
}
