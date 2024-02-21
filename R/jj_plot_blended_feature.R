#' Plot scatterplot of two features as combined colorscale
#'
#' jj_plot_blended_feature translates two continuous feature scales to a red and green colorscale and plots a blended colorscale a scatterplot. 
#' Useful e.g. to visualize co-expression patterns.
#'
#' @name jj_plot_blended_feature
#' @param df sparse matrix, vector, or data.frame
#' @param x vector with group annotation, has to have length equal to ncol(sparse_mat)/nrow(data.frame)/length(vector)
#' @param y provide `symbol_column` to add text to all points that are highlighted (when use_text=TRUE)
#' @param colfeat1 vector specifying range of x axis (and for fc_fc plot range of y axis)
#' @param colfeat2 define `marker_thres` to highlight all points with absolute value above
#' @param pt_size size of points
#' @returns ggplot scatterplot
#' @export
#' @examples
#' library(Seurat)
#' dr_df = jj_get_reduction_coords(pbmc_small, 'tsne')
#' dr_df = jj_bind_features_with_dr_df(pbmc_small, assay= 'RNA', slot = 'data', features = c('CD19', 'CD79A'), dr_df = dr_df, cap_top = 'q95')
#' jj_plot_blended_feature(dr_df, x = 'tSNE_1', y = 'tSNE_2', colfeat1 = 'CD19', colfeat2 = 'CD79A', pt_size = 3)

jj_plot_blended_feature = function(df, x, y, colfeat1, colfeat2, pt_size = 1){
  
  library(scales)
  library(ggnewscale)
  library(gridExtra)
  
  feat1 = pull(df, colfeat1)
  feat2 = pull(df, colfeat2)
  
  df$colour = rgb(
    rescale(feat1),
    rescale(feat2),
    0
  )
  
  gg = ggplot(df, aes_string(x = x, y = y, colour = 'colour')) +
    geom_point(size = pt_size) +
    scale_colour_identity() +
    new_scale_colour() +
    # shape = NA --> invisible layers
    geom_point(aes_string(colour = colfeat1), shape = NA) +
    scale_colour_gradient(low = "black", high = "red") +
    new_scale_colour() +
    geom_point(aes_string(colour = colfeat2), shape = NA) +
    scale_colour_gradient(low = "black", high = "green") + 
    theme(legend.position = 'none')
  
  
  #bicolor legend
  colsteps = seq(0, 1, 0.1)
  colmat = jj_initialize_df(length(colsteps),length(colsteps))
  for(i in seq_along(colsteps)){
    colmat[i, ] = rgb(colsteps[i],
                      seq(0, 1, 0.1),
                      0)
  }
  colnames(colmat) = colsteps
  colmat$x = factor(colsteps, levels = colsteps)
  colmatt = tidyr::pivot_longer(colmat, 1:(ncol(colmat)-1), names_to = 'y', values_to = 'colour')
  colmatt$y = factor(colmatt$y, levels = colsteps)
  
  
  new_vals = round(seq(min(feat1), max(feat1), length.out = length(colsteps)), 2)
  names(new_vals) = levels(colmatt$x)
  
  colmatt$x = from_to(vec= colmatt$x, new_vals)
  
  new_vals = round(seq(min(feat2), max(feat2), length.out = length(colsteps)),2)
  names(new_vals) = levels(colmatt$y)
  
  colmatt$y = from_to(vec= colmatt$y, new_vals)
  
  gg_leg = ggplot(colmatt, aes(x=x, y=y, fill=colour)) + geom_tile() + 
    scale_fill_manual(values = structure(colmatt$colour, names = colmatt$colour)) +
    theme(legend.position = 'none') + 
    labs(x = colfeat1, y = colfeat2)
  #cowplot::plot_grid(plotlist = list(gg, gg_leg), nrow = 1, rel_widths = c(0.8, 0.2))
  
  gg_grid = grid.arrange(
    grobs = list(gg, gg_leg),
    widths = c(5, 2),
    layout_matrix = rbind(c(1, 2),
                          c(1, NA)))
  
  gg_grid
}
