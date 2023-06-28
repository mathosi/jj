#' show results from differential expression analysis
#'
#' jj_plot_volcano volcano plot of pval (which is -log10 transformed) against logFC
#'
#' @name jj_plot_volcano
#' @param logfc_column sparse matrix, vector, or data.frame
#' @param pval_column vector with group annotation, has to have length equal to ncol(sparse_mat)/nrow(data.frame)/length(vector)
#' @param symbol_column provide `symbol_column` to add text to all points that are highlighted (when use_text=TRUE)
#' @param labs_range vector specifying range of x axis (and for fc_fc plot range of y axis)
#' @param marker_thres define `marker_thres` to highlight all points with absolute value above
#' @param col_vec define larger, smaller, intermediate and highlight colours
#' @param markers_highlight vector of genes from the symbol_column which should additionally be highlighted
#' @param markers_highlight_col named vector assigning a colour to each gene to be highlighed
#' @param only_highlight if TRUE, only show markers specified in `markers_highlight`
#' @param col_by_highlight if TRUE, only colour the markers specified in `markers_highlight`
#' @param use_text plot text if there are any points to highlight
#' @param add_thres_line if TRUE, plot dashed lines at the defined `marker_thres`
#' @param alpha alpha level for points (working?)
#' @returns ggplot scatterplot of log10(pval) vs FC or FC vs FC
#' @export
#' @examples
#' #Using Seurats FindMarkers output
#' library(Seurat)
#' marker_df = FindMarkers(pbmc_small, group.by = 'groups', ident.1 = 'g1', ident.2 = 'g2')
#' marker_df$symbol = rownames(marker_df)
#' #or using a random data.frame
#' marker_df = data.frame(avg_log2FC = rnorm(100, 0, 2), 
#' p_val = sample(c(runif(50, 0, 0.1), runif(50))), 
#' symbol = paste0(sample(LETTERS, 100, replace = T), sample(letters, 100, replace = T), 1:100))
#' jj_plot_volcano(marker_df, logfc_column = 'avg_log2FC', pval_column = 'p_val', marker_thres = Inf) + theme(legend.position = 'none')
#' #colour highlighting above a threshold and capping very high and low values, adding lines for the threshold
#' jj_plot_volcano(marker_df, logfc_column = 'avg_log2FC', pval_column = 'p_val',
#'  labs_range = c(-2,2), marker_thres = 1.5, add_thres_line=T)
#' #adding text to highlighted markers, possibility to show unbalanced x axis
#' jj_plot_volcano(marker_df, logfc_column = 'avg_log2FC', pval_column = 'p_val', symbol_column = 'symbol', 
#' use_text = T, labs_range = c(-2,4), marker_thres = 2)
#' #instead of highlighting everything above/below a threshold, specific markers can be highlighted
#' jj_plot_volcano(marker_df, logfc_column = 'avg_log2FC', pval_column = 'p_val', symbol_column = 'symbol',
#'                labs_range = c(-2.5, 2.5), marker_thres = 4,
#'                markers_highlight = marker_df$symbol[1:4], pt.size = 1.5,
#'                only_highlight = F, col_by_highlight = T)
#' #specific highlighting and general highlighting can be combined
#' jj_plot_volcano(marker_df, logfc_column = 'avg_log2FC', pval_column = 'p_val', symbol_column = 'symbol',
#'                labs_range = c(-2.5, 2.5), marker_thres = 2,
#'                markers_highlight = marker_df$symbol[1:4], pt.size = 1.5, highlight.pt.size = 3,
#'                only_highlight = F, col_by_highlight = F)
#' #only show highlighted markers
#' jj_plot_volcano(marker_df, logfc_column = 'avg_log2FC', pval_column = 'p_val', symbol_column = 'symbol',
#'                labs_range = c(-2.5, 2.5), marker_thres = 2,
#'                markers_highlight =  marker_df$symbol[1:4], pt.size = 1.5, 
#'                only_highlight = T, col_by_highlight = F)
#' #to change colour and size of highlighted markers use `highlight.pt.size` and `markers_highlight_cols`                
#' jj_plot_volcano(marker_df, logfc_column = 'avg_log2FC', pval_column = 'p_val', symbol_column = 'symbol',
#'                labs_range = c(-4, 4), marker_thres = Inf,
#'                markers_highlight =  marker_df$symbol[1:4], pt.size = 1.5, highlight.pt.size = 3,
#'                only_highlight = F, col_by_highlight =T,
#'                markers_highlight_col = structure(c('red','orange','red','blue'), names = marker_df$symbol[1:4])) + theme(legend.position='none')

#' @export 
jj_plot_volcano = function(plot_df, logfc_column, pval_column, symbol_column=NULL, 
                           marker_thres = c(0.5,Inf), labs_range = c(-1, 1), alpha = 1,                          
                           pt.size=2.5, highlight.pt.size = 2.5, 
                           col_vec = c("red", "blue", "black","green"),
                           markers_highlight = NULL, markers_highlight_col=NULL, only_highlight=FALSE, 
                           col_by_highlight=FALSE, use_text=TRUE, add_thres_line = FALSE, group_names=NULL,
                           arrow_pos=list(x=0.25, xend=1, y=40)){
  #TODO: sym() throws error if it is not used from the dplyr package: Error in local_error_context(dots = dots, .index = i, mask = mask) : 
  #  promise already under evaluation: recursive default argument reference or earlier problems?
  plot_df = as.data.frame(plot_df)
  stopifnot(is.data.frame(plot_df))
  if(!is.null(group_names)){
    stopifnot(length(group_names) == 2)
  }
  stopifnot(length(labs_range) == 2)
  stopifnot(all(is.numeric(labs_range)))
  plot_df = plot_df[, colnames(plot_df) %in% c(logfc_column, pval_column, symbol_column)]
  symbols_use = c('outside plotting area', 'inside plotting area')
  if(is.null(marker_thres)){
    marker_thres = c(Inf, Inf)
  }else if(length(marker_thres)==1){
    marker_thres = c(marker_thres, Inf)
  }else if(length(marker_thres) > 2){
    stop('marker_thres should specifiy thresholds on log2FC and adjusted p value (length 2)')
  }
  
  if(!is.null(labs_range)){
    #add a column that is used for different symbols (circle if inside plot, triangle if outside)
    #dplyr throws error 'attempt to use zero-length variable name' if there is a colname ''
    plot_df = dplyr::mutate(plot_df, localization = ifelse(!!sym(logfc_column) >= labs_range[2], symbols_use[1], 
                                          ifelse(!!sym(logfc_column) <= labs_range[1],
                                                 symbols_use[1], symbols_use[2])))
    plot_df$localization = factor(plot_df$localization, levels = c('inside plotting area', 'outside plotting area'))
    
    #cap the values outside of plot
    plot_df[, logfc_column] = ifelse(plot_df[, logfc_column] >= labs_range[2], labs_range[2],
                                     ifelse(plot_df[, logfc_column] <= labs_range[1],
                                            labs_range[1], plot_df[, logfc_column]))
  }else{
    plot_df$localization = ""
  }
  
  #give label for all genes > marker_thres
  mthres = as.character(marker_thres[1])
  small_val = paste0(">= ", mthres)
  middle_val = sprintf(">-%s & <%s", mthres, mthres)
  big_val = paste0("<= -", mthres)
  
  plot_df$logFC = ifelse(plot_df[, logfc_column] >= marker_thres[1], small_val, 
                         ifelse(plot_df[, logfc_column]<=-marker_thres[1] , big_val,
                                middle_val))
  
  plot_df[, pval_column] = -log10(plot_df[, pval_column])
  mthres2 = as.character(marker_thres[2])
  small_val2 = paste0(">= ", mthres2)
  big_val2 = paste0("<= -", mthres2)
  plot_df$pval_thres = ifelse(plot_df[, pval_column] >= marker_thres[2], small_val2, big_val2)
  
  if(!is.null(symbol_column)){
    plot_df$label_use <- ''
    plot_df$label_use[plot_df$logFC != middle_val] <- plot_df[, symbol_column][plot_df$logFC != middle_val]
    plot_df$label_use[plot_df$pval_thres != big_val2] <- plot_df[, symbol_column][plot_df$pval_thres != big_val2]
    
    if(!is.null(markers_highlight) & any(plot_df[, symbol_column] %in% markers_highlight)){
      if(only_highlight){
        plot_df = plot_df[plot_df[, symbol_column] %in% markers_highlight, ]
      }
      #markers_highlight = unique(plot_df[, symbol_column])
      plot_df$label_use[plot_df[, symbol_column] %in% markers_highlight] <- plot_df[, symbol_column][plot_df[, symbol_column] %in% markers_highlight]
      plot_df$logFC[plot_df[, symbol_column] %in% markers_highlight] = 'highlight'
      plot_df = plot_df[order(plot_df[, symbol_column] %in% markers_highlight, decreasing = F), ]
    }
  }
  
  cols_use = structure(col_vec, names = c(big_val, small_val, middle_val, "highlight"))
  cols_use = cols_use[names(cols_use) %in% plot_df$logFC]
  
  
  if(col_by_highlight){
    if(is.null(markers_highlight_col)){
      cols_use = jj_get_jj_colours(markers_highlight)
    }else{
      cols_use = markers_highlight_col
    }
    if(!is.null(labs_range)){
      g1 = ggplot() +
        geom_point(data=plot_df, 
                   aes_string(x=logfc_column, y=pval_column, shape='localization'),
                   alpha=alpha, size=pt.size) + 
        geom_point(data=plot_df[plot_df[, symbol_column] %in% markers_highlight, ], 
                   aes_string(x=logfc_column, y=pval_column, shape='localization',
                              colour='label_use'), 
                   size=highlight.pt.size) 
    }else{
      g1 =  ggplot() +
        geom_point(data=plot_df, 
                   aes_string(x=logfc_column, y=pval_column), 
                   alpha=alpha, size=pt.size) +
        geom_point(data=plot_df[plot_df[, symbol_column] %in% markers_highlight, ],
                   aes_string(x=logfc_column, y=pval_column, #label='label_use',
                              colour='label_use'), size=highlight.pt.size)
    }
  }else{
    if(!is.null(labs_range)){
      g1 = ggplot(plot_df, aes_string(x=logfc_column, y=pval_column)) +
        geom_point(aes(colour = logFC, shape=localization), size=pt.size, alpha=alpha) 
      # if(!is.null(markers_highlight[1])){
      #    g1 = g1 + geom_point(data=plot_df[plot_df[, symbol_column] %in% markers_highlight, ], 
      #               aes_string(x=logfc_column, y=pval_column, shape='localization'), 
      #               size=highlight.pt.size, colour = 'green') 
      #  }
    }else{
      g1 = ggplot(plot_df, aes_string(x=logfc_column, y=pval_column), alpha=alpha) +
        geom_point(aes(colour = logFC), size=pt.size)
      # if(!is.null(markers_highlight[1])){
      #   g1 = g1 + geom_point(data=plot_df[plot_df[, symbol_column] %in% markers_highlight, ],
      #              aes_string(x=logfc_column, y=pval_column), size=highlight.pt.size)
      # }
    }
  }
  
  if(!is.null(symbol_column) & use_text){
    g1 = g1 + 
      ggrepel::geom_text_repel(data=plot_df, 
                               aes_string(x=logfc_column,
                                          y=pval_column,
                                          label='label_use'), 
                               max.overlaps = 500) 
  }
  
  
  g1 = g1 +
    scale_colour_manual(values = cols_use) + 
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) + theme_minimal()
  
  if(!is.null(labs_range)){
    g1 = g1 + coord_cartesian(xlim=labs_range, expand = F, clip = 'off')
  }
  
  if(add_thres_line){
    if(!is.infinite(marker_thres[1])){
      g1 = g1 + 
        geom_vline(xintercept = -marker_thres[1], linetype="dashed", color ='black', size = 0.5) +
        geom_vline(xintercept = marker_thres[1], linetype="dashed", color ='black', size = 0.5)
    }
    if(! is.infinite(marker_thres[2])){
      g1 = g1 + 
        geom_hline(yintercept = marker_thres[2], linetype="dashed", color ='black', size = 0.5)
    }
  }
  
  label_df = data.frame(name = c('logfc_gr0', 'logfc_sm0'),
                        n = c(sum(plot_df[, logfc_column] > 0), sum(plot_df[, logfc_column] < 0)),
                        x_pos = labs_range[c(2,1)],
                        y_pos = c(1,1))
  # if(add_numbers){
  #   g1 = g1 + geom_text(data = label_df[1, ],hjust = 1, vjust=0, mapping = aes(x = x_pos, y =  y_pos, label = n)) +
  #     geom_text(data = label_df[2, ],hjust = 0, vjust=0, mapping = aes(x = x_pos, y =  y_pos, label = n))
  # }
  if(!is.null(group_names[1])){
    g1 = g1 + annotate('segment', x = arrow_pos$x, xend = arrow_pos$xend, y = arrow_pos$y, yend = arrow_pos$y,  arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
      annotate('segment', x = -arrow_pos$x, xend = -arrow_pos$xend, y = arrow_pos$y, yend = arrow_pos$y,  arrow = arrow(type = "closed", length = unit(0.02, "npc"))) +
      annotate('text', x = arrow_pos$xend+0.01, y = arrow_pos$y, label = sprintf('%s (%i genes)', group_names[1], label_df$n[1]), hjust = 0, size=3) + 
      annotate('text', x = -arrow_pos$xend-0.01, y = arrow_pos$y, label = sprintf('%s (%i genes)', group_names[2], label_df$n[2]), hjust = 1, size = 3)
  }
  
  g1 = g1 + labs(y = paste(pval_column, '(-log10)'))
  
  return(g1)
}
