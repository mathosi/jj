#' summarise values as mean, mode, or sum by a group vector 
#'
#' plot histogram of data excluding 0 (useful for showing gene expression in sparse single cell data)
#' use `thres` to define threshold for counting number of cells above and below this mark
#'
#' @name summarise_vals
#' @param feature_vec vector with values to plot
#' @param thres optional value threshold, for which the number of observations are counted above and below
#' @param title optional title
#' @param xlab xlab text
#' @export
#' @examples
#' cd8a_count = rnbinom(1e6, size = 1, prob = 0.7)
#' cell_nCount_RNA = rnbinom(1e6, size = 1000, prob=0.7)
#' cd8a_norm_count = cd8a_count / cell_nCount_RNA
#' jj_plot_hist_wo0(feature_vec = cd8a_norm_count, title = 'CD8A counts', thres = 0.015)
#' jj_plot_hist_wo0(feature_vec = cd8a_norm_count, title = 'CD8A counts', thres = 0.015, 
#'                 plot_numbers = c(3,2), cols = c('grey', 'black')) +
#'   theme(legend.position = 'none')

jj_plot_hist_wo0 = function(feature_vec, thres=NULL, title = '', xlab='value', 
                            cols = c('red', 'lightblue'), plot_numbers=1:3){
  if(!all(feature_vec >=0)){
    stop('Histogram currently only implemented for counts >= 0.')
  }
  if(!is.null(plot_numbers)[1]){
    stopifnot(all(plot_numbers %in% 1:3))
  }
  names(cols) = c('FALSE','TRUE')
  goi_df = data.frame(feature=feature_vec)
  if(is.numeric(thres)){
    goi_df[, 'above_threshold'] = goi_df[, 'feature'] > thres
    goi_gr_thres = sum(goi_df[, 'feature'] > thres)
    goi_smeq_thres = sum(goi_df[, 'feature'] <= thres)
    goi_smeq_thres_not_0 = sum(goi_df[, 'feature'] <= thres & goi_df[, 'feature'] != 0)
    
    strings_vec = vector()
    strings_vec[1] = sprintf('%-20s%d', paste0('>  ', as.character(thres),': '), goi_gr_thres)
    strings_vec[2] = sprintf('%-20s%d', paste0('<= ', as.character(thres),': ') , goi_smeq_thres)
    strings_vec[3] = sprintf('%-20s%d', paste0('<= ', as.character(thres), ' & > 0: '), goi_smeq_thres_not_0)
    string_plot = paste0(strings_vec[plot_numbers], collapse='\n')
    
    gg = ggplot(goi_df, aes_string(x='feature', fill= 'above_threshold')) +
      geom_histogram(bins = 300) + 
      scale_x_continuous(limits = c(0.01,NA)) + 
      labs(title = title, 
           x=xlab, y = 'n') +
      theme_minimal() + 
      scale_fill_manual(values = cols)
    if(!is.null(plot_numbers)){
      gg = gg + 
        labs(subtitle = string_plot) + 
        theme(text=element_text(family = 'mono')) 
    }
  }else{
    gg = ggplot(goi_df, aes_string(x='feature')) +
      geom_histogram(bins = 300) + 
      scale_x_continuous(limits = c(0.01,NA)) + 
      labs(title = title, x=xlab, y = 'n') + 
      theme_minimal()
  }
  gg
}
