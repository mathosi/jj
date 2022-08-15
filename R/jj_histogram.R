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
#' jj_plot_hist_wo0(feature_vec = pmax(rnorm(1000), 0), title = 'myfeat', thres = 1)

jj_plot_hist_wo0 = function(feature_vec, thres=NULL, title = '', xlab='value'){
  goi_df = data.frame(feature=feature_vec)
  if(is.numeric(thres)){
    goi_df[, 'above_threshold'] = goi_df[, 'feature'] > thres
    goi_gr_thres = sum(goi_df[, 'feature'] > thres)
    goi_smeq_thres = sum(goi_df[, 'feature'] <= thres)
    goi_smeq_thres_not_0 = sum(goi_df[, 'feature'] <= thres & goi_df[, 'feature'] != 0)
    ggplot(goi_df, aes_string(x='feature', fill= 'above_threshold')) +
      geom_histogram(bins = 300) + 
      scale_x_continuous(limits = c(0.01,NA)) + 
      labs(title = title, subtitle = paste0(sprintf('%-*s', 18, paste0('> ', thres,':')), goi_gr_thres, '\n',
                                          sprintf('%-*s', 17, paste0('<= ', thres,':')) , goi_smeq_thres, '\n',
                                          sprintf('%-*s', 14, paste0('<= ', thres, ' & > 0:')), goi_smeq_thres_not_0),
           x=xlab, y = 'n') + theme_minimal()
  }else{
    ggplot(goi_df, aes_string(x='feature')) +
      geom_histogram(bins = 300) + 
      scale_x_continuous(limits = c(0.01,NA)) + 
      labs(title = title, x=xlab, y = 'n') + 
      theme_minimal()
  }
}
