#' plot histogram of data
#'
#' plot histogram of data
#'
#' @name plot_hist
#' @param feature vector with values to plot
#' @param fill optional value threshold, for which the number of observations are counted above and below
#' @param colour optional title
#' @param group xlab text
#' @param title colours for bars above and below threshold
#' @param nbins NULL (do not plot numbers) or vector containing any combination of 1 (observations above threshold), 2 (observations below threshold) and 3 (observations below threshold that are not 0)
#' @param size border size for the histogram stacks
#' @param exclude_zeros if TRUE, omit all 0 values from the plot (useful e.g. to plot sparse gene expression)
#' @param thres Value threshold used to fill the histogram. Also used to count the number of observations above and below in `include_summary`
#' @param include_summary Above the plot, include the  number of observations within some ranges. Vector with the possibilties c('n','>=', '<', '>0<'), where n is the total number of observations, >= and < are those above and below the threshold and >0< are those greater than zero, but below the threshold
#' @export
#' @examples
#' df = data.frame(letters = sample(letters[1:4], 100, replace = T),
#'                 groups = paste0('G_', sample(1:2, 100, replace = T)),
#'                 treatment = factor(sample(c('no','yes'), 100, replace = T, prob = c(0.9,0.1)), levels=c('yes','no')),
#'                 values = sample(c(rep(0, 20), rnorm(80, 4, 2))))
#' jj_plot_hist(df$values)
#' jj_plot_hist(df, feature = 'values', group = 'groups')
#' jj_plot_hist(df, feature = 'values', group = 'groups', fill = 'letters', nbins = 30)
#' jj_plot_hist(df, feature = 'values', colour='treatment', fill = 'letters', title='my title', group = 'groups', size = 1.3, binwidth = 0.4) +
#'   scale_fill_manual(values = jj_get_jj_colours(df$letters)) + scale_colour_manual(values = c(yes='black', no = "white"))
#' jj_plot_hist(df, feature = 'values', thres = 4)
#' jj_plot_hist(df, feature = 'values', thres = 4, exclude_zeros = T, include_summary = c('<', '>='))

#' @export
jj_plot_hist = function(data, feature = NULL, fill = NULL, thres = NULL, colour= NULL, group = NULL,
                        title = '',  nbins = 50, size = 1, exclude_zeros = FALSE, include_summary = c('n','>=', '<', '>0<'), ...){
  #helper function to plot histograms providing either a data.frame as `data` or a vector of values
  if(is.vector(data)){
    data = data.frame(feature=data)
    feature = 'feature'
  }else{
    data = as.data.frame(data)
  }
  
  if(!is.null(thres)){
    stopifnot(is.numeric(thres))
    if(!is.null(fill)) warning('Both `fill` and `thres` are specified. Using fill by thres.')
    fill = 'above_threshold'
    data$above_threshold = data[, feature] > thres
    if(!is.null(include_summary[1])){
      goi_gr_thres = sum(data[, feature] >= thres)
      goi_sm_thres = sum(data[, feature] < thres)
      goi_sm_thres_not_0 = sum(data[, feature] < thres & data[, feature] != 0)
      strings_vec = vector()
      strings_vec['n'] = paste0(nrow(data), ' ( n )')
      strings_vec['>='] = paste0(goi_gr_thres, ' ( n >= ', as.character(thres), ' )')
      strings_vec['<'] = paste0(goi_sm_thres, ' ( n < ', as.character(thres), ' )')
      strings_vec['>0<'] = paste0(goi_sm_thres_not_0, ' ( n > 0 & < ', as.character(thres), ' )') 
      string_plot = paste0(paste0(strings_vec[include_summary], collapse='\n'), paste0(rep('\n', length(strings_vec) - length(include_summary)), collapse = ''), collapse = '')
    }
  }
  
  if(exclude_zeros){
    data = data[data[, feature] != 0,, drop = F]
  }
  
  gg = ggplot(data, aes_string(x=feature, fill= fill, colour=colour)) +
    geom_histogram(bins = nbins, size = size, ... ) + 
    labs(title = title) + 
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5))
  
  if(!is.null(thres)){
    gg = gg + labs(fill = paste0('Value >= ', thres))
    if(!is.null(include_summary)) gg = gg + labs(subtitle = string_plot) #+ theme(text=element_text(family = 'mono'))  + 
  }
  
  if(!is.null(group)){
    gg = gg + facet_wrap(as.formula(paste("~", group)))
  }
  gg
}
