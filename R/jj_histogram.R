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
#' @export
#' @examples
#' df = data.frame(letters = sample(letters[1:4], 100, replace = T),
#'                 groups = paste0('G_', sample(1:2, 100, replace = T)),
#'                 treatment = factor(sample(c('no','yes'), 100, replace = T, prob = c(0.9,0.1)), levels=c('no','yes')),
#'                 values = rnorm(100, 4, 2))
#' jj_plot_hist(df$values)
#' jj_plot_hist(df, feature = 'values', group = 'groups')
#' jj_plot_hist(df, feature = 'values', group = 'groups', fill = 'letters', nbins = 30)
#' jj_plot_hist(df, feature = 'values', colour='treatment', fill = 'letters', title='my title', group = 'groups', size = 1.3, binwidth = 0.1) + 
#'   scale_fill_manual(values = jj_get_jj_colours(df$letters)) + scale_colour_manual(values = c(yes='black', no = "white"))

#' @export
jj_plot_hist = function(data, feature = NULL, fill = NULL, colour= NULL, group = NULL,
                        title = '',  nbins = 50, size = 1, ...){
  #helper function to plot histograms providing either a data.frame as `data` or a vector of values
  if(is.vector(data)){
    data = data.frame(feature=data)
    feature = 'feature'
  }
  
  gg = ggplot(data, aes_string(x=feature, fill= fill, colour=colour)) +
    geom_histogram(bins = nbins, size = size, ... ) + 
    labs(title = title) + 
    theme_minimal() + 
    theme(plot.title = element_text(hjust = 0.5))
  
  if(!is.null(group)){
    gg = gg + facet_wrap(as.formula(paste("~", group)))
  }
  gg
}
