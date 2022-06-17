#' arrange ggplots in a grid
#'
#' plot multiple plots together in a grid based on a list of ggplot objects
#' The following error can occur and is related to problems with the legend. Try to change legend.pos to solve it:
#' Error in grid.Call(C_convert, x, as.integer(whatfrom), as.integer(whatto),  : Viewport has zero dimension(s)
#'
#' @name jj_arrange_ggplots
#' @param gg_list list of ggplot objects
#' @param nplots number of plots to combine together
#' @param cols number of columns for the plot layout
#' @param legend.pos position of legends. One of 'left','right','top','bottom','none'. If null, do not change legend positions
#' @param equal_size_plots if TRUE, make all individual plots the same size
#' @keywords ggplot,plot_grobs,plot_grid
#' @export
#' @examples
#' df = data.frame(a=rnorm(100, 0, 5), b=rnorm(100, 0, 5), d=rbinom(100, 50, 0.3), e = sample(c('A','B', 'C'), 100, replace=T))
#' ggl = list()
#' ggl[[1]] = ggplot(df, aes(x = a, y = b, colour=d)) + geom_point()
#' ggl[[2]] = ggplot(df, aes(x = a, y = b, colour=e)) + geom_point()
#' ggl[[3]] = ggplot(df, aes(x = a)) + geom_histogram()
#' jj_arrange_ggplots(ggl, nplots = 4, cols = 2)
#' jj_arrange_ggplots(ggl, nplots = 4, cols = 2, legend.pos = 'none')

jj_arrange_ggplots = function(gg_list, nplots=4, cols=2, legend.pos = NULL, equal_size_plots=TRUE){
  #plot `nplots` ggplots together in one grid with `cols` columns
  #legend.pos is one of 'left','right','top','bottom','none'. If null, do not change legend positions
  #equal_size_plots if TRUE, make plot size equally large independent from the legend
  if(!all(sapply(gg_list, function(x) 'ggplot' %in% class(x)))){
    stop('All elements in gg_list need to be ggplot objects')
  }
  if(equal_size_plots){
    align = 'hv'
  }else{
    align = 'none'
  }
  if(!is.null(legend.pos)){
    gg_list = lapply(gg_list, function(x) x + theme(legend.position = legend.pos))
  }
  stopifnot(is.numeric(nplots))
  nplots = as.integer(nplots)
  grobs_plot = seq(1, length(gg_list), nplots)
  grob_list = list()
  for(i in grobs_plot){
    grob_list[[as.character(i)]] = do.call(cowplot::plot_grid, list(plotlist=gg_list[i:min((i+nplots-1), length(gg_list))], ncol=cols, align=align))
    #do.call(gridExtra::grid.arrange, list(nplots=gg_list[i:min((i+nplots-1), length(gg_list))], ncol=cols))
  }
  return(grob_list)
}