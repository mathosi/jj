#' arrange ggplots into combined plots of n plots
#'
#' arrange ggplots into combined plots of `n` plots and `cols` columns
#'
#' @param gg_list list of ggplot objects
#' @param nplots number of ggplots to plot together in one combined plot
#' @param cols number of columns for each combined plot
#' @keywords plot
#' @export
#' @examples
#' gg_list = lapply(1:6, function(x) ggplot(data = data.frame(a=rnorm(10, 0, 5), b=rnorm(10, 0, 5)), aes(x=a,y=b))+geom_point(colour=x))
#' jj_combine_plots(gg_list, nplots=6, cols = 4)


jj_combine_plots = function(gg_list, nplots=4, cols=2){
  stopifnot(is.numeric(nplots))
  nplots = as.integer(nplots)
  n_plot = seq(1, length(gg_list), nplots)
  for(i in n_plot){
    do.call(gridExtra::grid.arrange, list(grobs=gg_list[i:min((i+nplots-1), length(gg_list))], ncol=cols))
  }
}
