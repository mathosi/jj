#' arrange ggplots into combined plots of n plots
#'
#' arrange ggplots into combined plots of `n` plots and `cols` columns
#'
#' @param gg_list list of ggplot objects
#' @param n number of ggplots to plot together in one combined plot
#' @param cols number of columns for each combined plot
#' @keywords plot
#' @export
#' @examples
#' gg_list = lapply(1:10, function(x) ggplot(data = data.frame(a=rnorm(10, 0, 5), b=rnorm(10, 0, 5)), aes(x=a,y=b))+geom_point(colour=x))
#' jj_plot_n(gg_list,  n=10, cols = 4)


jj_plot_n = function(gg_list, n=4, cols=2){
  #plot `n` ggplots together in one grid with `cols` columns
  stopifnot(is.numeric(n))
  n = as.integer(n)
  n_plot = seq(1, length(gg_list), n)
  for(i in n_plot){
    do.call(gridExtra::grid.arrange, list(n=gg_list[i:min((i+n-1), length(gg_list))], ncol=cols))
  }
}
