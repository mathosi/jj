#' encirle fraction of points that contain x percent of the density mass
#'
#' estimate densities to include density area in tsne/umap
#'
#' @param df data.frame containing group_column, value_column, sample_column
#' @param group_column column containing the grouping variable (i.e. the x-axis)
#' @param value_column column containing the categorical data to summarize per group. Within each group, the relative count of observations for each level in the value_column is summarized
#' @param sample_column column containing annotation for the sample from which the observations are derived. Used to calculate the SEM
#' @param custom_colors optional named vector with custom colours
#' @param return_df return the summarized data as data.frame instead of a ggplot 
#' @param show_error_bars if TRUE, plot the SEM error bars
#' @param show_fraction if TRUE, add the fraction of observations as a number to each bar
#' @param error_bar_color if not NULL, overwrite the standard colour of the error bars
#' @return Returns a stacked barplot with error bars indicating the SEM
#' @export
#' @examples
#' pbmc_small$sample = sample(1:3, ncol(pbmc_small), replace = T)
#' jj_plot_stacked_barplot(df = pbmc_small[[]], group_column = 'groups', 
#'                         value_column = 'RNA_snn_res.1', sample_column = 'sample', 
#'                         show_error_bars = T, show_fraction = T, return_df = F, 
#'                         error_bar_color = 'black')
#' jj_plot_stacked_barplot(df = pbmc_small[[]], group_column = 'groups',
#'                         value_column = 'RNA_snn_res.1', sample_column = 'sample',
#'                         show_error_bars = T, show_fraction = F, return_df = F,
#'                         custom_colors = jj_get_jj_colours(levels(pbmc_small$RNA_snn_res.1)))

jj_plot_stacked_barplot = function(df, group_column, value_column, sample_column, custom_colors=NULL,
                           return_df = FALSE, show_error_bars=TRUE, show_fraction =TRUE, error_bar_color=NULL){
  
  n_sample_df = df %>%
    dplyr::select(!!rlang::sym(group_column), !!rlang::sym(sample_column)) %>% 
    dplyr::group_by(!!rlang::sym(group_column)) %>% 
    dplyr::summarise(n_samples= n_distinct(!!rlang::sym(sample_column)))
  #return(n_sample_df)
  
  cell_sem_df = df %>%
    dplyr::select(!!rlang::sym(group_column), !!rlang::sym(sample_column),!!rlang::sym(value_column)) %>% 
    dplyr::left_join(n_sample_df, by = group_column) %>% 
    dplyr::group_by_all() %>% 
    dplyr::summarise(ncells= n()) %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(!!rlang::sym(group_column),  !!rlang::sym(sample_column)) %>% 
    dplyr::mutate(fraction = ncells / sum(ncells)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!rlang::sym(group_column), !!rlang::sym(value_column)) %>% 
    dplyr::summarise(sem_fraction=sd(fraction)/sqrt(n_samples),
                     fraction=sum(fraction) / n_samples,
                     n_samples=n_samples) %>% 
    dplyr::distinct() %>% 
    dplyr::ungroup() %>% 
    dplyr::group_by(!!rlang::sym(group_column)) %>%
    dplyr::arrange(desc(!!rlang::sym(value_column))) %>% 
    dplyr::mutate(cumfrac=cumsum(fraction)) %>% 
    dplyr::ungroup()
  
  cell_sem_df = as.data.frame(cell_sem_df)
  cell_sem_df$cumfrac_sem_frac = cell_sem_df$cumfrac + cell_sem_df$sem_fraction
  
  if(return_df){
    return(cell_sem_df)
  }
  
  gg = ggplot() +
    geom_bar(data = cell_sem_df, mapping = aes_string(x=group_column, fill = value_column, y = 'fraction'), 
             stat='identity', position="fill") + 
    theme_minimal() + 
    scale_y_continuous(breaks=seq(0,1,0.2)) + 
    labs(fill = value_column, colour=value_column) + 
    theme(panel.grid.major.x = element_blank())
  
  if(show_error_bars){
    for(i in unique(cell_sem_df[, group_column])){
      dat_use = cell_sem_df[cell_sem_df[, group_column]== i, ]
      dat_use$gr = i
      if(!is.null(error_bar_color)){
        gg = gg + geom_segment(data  = dat_use,
                               mapping = aes_string(x='gr',xend='gr', y='cumfrac',yend='cumfrac_sem_frac'), color=error_bar_color, 
                               arrow = arrow(angle=90, length = unit(0.05, "inches")))
      }else{
        gg = gg + geom_segment(data  = dat_use,
                               mapping = aes_string(x='gr',xend='gr', y='cumfrac',yend='cumfrac_sem_frac', color=value_column), 
                               arrow = arrow(angle=90, length = unit(0.05, "inches")))
      }
      
    }
  }
  
  if(show_fraction){
    cell_sem_df$rounded_fraction = round(cell_sem_df$fraction, 2)
    gg = gg + geom_text(data = cell_sem_df, aes_string(x=group_column, y='fraction', label = 'rounded_fraction'), 
                        position = position_stack(vjust = 0.5), size = 4) 
  }
  
  if(!is.null(custom_colors)){
    gg = gg + 
      scale_fill_manual(values=custom_colors) + 
      scale_color_manual(values=custom_colors) 
    
  }
  return(gg)
}

