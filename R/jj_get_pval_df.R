#' get p values for group comparisons with wilcox test
#'
#' get a data.frame of pvalues which can be passed to ggpubr to include in a barplot/boxplot...
#' Each donor can have multiple samples. Each donor should uniquely belong to one group
#' 
#' @param sample_column column containing annotation for the sample from which the observations are derived.  number of columns
#' @param donor_column column specifying the donor from which the observations are derived. If all samples are unpaired, it is the same as `sample_column` but must be passed as a second column to dplyr number of rows
#' @param group_column column specifying the groups that should be used for the statistical comparison
#' @param score_column column containing the values that are used for the comparison between groups (wilcox.test)
#' @param comparisons_paired list of vectors of length 2 specifying the desired paired comparisons between levels in the `group_column`
#' @param comparisons_unpaired list of vectors of length 2 specifying the desired unpaired comparisons between levels in the `group_column`
#' @return Returns a tibble of pvalues
#' @export
#' @examples
#' #make some example data: 10 samples, 4 healthy donors and 3 patients. From each patients, values are available from
#' #two timepoints: diagnosis and relapse. We want to quantify the difference in a `Disease_score` between Healthy, Diagnosis and 
#' #relapse state on the sample level (not on single-cell level). Both paired (Diagnosis - Relapse) and unpaired
#' #(Healthy-Diagnosis, Healthy-Relapse) comparisons need to be performed.
#' data_df = data.frame(Sample = rep(1:10, each = 10),
#'                      Donor = c(rep(c('Healthy1','Healthy2','Healthy3','Healthy4'), each=10), rep(c('Disease1','Disease2','Disease3'), each=20)),
#'                      Group = c(rep('Healthy', 40), rep(rep(c('Diagnosis','Relapse'), each=10), 3)),
#'                      Disease_score = c(rnorm(40, mean = 3, sd = 0.5),  rnorm(20, mean = 4, sd = 0.1),
#'                                        rnorm(10, mean = 6, sd = 0.1), rnorm(10, mean = 4.5, sd = 0.1),
#'                                        rnorm(10, mean = 4.5, sd = 0.1), rnorm(10, mean = 3, sd = 0.1)))
#' #plot the data per group as violins
#' gg = jj_plot_numeric_by_group(data_df, feature_column = 'Disease_score', group_column = 'Group')
#' gg
#' #get mean score per sample
#' sample_df = get_pval_df(data_df,
#'                       sample_column = 'Sample',
#'                       donor_column =  'Donor',
#'                       group_column = 'Group',
#'                       score_column =  'Disease_score',
#'                       return_per_sample_df = T)
#' #i.e. comparisons which are calculated in the function 'get_pval_df'
#' wilcox.test(x=sample_df$score[sample_df$group == 'Healthy'], 
#'             y=sample_df$score[sample_df$group == 'Diagnosis'],
#'             paired = F)
#' wilcox.test(x=sample_df$score[sample_df$group == 'Relapse'], 
#'             y=sample_df$score[sample_df$group == 'Diagnosis'],
#'             paired = T)
#' #get the full pvalue tibble
#' pval_df = jj_get_pval_df(data_df,
#'                       sample_column = 'Sample',
#'                       donor_column =  'Donor',
#'                       group_column = 'Group',
#'                       score_column =  'Disease_score',
#'                       comparisons_paired = list(c('Diagnosis','Relapse')),
#'                       comparisons_unpaired = list(c('Healthy', 'Diagnosis'), c('Healthy', 'Relapse')),
#'                       )
#' pval_df
#' #add suitable y positions for p value visualization in the violin plot
#' pval_df$ypos = c(7, 6.5, 6.75)
#' #plot everything together
#' gg + geom_point(data = sample_df, mapping = aes(x = group, y = score, colour=donor), size=3) + 
#'   scale_colour_manual(values=jj_get_jj_colours(sample_df$donor)) +
#'   ggpubr::stat_pvalue_manual(data = pval_df, label='p_val_adj',
#'                                 tip.length = 0.01,
#'                                 y.position = 'ypos') 

jj_get_pval_df = function(df, sample_column, donor_column, group_column, score_column,
                       comparisons_paired, comparisons_unpaired,
                       p_adj_method='BH', p_adj_round=3, return_per_sample_df = FALSE){
  
  df = df %>% dplyr::select(sample_column, donor_column, group_column, score_column) %>% 
    dplyr::rename(sample = sample_column, donor=donor_column, group=group_column, score=score_column) %>% 
    dplyr::group_by(sample, group) %>% 
    dplyr::summarise(donor = get_mode(donor), #group=get_mode(group), 
                     score = mean(score, na.rm = T))
  if(return_per_sample_df){
    return(df)
  }
  
  #df = df %>% dplyr::select(sample_column, donor_column, group_column, score_cols) %>% as.data.frame
  #stopifnot(all(c(unlist(comparisons_paired), unlist(comparisons_unpaired)) %in% df[, group_column]))

  # df = df %>% dplyr::group_by(!!rlang::sym(sample_column)) %>%
  #   dplyr::summarise(across(c(!!rlang::sym(group_column), !!rlang::sym(donor_column)), head, 1),
  #                    across(score_cols, mean)) #todo: code cleaner

  #initialize empty tibble
  pval_df = tibble(group1='', group2='', p=0) %>% dplyr::filter(p!=0)

  #compute pvalues for unpaired comparisons
  for(j in seq_along(comparisons_unpaired)){
      g1 = df$score[df$group == comparisons_unpaired[[j]][1]]
      g2 = df$score[df$group == comparisons_unpaired[[j]][2]]
      pval_df = pval_df %>% dplyr::add_row(group1=comparisons_unpaired[[j]][1],
                                           group2=comparisons_unpaired[[j]][2],
                                           p=wilcox.test(g1, g2, paired = F)$p.value)
  }
  
  for(j in seq_along(comparisons_paired)){
      gl = df %>% dplyr::filter(group %in% comparisons_paired[[j]]) %>%
        dplyr::select(donor, group, score) %>%
        dplyr::group_by(donor) %>%  dplyr::filter(n() > 1) %>%
        dplyr::arrange(donor) %>% dplyr::ungroup() %>%
        dplyr::group_by(group) %>% dplyr::group_split()
      stopifnot(identical(gl[[1]]$donor, gl[[2]]$donor))
      g1 = gl[[1]]$score
      g2 = gl[[2]]$score
      pval_df = pval_df %>% dplyr::add_row(group1=comparisons_paired[[j]][1],
                                           group2=comparisons_paired[[j]][2],
                                           p=wilcox.test(g1, g2, paired = T)$p.value)
  }

  pval_df = dplyr::mutate(pval_df,
                          p_val_adj = round(p.adjust(p, method=p_adj_method),
                                            p_adj_round))

  return(pval_df)
}
 

# score_boxplot = function(summary_df, sample_column, donor_column, group_column, score_columns, sd_suffix = NULL, paired_comparisons=NULL, unpaired_comparisons=NULL){
#   colour_by = group_column
#   stopifnot(all(c(group_column, donor_column, sample_column, score_columns) %in% colnames(summary_df)))
#   g_list = list()
#   for(i in seq_along(score_columns)){
#     message(i)
#     if(!is.null(sd_suffix)){
#       frdf = summary_df %>% dplyr::select(group_column, donor_column, sample_column, score_columns[i], paste0(score_columns[i], sd_suffix)) %>% 
#         dplyr::rename(score = score_columns[i], sd =  paste0(score_columns[i], sd_suffix))
#     }else{
#       frdf = summary_df %>% dplyr::select(group_column, donor_column, sample_column, score_columns[i]) %>% 
#         dplyr::rename(score = score_columns[i])
#     }
#     
#     #initialize empty tibble
#     pval_df = tibble(group1='', group2='', p=0) %>% dplyr::filter(p!=0)
#     #compute pvalues for unpaired comparisons
#     if(!is.null(unpaired_comparisons)){
#       for(j in seq_along(unpaired_comparisons)){
#         g1 = frdf %>% dplyr::filter(!!sym(group_column) == unpaired_comparisons[[j]][1]) %>% pull(score)
#         g2 = frdf %>% dplyr::filter(!!sym(group_column) == unpaired_comparisons[[j]][2]) %>% pull(score)
#         pval_df = pval_df %>% dplyr::add_row(group1=unpaired_comparisons[[j]][1], 
#                                              group2=unpaired_comparisons[[j]][2],
#                                              p=wilcox.test(g1, g2, paired = F)$p.value)
#       }
#     }
#     #compute pvalues for paired comparisons (only keep paired samples)
#     if(!is.null(paired_comparisons)){
#       for(j in seq_along(paired_comparisons)){
#         gl = frdf %>% dplyr::filter(!!sym(group_column) %in% paired_comparisons[[j]]) %>% 
#           dplyr::select(donor_column, group_column, score) %>% 
#           dplyr::group_by(donor_column) %>%  dplyr::filter(n() > 1) %>% 
#           dplyr::arrange(donor_column) %>% dplyr::ungroup() %>%
#           dplyr::group_by(group_column) %>% dplyr::group_split()
#         stopifnot(identical(gl[[1]]$donor_column, gl[[2]]$donor_column))
#         g1 = gl[[1]]$score 
#         g2 = gl[[2]]$score 
#         pval_df = pval_df %>% dplyr::add_row(group1=paired_comparisons[[j]][1], 
#                                              group2=paired_comparisons[[j]][2],
#                                              p=wilcox.test(g1, g2, paired = T)$p.value)
#       }
#     }
#     
#     # #add x and y positions of pvalues in the plot using rstatix (dummy comparison)
#     # pval_df$comp = paste0(pval_df$group1, pval_df$group2)
#     # pdf_dummy = summary_df %>% t_test(score~group_column) %>% #make dummy df
#     #   add_xy_position(x='group_column', fun="max") %>% dplyr::mutate(comp=paste0(group1,group2)) %>% 
#     #   dplyr::select(comp, xmin, xmax, y.position, groups) 
#     # pdf_dummy2 = summary_df %>% t_test(score~group_column) %>% #make dummy df
#     #   add_xy_position(x='group_column', fun="max") %>% dplyr::mutate(comp=paste0(group2,group1)) %>% 
#     #   dplyr::select(comp, xmin, xmax, y.position, groups)
#     # pdf_dummy = rbind(pdf_dummy, pdf_dummy2)
#     pval_df = pval_df %>% #dplyr::left_join(pdf_dummy, by = 'comp') %>% 
#       dplyr::mutate(p = round(p, 3))
#     return(pval_df)
#     
#     score_range = diff(range(frdf$score))*0.1
#     
#     g_list[[i]] =  ggplot(summary_df, aes_string(group = 'group_column' , x= 'group_column' , y=score_columns[i])) +
#       #geom_boxplot(width=0.3) +
#       geom_violin(width=0.3) + 
#       geom_point(aes_string(colour=colour_by), size=2) + 
#       geom_line(aes_string(x='group_column', y = score_columns[i], group = 'donor_column'), alpha=0.4, color='black') + 
#       #stat_compare_means() + 
#       ggpubr::stat_pvalue_manual(pval_df, label='p', tip.length = 0.01,  
#                                  y.position = c(max(frdf$score)+ 1:nrow(pval_df)*score_range)) +
#       labs(title = score_columns[i],  x = 'group_column group', y = 'score') + 
#       NoLegend() + theme_minimal() + theme(plot.title = element_text(hjust = 0.5, size=10)) +
#       scale_x_discrete(guide = guide_axis(n.dodge = 2))
#     
#     # if(colour_by == 'donor_column'){
#     #   g_list[[i]] = g_list[[i]] + scale_colour_manual(values=jj_get_jj_colours(summary_df$donor_column))
#     # }else if(colour_by == 'group_column'){
#     #   g_list[[i]] = g_list[[i]] + scale_color_manual(values= load_colours(summary_df$group_column, colour_csv = '/omics/groups/OE0436/data2/simonma/projects/mm_scrna/scripts/colour_map.csv'))#c(Healthy='#4daf4a',ID='#e41a1c',CR='#a6cee3',nonCR='#1f78b4',LTS='#377eb8'))
#     # }
#     #stat_compare_means(comparisons = my_comparisons_unpaired, method = "wilcox", paired = F) + 
#     #stat_compare_means(comparisons = my_comparisons_paired, method = "wilcox", paired = T) + 
#     #stat_compare_means() 
#   }
#   return(g_list)
# }
