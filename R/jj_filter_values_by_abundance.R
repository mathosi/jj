#' fitler a list of vectors for those values having min/max abundance
#'
#'
#' @name jj_filter_values_by_abundance
#' @param list_of_vectors list of vectors, which are compared against each other
#' @param max_abundance only return elements that occur maximum 'max_abundance' times in the list (only counting once per vector)
#' @param min_abundance only return elements that occur minimum 'min_abundance' times in the list (only counting once per vector)
#' 
#' example_list = list(l1 = letters[1:7], l2 = letters[1:4], l3 = c( 'd', 'd', letters[4:6]))
#' jj_filter_values_by_abundance(example_list, min_abundance = 3)
#' jj_filter_values_by_abundance(example_list, max_abundance = 2)
#' jj_filter_values_by_abundance(example_list, min_abundance = 2, max_abundance = 3)
 
jj_filter_values_by_abundance = function(list_of_vectors, max_abundance=NULL, min_abundance=NULL){
  if(is.null(max_abundance) & is.null(min_abundance)) stop('Either max_abundance or min_abundance must be specified')
  #if(!is.null(max_abundance) & !is.null(min_abundance)) stop('Either max_abundance or min_abundance must be NULL')
  if(!is.null(max_abundance)){
    if(max_abundance < 1 | max_abundance > length(list_of_vectors)) stop('max_abundance value cannot be larger than the length of the list')
    # olap_df = jj_make_overlap_df(list_of_vectors)$overlaps
    # olap_intersects = olap_df[, nchar(colnames(olap_df)) <= max_abundance ]
  }else{
    if(min_abundance < 1 | min_abundance > length(list_of_vectors)) stop('min_abundance value cannot be larger than the length of the list')
    # olap_df = jj_make_overlap_df(list_of_vectors)$overlaps
    # olap_intersects = olap_df[, nchar(colnames(olap_df)) >= min_abundance ]
  }
  
  all_features = unique(unlist(list_of_vectors))
  features_in_all = lapply(list_of_vectors, function(x) as.integer(all_features %in% 
                                                                     x))
  upsetr_tb = jj_list_to_df(features_in_all, long = F)
  rownames(upsetr_tb) = all_features
  rsums = rowSums(upsetr_tb)
  
  if(!is.null(max_abundance) & !is.null(min_abundance)){
    values_keep = rownames(upsetr_tb)[rsums >= min_abundance & rsums <= max_abundance]
  }else if(!is.null(max_abundance)){
    values_keep = rownames(upsetr_tb)[rsums <= max_abundance]
  }else if(!is.null(min_abundance)){
    values_keep = rownames(upsetr_tb)[rsums >= min_abundance]
  }
  
  olap_list = lapply(list_of_vectors, function(x) x[x %in% values_keep])
  # value_sets = colnames(olap_intersects)[nchar(colnames(olap_intersects)) == 1]
  # olap_list = list()
  # for(i in seq_along(value_sets)){
  #   olap_list[[i]] = unique(na.omit(unlist(as.list(olap_intersects[, grep(value_sets[i], colnames(olap_intersects)), drop = F]))))
  # }
  #names(olap_list) = names(value_list)
  olap_list
}
