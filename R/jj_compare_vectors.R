#' quantify number of common elements between two or more vectors
#'
#' jj_compare_vectors returns a data.frame specifying how many elements of vector A are present in B and vice versa
#' jj_plot_upsetr uses the UpSetR:upset function to plot overlaps between vectors passed as a list
#' jj_make_overlap_df returns a list: key_table mapping the names of the list to letters and overlaps, a data.frame with the elements that are present in each combination of vectors passed as a list
#'
#' @param A first vector
#' @param B second vector
#' @param list_of_vectors list of vectors to be compared
#' @param mainbar.y.label Passed to upset function, default: 'Intersection Size',  
#' @param mb.ratio Passed to upset function, default: c(0.55, 0.45)
#' @export
#' @examples
#' vec_list = list(group1= 1:5, group2= 2:4, group3 = c(3,6), group4=5:10, group5=12)
#' jj_compare_vectors(A=vec_list[[1]], B=vec_list[[2]])
#' jj_plot_upsetr(vec_list)
#' jj_make_overlap_df(vec_list)

#' @export
jj_compare_vectors = function(A, B){
  a_in_b = sum(A %in% B)
  a_not_in_b = length(A) - a_in_b
  b_in_a = sum(B %in% A)
  b_not_in_a = length(B) - b_in_a
  name_a = deparse(substitute(A))
  name_b = deparse(substitute(B))
  comp1 = sprintf('%s in %s', name_a, name_b)
  comp2 = sprintf('%s in %s', name_b, name_a)
  df = data.frame(c(a_in_b, a_not_in_b), c(b_in_a, b_not_in_a), row.names = c('TRUE', 'FALSE'))
  colnames(df) = c(comp1, comp2)
  return(df)
}

#' @export
jj_plot_upsetr = function(list_of_vectors, mainbar.y.label =  'Intersection Size',  mb.ratio = c(0.55, 0.45), 
                          nsets=100, nintersects=40, ...){
  library(UpSetR)
  all_features = unique(unlist(list_of_vectors))
  features_in_all = lapply(list_of_vectors, function(x) as.integer(all_features %in% x))
  upsetr_tb = jj_list_to_df(features_in_all,long = F)
  upsetr = UpSetR::upset(upsetr_tb, nsets=nsets, nintersects = nintersects,  mb.ratio = mb.ratio,
                         sets = rev(colnames(upsetr_tb)), keep.order = T, 
                         mainbar.y.label =  mainbar.y.label, 
                         order.by='freq', set_size.show = TRUE,   
                         set_size.scale_max= nrow(upsetr_tb) +1.5, ...)
  return(upsetr)
}

#' @export
jj_make_overlap_df = function(list_of_vectors){
  #find intersects between different vectors
  if(length(list_of_vectors) > 10){
    stop('Currently a safety stop is implemented, which avoids lists with more than 10 elements resulting in 61 comparisons')
  }
  ##example
  ##returns list with two elements
  # $key_table
  # signature_name key
  # 1         group1   a
  # 2         group2   b
  # 3         group3   c
  # 4         group4   d
  # 5         group5   e
  # 
  # $overlaps
  # a  b  c  d  e ab ac ad bc cd abc
  # 1  1  2  3  5 12  2  3  5  3  6   3
  # 2  2  3  6  6 NA  3 NA NA NA NA  NA
  # 3  3  4 NA  7 NA  4 NA NA NA NA  NA
  # 4  4 NA NA  8 NA NA NA NA NA NA  NA
  # 5  5 NA NA  9 NA NA NA NA NA NA  NA
  # 6 NA NA NA 10 NA NA NA NA NA NA  NA
  
  all_levels_list = list_of_vectors
  names(all_levels_list) = letters[1:length(list_of_vectors)]
  key_df = data.frame(signature_name = names(list_of_vectors), key = names(all_levels_list))
  message(sprintf('%s is encoded by %s.', paste(names(list_of_vectors), collapse = ', '), paste(names(all_levels_list), collapse = ', ')))
  #trib = paste0("~signature_name, ~key,'", paste(paste(names(list_of_vectors), names(all_levels_list), sep = "','"), collapse = "','"),  "'"))
  # last_combo = paste0(names(all_levels_list), collapse = '')
  # while(!last_combo %in% names(all_levels_list)){
  #   all_pairs = combn(names(all_levels_list), 2)
  #   all_pairs = all_pairs[, !apply(all_pairs, 2, function(x) any(unlist(strsplit(x[1], '')) %in% unlist(strsplit(x[2], ''))))]
  #   for(i in seq(ncol(all_pairs))){
  #     all_levels_list[[paste0(all_pairs[1, i], all_pairs[2, i], collapse = '')]] = intersect(all_levels_list[names(all_levels_list)==all_pairs[1, i]][[1]],
  #                                                                                            all_levels_list[names(all_levels_list)==all_pairs[2, i]][[1]])
  #     if(last_combo %in% names(all_levels_list)) break
  #   }
  # }
  
  for(n in 1:length(list_of_vectors)){
    all_pairs = combn(key_df$key, n)
    for(i in seq(ncol(all_pairs))){
      levels_use_list = all_levels_list[names(all_levels_list) %in% all_pairs[, i]]
      isect = levels_use_list[[1]]
      if(length(levels_use_list)>1){
        for(j in 2:length(levels_use_list)){
          isect = intersect(isect, levels_use_list[[j]])
        }
      }
      all_levels_list[[paste0(all_pairs[, i], collapse = '')]] = isect
    }
  }
  
  #make data.frame out of list
  maxLength <- max(unlist(lapply(all_levels_list, length)))
  eqLengthVectorlist <- lapply(all_levels_list, 
                               function(x) c(x, rep(NA, maxLength - length(x))))
  filledDf <- do.call(cbind, eqLengthVectorlist)
  filledDf[is.na(filledDf)] <- NA
  Df <- as.data.frame(filledDf)
  Df = Df[, !apply(Df,2,function(x) all(is.na(x)))]
  res_list = list(key_table = key_df, overlaps = Df)
  return(res_list)
}




# plot_upsetr_from_df = function(zero_one_integer_df){
#   library(UpSetR)
#   #zero_one_integer_df = apply(zero_one_integer_df, 2, as.integer)
#   upsetr = UpSetR::upset(as.data.frame(zero_one_integer_df), nsets=100,  mb.ratio = c(0.55, 0.45),
#                          sets = rev(colnames(zero_one_integer_df)), keep.order = T, 
#                          mainbar.y.label =  'Intersection Size', 
#                          order.by='freq', set_size.show = TRUE,   
#                          set_size.scale_max= nrow(zero_one_integer_df) +1.5 )
#   return(upsetr)
# }

