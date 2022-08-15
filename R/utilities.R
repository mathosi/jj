

min_max_normalize <- function(x)(x- min(x)) /(max(x)-min(x))

#' Convert a genomic coordinate string to a GRanges object
#'
#' @param regions Vector of genomic region strings
#' @param sep Vector of separators to use for genomic string. First element is
#' used to separate chromosome and coordinates, second separator is used to
#' separate start and end coordinates.
#' @param ... Additional arguments passed to
#' \code{\link[GenomicRanges]{makeGRangesFromDataFrame}}
#' @return Returns a GRanges object
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom tidyr separate
#' @examples
#' regions <- c('chr1-1-10', 'chr2-12-3121')
#' StringToGRanges(regions = regions)
#' @export
#' @concept utilities
from_to = function(annot_vec, old_new_map_vec){
  #map values from old to new
  #from_to(annot_vec = c('T_cell', 'B_cells_2', 'B_cell_1'), 
  #old_new_map_vec = c('B_cells_1' = 'B cell','B_cells_2' = 'B cell'))
  new_annot = plyr::mapvalues(annot_vec, 
                              from=names(old_new_map_vec),
                              to=old_new_map_vec)
  return(new_annot)
}


zip_lists = function(list1, list2){
  #zip two lists of equal length similar to pythons zip()
  #zip_lists(list('A', 'B'), list(1,2))
  #[[1]][1] "A" [[2]][1] 1 [[3]][1] "B" [[4]][1] 2
  stopifnot(length(list1)==length(list2))
  comb_list = list()
  i=j=1
  while(j <= length(list1)) {
    comb_list[[i]] = list1[[j]]
    i=i+1
    comb_list[[i]] = list2[[j]]
    i=i+1
    j=j+1
  }
  return(comb_list)
}

zip_vectors = function(v1, v2){
  stopifnot(length(v1)==length(v2))
  comb_v = vector()
  i=j=1
  while(j <= length(v1)) {
    comb_v[i] = v1[j]
    i=i+1
    comb_v[i] = v2[j]
    i=i+1
    j=j+1
  }
  return(comb_v)
}

zip_lists2 = function(listoflists){
  #zip multiple lists of equal length into one
  #test_list = list(a=list(c('A'), c('AA'), c('AAA')), b=list(c('B'), c('BB'), c('BBB')), c=list(c('C'), c('CC'), c('CCC')))
  #zip_lists2(test_list)
  #result: c('A'), c('B'), c('C'), c('AA'), ...
  lllen = length(listoflists)
  llen = unique(sapply(listoflists, length))
  stopifnot(length(llen)==1)
  comb_list = vector(mode = "list", length = lllen*llen)
  j = 1
  for(i in 1:llen){
    k = j + lllen -1
    comb_list[j:k] =  lapply(listoflists, '[[', i)
    j = k + 1
  }
  return(comb_list)
}

matrix_to_vector = function(mat, by = 'column'){
  #convert matrix to vector (columnwise)
  vals_vec = c()
  by = match.arg(by, choices = c('column', 'row'))
  if(by =='row'){
    mat = t(mat)
  }
  for (i in seq(1:ncol(mat))){
    vals_vec = c(vals_vec, mat[, i])
  }
  return(vals_vec)
}

bin_vector = function (numeric_vector, breaks = 10,
                       binLabels = NULL, return_df = TRUE, ...) 
{
  if (is.null(binLabels)) {
    if (length(breaks) == 1) {
      breaks <- seq(min(numeric_vector), max(numeric_vector), 
                    length.out = breaks + 1)
    }
    binLabels <- character()
    for (i in 1:(length(breaks) - 1)) {
      binLabels[i] <- paste0(breaks[i], "-", breaks[i + 
                                                      1])
    }
  }
  numeric_vector <- as.numeric(as.character(numeric_vector))
  binned <- cut(numeric_vector, breaks = breaks, labels = binLabels, 
                include.lowest = TRUE)
  
  if(return_df){
    binned_df = data.frame(lower_end=as.numeric(gsub('^([0-9.]+)-[0-9.]+$', '\\1', as.character(binned))), 
                           upper_end=as.numeric(gsub('^[0-9.]+-([0-9.]+)$', '\\1', as.character(binned))), 
                           value=numeric_vector,
                           bin=binned)
    return(binned_df)
  }
  return(binned)
}

mixsort_factor = function(vec){
  vec = as.character(vec)
  vec = factor(vec, levels = gtools::mixedsort(unique(vec)))
  return(vec)
}

get_mode <- function(v) {
  uniqv <- unique(na.omit(v))
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

#example_vec = sample(letters[1:4], 1000,replace=T, prob=c(0.5,0.3,0.1,0.1))
#table(example_vec)
#replaced_vec = replace_if(example_vec, count_below=200, replacement = 'rare letter')
#head(replaced_vec)
#table(replaced_vec)
replace_if = function(vec, count_below=NULL, replacement){
  count_df = as.data.frame(table(vec))
  colnames(count_df) = c('Var', 'Freq')
  if(!is.null(count_below)){
    stopifnot(is.numeric(count_below))
    replace_groups = as.character(count_df$Var[count_df$Freq < count_below])
    n_replace = sum(count_df$Freq[count_df$Freq < count_below])
  }
  if(length(replace_groups)>0){
    message(sprintf('Replacing: %s (%i values) by %s', paste(replace_groups, collapse=', '), n_replace, replacement))
    vec[vec %in% replace_groups] = replacement
  }
  return(vec)
}

