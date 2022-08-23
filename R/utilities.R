#' remap values
#'
#' pass a named vector to replace values in a second vector
#'
#' @param vec Vector of genomic region strings
#' @param old_new_map_vec named vector with names and values corresponding to old + new pairs of values
#' @param ... Additional arguments passed to
#' @export
#' @examples
#' #map values from old to new
#' from_to(vec = c('T_cell', 'B_cells_2', 'B_cell_1'),  
#' old_new_map_vec = c('B_cells_1' = 'B cell','B_cells_2' = 'B cell'))

from_to = function(vec, old_new_map_vec){

  new_annot = plyr::mapvalues(vec, 
                              from=names(old_new_map_vec),
                              to=old_new_map_vec)
  return(new_annot)
}

#' min-max normalize
#'
#' @param x Vector of genomic region strings
#' @export

min_max_normalize <- function(x)(x- min(x)) /(max(x)-min(x))

#' zip together two lists into one
#'
#' zip two (or more) lists of equal length similar to pythons zip()
#'
#' @param list1 list 1
#' @param list2 list 2
#' @param listoflists for zip_lists_multi: a list containing multiple lists
#' @returns One list with alternating values from either list
#' @export
#' @examples
#' zip_lists(list('A', 'B'), list(1,2))
#' zip_lists_multi(list(list('A', 'B'), list(1,2), list('aa','bb')))

zip_lists = function(list1, list2){
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

#' @export
zip_lists_multi = function(listoflists){
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

#' zip vectors into one
#'
#' Combine two (or more) vectors into one in alternating order
#'
#' @param v1 vector 1
#' @param v2 vector 2
#' @param listofvectors list of vectors of equal length to zip
#' @export
#' @examples
#' zip_vectors(v1 = letters[1:3], v2 = LETTERS[1:3])
#' test_list = list(a=c('A','AA','AAA'), b=c('B','BB','BBB'), c=c('C', 'CC', 'CCC'))
#' zip_vectors_multi(test_list)

zip_vectors = function(v1, v2){
  stopifnot(identical(length(v1), length(v2)))
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

#' @export
zip_vectors_multi = function(listofvectors){
  lllen = length(listofvectors)
  llen = unique(sapply(listofvectors, length))
  stopifnot(length(llen)==1)
  comb_vec = vector()
  j = 1
  for(i in 1:llen){
    k = j + lllen -1
    comb_vec[j:k] =  sapply(listofvectors, '[[', i)
    j = k + 1
  }
  return(comb_vec)
}

#' convert a matrix to a 1d vector
#'
#' @param mat matrix
#' @param by one of `column`, `row`: concatenate column- or rowwise
#' @export
#' @examples
#' test_mat = matrix(1:12, ncol=4, byrow = FALSE)
#' test_mat
#' matrix_to_vector(test_mat)
#' matrix_to_vector(test_mat, by ='row')

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

#' bin a numeric vector and return a factor of bins
#'
#' @param numeric_vector 
#' @param breaks either a single number specifying the number of equally widthed breaks, or a vector of break boundaries
#' @param binLabels optional labels for the new bins that are returned instead of the bin range
#' @param return_df return a data.frame with boundaries, the value and the bin level
#' @export
#' @examples
#' vec = sample(1:100, size = 100)
#' bin_vector(vec, breaks = 10)
#' bin_vector(vec, breaks = c(0, 33, 66, 100))
#' bin_vector(vec, breaks = c(0, 33, 66, 100), 
#'            binLabels = c('lower tertile', 'middle tertile', 'upper tertile'),
#'            return_df = FALSE)

bin_vector = function (numeric_vector, breaks = 10,
                       binLabels = NULL, return_df = TRUE){
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

#' alphanumerically sort levels of a factor
#'
#' Return a factor with alphabetically and numerically sorted levels
#'
#' @param vec vector to sort
#' @export
#' @examples
#' mixsort_factor(c('Sample_3', 'Sample_1', 'Sample_3', 'Sample_2'))

mixsort_factor = function(vec){
  vec = as.character(vec)
  vec = factor(vec, levels = gtools::mixedsort(unique(vec)))
  return(vec)
}

#' get the mode
#'
#' On equal numbers of occurrence return first item
#' 
#' @param vec 
#' @export
#' @examples
#' get_mode(c('A','B','D','D','A'))

get_mode <- function(vec) {
  uniqv <- unique(na.omit(vec))
  uniqv[which.max(tabulate(match(vec, uniqv)))]
}

#' replace values that meet a condition
#'
#' Currently only replacing values below a specified abundance is implemented
#'
#' @param vec vector with values to replace
#' @param count_below replace values that occur less often than this threshold
#' @param replacement value that is used for replacement
#' @export
#' @examples
#' example_vec = sample(letters[1:4], 1000,replace=T, prob=c(0.5,0.3,0.1,0.1))
#' table(example_vec)
#' replaced_vec = replace_if(example_vec, count_below=200, replacement = 'rare letter')
#' head(replaced_vec)
#' table(replaced_vec)

replace_if = function(vec, count_below=NULL, replacement){
  count_df = as.data.frame(table(vec))
  colnames(count_df) = c('Var', 'Freq')
  replace_groups = NULL
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

#' Prime factorization
#'
#' Returns the prime numbers as vector
#'
#' @param x single number, which is to be split into prime numbers
#' @export
#' @examples
#' jj_prime_factorization(45)

jj_prime_factorization=function(x){
  if(length(x) > 1){
    warning('Only first entry of x will be used')
    x = x[1]
  }
  stopifnot(is.numeric(x))
  n=c()
  i=2
  r=x
  while(prod(n)!=x){
    if(!r%%i) {n=c(n,i);r=r/i;i=1}
    i=i+1
  }
  n
}