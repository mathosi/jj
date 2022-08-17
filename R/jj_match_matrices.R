#' match two matrices/data.frames by their row names/column names
#'
#' return the subset of two matrices that shares the same rownames, colnames or both.
#'
#' @param mat_a matrix 1
#' @param mat_b matrix 2
#' @param match_by match by `columns`, `rows` or `all`
#' @param print_only only print out whether rownames or colnames are identical between two matrices
#' @param verbose rownames to set
#' @return Returns a list with the matched features/samples in the same order
#' @export
#' @examples
#' matrix1 = matrix(1:80, nrow = 8, dimnames = list(paste0('gene_',1:8), paste0('cell_', 1:10)))
#' matrix2 = matrix(1:50, nrow = 5, dimnames = list(paste0('gene_',c(1,4,8,2,3)), paste0('cell_', c(2:5, 8:13))))
#' jj_match_matrices(matrix1, matrix2, print_only = T)
#' jj_match_matrices(matrix1, matrix2)

jj_match_matrices = function(mat_a, mat_b, match_by = c('all', 'columns', 'rows')[1], print_only=FALSE, verbose = T){
  names_match = function(a, b){
    if(is.null(rownames(a))) message('a does not have rownames')
    if(is.null(rownames(b))) message('b does not have rownames')
    if(is.null(colnames(a))) message('a does not have colnames')
    if(is.null(colnames(b))) message('b does not have colnames')
    c_a = identical(rownames(a), rownames(b))
    c_b = identical(colnames(a), colnames(b))
    c_c = identical(rownames(a), colnames(b))
    c_d = identical(colnames(a), rownames(b))
    if(c_a) message(sprintf('rownames in a match rownames in b: %s, ...', paste(head(rownames(a), 2), collapse = ', ')))
    if(c_b) message(sprintf('colnames in a match colnames in b: %s, ...', paste(head(colnames(a), 2), collapse = ', ')))
    if(c_c) message(sprintf('rownames in a match colnames in b: %s, ...', paste(head(rownames(a), 2), collapse = ', ')))
    if(c_d) message(sprintf('colnames in a match rownames in b: %s, ...', paste(head(colnames(a), 2), collapse = ', ')))
    if(!any(c(c_a, c_b, c_c, c_d))) message('a and b do not have identical columnnames or rownames')
  }
  if(print_only){
    names_match(mat_a, mat_b)
    return()
  }
  
  if(match_by %in% c('all','columns')){
    keep_a = colnames(mat_a) %in% colnames(mat_b)
    if(verbose) message('Dropping ', sum(!keep_a), ' columns from mat_a')
    mat_a = mat_a[, keep_a]
    keep_b = colnames(mat_b) %in% colnames(mat_a)
    if(verbose) message('Dropping ', sum(!keep_b), ' columns from mat_b')
    mat_b = mat_b[, keep_b]
    mat_a = mat_a[, match(colnames(mat_b), colnames(mat_a))]
    stopifnot(identical(colnames(mat_a), colnames(mat_b)))
  }
  if(match_by %in% c('all','rows')){
    keep_a = rownames(mat_a) %in% rownames(mat_b)
    if(verbose) message('Dropping ', sum(!keep_a), ' rows from mat_a')
    mat_a = mat_a[keep_a, ]
    keep_b = rownames(mat_b) %in% rownames(mat_a)
    if(verbose) message('Dropping ', sum(!keep_b), ' rows from mat_b')
    mat_b = mat_b[keep_b, ]
    mat_a = mat_a[match(rownames(mat_b), rownames(mat_a)), ]
    stopifnot(identical(rownames(mat_a), rownames(mat_b)))
  }
  out = list(mat_a = mat_a, mat_b = mat_b)
  if(verbose) {
    j(out)
  }
  return(out)
}
