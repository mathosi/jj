#' create data.frame of defined size with initial value in each cell
#'
#' Creates a data.frame of size `ncol` * `nrow` filled with `init`. `row.names` and `column.names` may also be specified.
#'
#' @param ncol number of columns
#' @param nrow number of rows
#' @param init value to put into every cell
#' @param col.names colnames to set
#' @param row.names rownames to set
#' @param return_matrix return matrix instead of data.frame, default=FALSE
#' @keywords create data.frame, initialize data.frame, dummy data.frame
#' @export
#' @examples
#' jj_initialize_df(ncol = 3, nrow = 4, init = 0, col.names = paste0("col", 1:3), row.names = paste0("row", 1:4))


jj_initialize_df <- function(ncol, nrow, init=NA, col.names=NULL,return_matrix=FALSE, sparse=FALSE, ...){
  init_df <- data.frame(matrix(init, ncol=ncol, nrow=nrow),stringsAsFactors = FALSE,...)
  if(!is.null(col.names) & length(col.names) == ncol){
    colnames(init_df) <- col.names
  }
  if(return_matrix){
    init_df = as.matrix(init_df)
    if(sparse){
      library(Matrix)
      init_df = as(df, 'sparseMatrix')
    }
  }
  return(init_df)
}
