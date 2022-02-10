#' Load/save excel workbooks using openxlsx library
#'
#' load excel workbook as list of data.frames or save list of data.frames as excel
#'
#' @name load_save_excel
#' @param file_name Full path to the file to read or write
#' @param list_of_df list of data.frames to write into excel workbook
#' @param sheet_names Names of sheets to read or write, default: all sheets when loading, names(list) when saving
#' @param row_names Set first column as rownames or write rownames from each df, default: FALSE
#' @param ... further arguments passed to read.xlsx or write.xlsx
#' @keywords excel
#' @export
#' @examples
#' 

#' @rdname load_save_excel
#' @export
jj_load_excel = function(file_name, sheet_names= NULL, row_names=FALSE, ...){
  library(openxlsx)
  if(is.null(sheet_names)){
    sheet_names = openxlsx::getSheetNames(file_name)
  }
  excel_list <- list()
  for (i in seq_along(sheet_names)) {
    excel_list[[i]] <- openxlsx::read.xlsx(file_name, sheet = sheet_names[i], ...)
  }
  names(excel_list) <- sheet_names
  return(excel_list)
}

#' @rdname load_save_excel
#' @export
jj_save_excel = function(list_of_df, file_name, sheet_names = NULL, row_names=FALSE, ...){
  library(openxlsx)
  stopifnot(is.list(list_of_df))
  stopifnot(all(sapply(list_of_df, is.data.frame)))
  if(!is.null(sheet_names)){
    stopifnot(length(sheet_names) == length(list_of_df))
    names(list_of_df) <- sheet_names
  } 
  if(!endsWith(file_name, '.xlsx')){
    file_name = paste0(file_name, '.xlsx')
  }
  write.xlsx(x = list_of_df, file = file_name, rowNames = row_names, ...)
}