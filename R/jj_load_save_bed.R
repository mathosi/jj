#' write Granges or region string vector to bed file
#'
#' estimate densities to include density area in tsne/umap
#'
#' @param regions regions to save to bed file, either as GRanges or as strings in the format 'chr1-1-99'
#' @param file_name file name to save to or load from
#' @param return_granges if TRUE, return a GRanges object instead of a tibble when loading data
#' @export
#' @examples
#' jj_save_bed(c('chr1-234-3488', 'chrX-300-1500'), 'test.bed')
#' jj_load_bed('test.bed', return_granges = TRUE)


jj_save_bed <- function(regions, file_name){
  if(!is(regions, 'GRanges')){
    regions <- data.frame(regions = regions)
    regions <- tidyr::separate(
      data = regions,
      col = "regions",
      sep = "-",
      into = c("chr", "start", "end")
    )
    regions <- GenomicRanges::makeGRangesFromDataFrame(regions)
  }
  regions_df <- as.data.frame(regions)
  regions_df$name <- with(regions_df, paste(seqnames, start, end, sep = '-'))
  regions_df <- dplyr::select(regions_df, seqnames, start, end, name)
  readr::write_tsv(regions_df, file = file_name, col_names = F)
}

#' @export
jj_load_bed <- function(file_name, return_granges = FALSE, col_names = c('seqnames','start','end','name')){
  regions = readr::read_tsv(file = file_name, col_names = col_names)
  if(return_granges){
    regions = GenomicRanges::makeGRangesFromDataFrame(regions, keep.extra.columns = TRUE)
  }
  return(regions)
}
