#' get colours as named vector from file
#'
#'
#' @param annotation_vector vector with the groups for which colours are needed
#' @param colour_csv two column csv with annotation in first column and colour in second column (alternatively a data.frame with columns `group` and `colour` or a named vector of colours)
#' @param plot_colours visualize colours as ggplot heatmap, default = FALSE
#' @param comment_char Lines in the colour file starting with this character are ignored, default = '/'
#' @export
#' @examples
#' #Save colour mapping in file (may contain comment lines starting with `/`)
#' colour_df  = data.frame(group= c('g1','g2','g3'),colour= c('yellow', 'darkred','#D7B5D8'))
#' write.table(colour_df, file = 'colour_map.csv',
#'             sep = ',', row.names = F, col.names=F)
#' # just show the colours in the file (also works passing the data.frame)
#' jj_get_colours(colour_csv='colour_map.csv', plot_colours=TRUE)
#' jj_get_colours(colour_csv=colour_df, plot_colours=TRUE)
#' # use the colours in ggplot (automatically selects the groups which are needed)
#' custom_cols = jj_get_colours(pbmc_small$groups, colour_csv='colour_map.csv')
#' Seurat::DimPlot(pbmc_small,  group.by = 'groups') + scale_colour_manual(values=custom_cols)
#' #or using jj_plot_features
#' jj_plot_features(pbmc_small, reduction='tsne', meta_features='groups', custom_colors = custom_cols, pt.size = 3)

jj_get_colours = function(annotation_vector=NULL, colour_csv, plot_colours=FALSE, comment_char="/"){
  if(is.data.frame(colour_csv)){
    stopifnot(identical(colnames(colour_csv), c('group', 'colour')))
    col_df = as.data.frame(colour_csv)
  }else if(is.character(colour_csv) & !is.null(names(colour_csv))){
    col_df = data.frame(group=names(colour_csv), colour=colour_csv)
  }else{
    col_df = read.csv(colour_csv, header = F, col.names = c('group', 'colour'), strip.white = T, blank.lines.skip = T, comment.char = comment_char)
  }
  stopifnot(!anyDuplicated(col_df$group))
  col_df$colour = toupper(col_df$colour)
  col_vec = structure(col_df$colour, names=col_df$group)
  if(plot_colours){
    colnr = floor(sqrt(nrow(col_df)))
    rownr = ceiling(nrow(col_df)/colnr)
    col_df$rownr = rep(1:rownr, each = colnr)[1:nrow(col_df)]
    col_df$colnr = rep(1:colnr, length.out = nrow(col_df))
    col_df$group_label = paste(col_df$group, col_df$colour, sep='\n')
    print(ggplot2::ggplot(col_df, aes(x=colnr, y=rownr, fill = group)) + 
            geom_tile() + scale_y_reverse() +
            scale_fill_manual(values= col_vec) + geom_text(aes(label=group_label)) + theme_void() + theme(legend.position = 'none'))
  }
  if(is.null(annotation_vector)){
    return(col_vec)
  }
  annotation_vector = unique(annotation_vector)
  if(!all(annotation_vector %in% names(col_vec))){
    stop(sprintf('%s not present in the colour_csv', paste(annotation_vector[!annotation_vector %in% names(col_vec)], collapse = ', ')))
  }
  col_vec = col_vec[names(col_vec) %in% annotation_vector]
  return(col_vec)
}
