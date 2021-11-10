#' get a set of predefined colours for a vector
#'
#' Returns a vector of colors with names according to input vector levels.
#'
#' @param vec vector containing the names, for which colours should be assinged
#' @param symbols if TRUE, assign numbers instead of colors (useful if groups should be encoded by symbols)
#' @keywords colours
#' @export
#' @examples
#' jj_get_jj_colours(c("a","d","c","a"))
#' jj_get_jj_colours(c("a","d","c","a"), symbols=TRUE)

jj_get_jj_colours <- function(vec, symbols=FALSE){
  if(symbols){
    myColors <- c(0, 4, 2, 7,10,
                  16, 12,13, 14, 11,
                  8, 1, 15,17,9,
                  23,6,20, 5,18,
                  19,21,22,24,3)
  }else{
    # qual_col_pals = RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
    # col_vector = unlist(mapply(RColorBrewer::brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    # greys <- grDevices::colorRampPalette(c("grey20", "grey80"))
    # myColors <- c(RColorBrewer::brewer.pal(9, "Set1")[-6], RColorBrewer::brewer.pal(8, "Set2")[c(1,7)],
    #               RColorBrewer::brewer.pal(9, "YlOrRd")[9], RColorBrewer::brewer.pal(9,"RdPu")[9], 
    #               RColorBrewer::brewer.pal(9,"BuGn")[9],
    #               "darkolivegreen2", "cyan", "green", "magenta", "gold", "pink", "black",
    #               RColorBrewer::brewer.pal(9,"BrBG")[2:4], RColorBrewer::brewer.pal(9,"YlOrRd")[c(1,3)], "yellow", "blueviolet",
    #               "cyan4", "burlywood4", "deeppink3", "lightcyan2", "mediumspringgreen",
    #               "darkblue", "orangered3", "royalblue", "plum","turquoise4", "yellowgreen",
    #               "yellow4", "thistle", "lightslateblue",grDevices::rainbow(200))
    myColors = c('#E41A1C', '#377EB8', '#4DAF4A', '#984EA3', '#FF7F00', '#A65628', '#F781BF', '#999999', '#66C2A5', '#E5C494', '#800026', 
                 '#49006A', '#00441B', 'darkolivegreen2', 'cyan', 'green', 'magenta', 'gold', 'pink', 'black', '#BF812D', '#DFC27D', 
                 '#F6E8C3', '#FFFFCC', '#FED976', 'yellow', 'blueviolet', 'cyan4', 'burlywood4', 'deeppink3', 'lightcyan2', 'mediumspringgreen', 
                 'darkblue', 'orangered3', 'royalblue', 'plum', 'turquoise4', 'yellowgreen', 'yellow4', 'thistle', 'lightslateblue', '#FF0000', 
                 '#FF0800', '#FF0F00', '#FF1700', '#FF1F00', '#FF2600', '#FF2E00', '#FF3600', '#FF3D00', '#FF4500', '#FF4D00', '#FF5400', '#FF5C00',
                 '#FF6300', '#FF6B00', '#FF7300', '#FF7A00', '#FF8200', '#FF8A00', '#FF9100', '#FF9900', '#FFA100', '#FFA800', '#FFB000', '#FFB800',
                 '#FFBF00', '#FFC700', '#FFCF00', '#FFD600', '#FFDE00', '#FFE500', '#FFED00', '#FFF500', '#FFFC00', '#FAFF00', '#F2FF00', '#EBFF00', 
                 '#E3FF00', '#DBFF00', '#D4FF00', '#CCFF00', '#C4FF00', '#BDFF00', '#B5FF00', '#ADFF00', '#A6FF00', '#9EFF00', '#96FF00', '#8FFF00', 
                 '#87FF00', '#80FF00', '#78FF00', '#70FF00', '#69FF00', '#61FF00', '#59FF00', '#52FF00', '#4AFF00', '#42FF00', '#3BFF00', '#33FF00', 
                 '#2BFF00', '#24FF00', '#1CFF00', '#14FF00', '#0DFF00', '#05FF00', '#00FF03', '#00FF0A', '#00FF12', '#00FF1A', '#00FF21', '#00FF29', 
                 '#00FF30', '#00FF38', '#00FF40', '#00FF47', '#00FF4F', '#00FF57', '#00FF5E', '#00FF66', '#00FF6E', '#00FF75', '#00FF7D', '#00FF85', 
                 '#00FF8C', '#00FF94', '#00FF9C', '#00FFA3', '#00FFAB', '#00FFB3', '#00FFBA', '#00FFC2', '#00FFC9', '#00FFD1', '#00FFD9', '#00FFE0',
                 '#00FFE8', '#00FFF0', '#00FFF7', '#00FFFF', '#00F7FF', '#00F0FF', '#00E8FF', '#00E0FF', '#00D9FF', '#00D1FF', '#00C9FF', '#00C2FF',
                 '#00BAFF', '#00B2FF', '#00ABFF', '#00A3FF', '#009CFF', '#0094FF', '#008CFF', '#0085FF', '#007DFF', '#0075FF', '#006EFF', '#0066FF', 
                 '#005EFF', '#0057FF', '#004FFF', '#0047FF', '#0040FF', '#0038FF', '#0030FF', '#0029FF', '#0021FF', '#0019FF', '#0012FF', '#000AFF', 
                 '#0003FF', '#0500FF', '#0D00FF', '#1400FF', '#1C00FF', '#2400FF', '#2B00FF', '#3300FF', '#3B00FF', '#4200FF', '#4A00FF', '#5200FF', 
                 '#5900FF', '#6100FF', '#6900FF', '#7000FF', '#7800FF', '#8000FF', '#8700FF', '#8F00FF', '#9600FF', '#9E00FF', '#A600FF', '#AD00FF', 
                 '#B500FF', '#BD00FF', '#C400FF', '#CC00FF', '#D400FF', '#DB00FF', '#E300FF', '#EB00FF', '#F200FF', '#FA00FF', '#FF00FC', '#FF00F5', 
                 '#FF00ED', '#FF00E6', '#FF00DE', '#FF00D6', '#FF00CF', '#FF00C7', '#FF00BF', '#FF00B8', '#FF00B0', '#FF00A8', '#FF00A1', '#FF0099', 
                 '#FF0091', '#FF008A', '#FF0082', '#FF007A', '#FF0073', '#FF006B', '#FF0063', '#FF005C', '#FF0054', '#FF004C', '#FF0045', '#FF003D',
                 '#FF0036', '#FF002E', '#FF0026', '#FF001F', '#FF0017', '#FF000F', '#FF0008')
  }
  unique_vals <- as.character(unique(vec))
  n_unique_vals <- length(unique_vals)
  if(n_unique_vals > 50){
    warning(paste0("number of levels is ", n_unique_vals, ", but only 50 distinguishable colors are available."))
  }else{
    names(myColors) <- unique_vals
    myColors <- myColors[1:n_unique_vals]
  }
  #quick fix: define mistyrose for 'NA' values
  #myColors <- c("NA"="mistyrose", myColors)[!duplicated(c('NA',names(myColors)))]
  return(myColors)
}
