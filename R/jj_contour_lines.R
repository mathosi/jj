#' encirle fraction of points that contain x percent of the density mass
#'
#' estimate densities to include density area in tsne/umap.
#' May fail for small datasets
#'
#' @param reduction_df number of columns
#' @param grouping_var number of rows
#' @param cont_thres value to put into every cell
#' @returns 
#' @export
#' @examples
#' dr_df = data.frame(dim1 = c(rnorm(500, 3, 0.5), rnorm(200, 4, 0.1)),
#'                   dim2 = c(rnorm(500, 3, 0.5), rnorm(100, 4, 1), rnorm(100, 3.5, 1)),
#'                   group = rep(LETTERS[1:2], c(500,200)))
#' cont_l = jj_get_contour_lines(reduction_df = dr_df, grouping_var = 'group', cont_thres = 0.75)
#' gg = jj_plot_features(reduction = dr_df, meta_features = 'group', return_gg_object = T)
#' gg[[1]] + geom_path(data = cont_l,
#'                    aes(x=x, y=y, group=cont_group, linetype=group, colour=group),
#'                    size=1) + 
#'  labs(linetype='')

jj_get_contour_lines <- function(reduction_df, grouping_var, cont_thres){
  
  reduction_df <- reduction_df[!is.na(reduction_df[, grouping_var]), ]
  ls <- getContourLines(reduction_df, grouping_var, cont_thres)
  cont_df <- getContourCoordinates(ls)
  names(cont_df)[5] <- grouping_var
  return(cont_df)
}

getLevel2 <- function(x,y,prob){
  kd <- ks::kde(data.frame(x=x,y=y), compute.cont=TRUE)
  kd$cont[paste0(100-prob*100,'%')]
}

#get the coordinates of the contours containing cont_thres fraction of cells
getContourLines <- function(reduction_df, grouping_var, cont_thres){
  ls <- list()
  for(pop in unique(reduction_df[,grouping_var])){
    reduction_group_df <- reduction_df[reduction_df[, grouping_var] == pop, ]
    cont_level <- getLevel2(reduction_group_df[, 1],
                            reduction_group_df[, 2],
                            cont_thres)
    
    dens <- MASS::kde2d(reduction_group_df[, 1],
                        reduction_group_df[, 2],
                        n=300, #more lines, more accurate, but longer run time
                        lims=c(c(-12, 12), c(-12, 12)))  # don't clip the contour
    
    ls[[as.character(pop)]] <- contourLines(dens, levels=cont_level)
  }
  #names(ls) <- names(reduction_df)
  return(ls)
}

#arrange ls as a dataframe for plotting
getContourCoordinates <- function(ls){
  cont_df_list <- list()
  for(i in seq_along(ls)){
    x1 <- unlist(sapply(ls[[i]], '[[', 2))
    y1 <- unlist(sapply(ls[[i]], '[[', 3))
    length_contours <- sapply(ls[[i]], function(x) length(x[[2]]))
    cont_name <- unlist(mapply(rep, times=length_contours, x=seq_along(ls[[i]])))
    cont_df_list[[i]] <- data.frame(x=x1, y=y1, cont_name, grouping_var = names(ls)[i])
  }
  cont_df <- do.call(rbind, cont_df_list)
  cont_df <- cont_df %>% unite('cont_group', c('grouping_var', 'cont_name'), remove = F)
  return(cont_df)
}
 

# get_density <- function(x, y, ...) {
#   dens <- MASS::kde2d(x, y, ...)
#   ix <- findInterval(x, dens$x)
#   iy <- findInterval(y, dens$y)
#   ii <- cbind(ix, iy)
#   return(dens$z[ii])
# }
# 
# #required to estimate density separately by sample to obtain correct estimate for each subplot
# get_density_by_groups <- function(df, group, x, y, z, ...) {
#   lDf <- split(df, df[[group]])
#   lDf <- lapply(lDf, function(mydf){mydf[[z]] <- get_density(mydf[[x]], mydf[[y]], ...); return(mydf)})
#   lDf <- do.call(rbind, lDf)
#   return(lDf)
# }
