#' plot features
#'
#' Extension of the functionalities provided by Seurats FeaturePlot/DimPlot functions
#' Can plot scatterplots of multiple features provided under `features` (from Seurat assay) or `meta_features` (from Seurat meta.data or data.frame)
#' Embeddings from Seurat need to be specified under `reduction`. Alternatively, a data.frame with the embedding coordinates in the first two columns 
#' can be passed to `reduction` and seurat_obj can be set to NULL.
#' Loads ggplot2 and dplyr libraries.
#' 
#' @param seurat_obj optional Seurat object
#' @param reduction either a data.frame with embeddings in column 1+2 and further meta data or a string specifying a dimensionality reduction from the Seurat object
#' @param features If Seurat object is provided, extract features from the specified assay and slot
#' @param meta_features if Seurat object is provided, extract features from the meta.data. Otherwise extract features from the data.frame provided in `reduction`
#' @param assay assay to use from the Seurat object, default: DefaultAssay(seurat_obj)
#' @param slot slot to use from the Seurat object, default: data
#' @param colorScale color scale to use for continuous features, one of 'wbr', 'bry', 'seurat', 'viridis'
#' @param cap_top upper value threshold passed to jj_cap_vals
#' @param cap_bottom lower value threshold passed to jj_cap_vals
#' @param custom_colors named vector of colours used to colour categorical features
#' @param custom_theme custom theme to be used, default: theme_minimal()
#' @param shape point shape, default: 16
#' @param alpha alpha value, default: 1
#' @param pt.size size of points, default: 0.1
#' @param return_gg_object return ggplot object instead of plotting, default: FALSE
#' @param my_title optional title for the ggplot
#' @param no_legend omit legend, default: FALSE
#' @param facet_by string specifying a meta_feature to facet by
#' @param n_facet_rows number of rows for facetted plots
#' @param cont_or_disc string of length 1 or length n features indicating whether the features are continuous 'c' or discrete 'd'. Try to set this manually, when the function fails. Otherwise, set to 'a' to automatically determine c or d for each feature.
#' @param use_pointdensity colour by pointdensity using the ggpointdensity package. 
#' @param pointdensity_subset Only used if use_pointdensity=T. If NULL, use all cells. If set to groups within meta_features, only calculate density for these subgroups
#' @param facet_subset Only used when facet_by not FALSE. Only plot the facets for the groups supplied here. Background cells are still shown for the whole dataset 
#' @param order order points so that largest values are on top
#' @param background_cells when using facets, include the cells not part of the facet as grey background 
#' @param label add boxes with labels to the discrete variable
#' @param box_col colour to fill boxes, if label = T. If NULL, use colours from the respective groups
#' @param xlabel x axis label, default: UMAP 1
#' @param ylabel y axis label, default: UMAP 2
#' @keywords plot
#' @export
#' @examples
#' df = data.frame(umap1=c(1,5,3), umap2=c(2,5,2), fruit=c('Apple', 'Banana', 'Apple'), dish=c('Salad','Salad','Snack'))
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, shape=15, my_title='fruits')
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, facet_by='fruit', n_facet_rows=2)
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, facet_by='dish')
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, custom_colors=c(Apple='green', Banana='yellow'), label=T)
#' jj_plot_features(seurat_rna, features=c('CD4', 'CD8A'), cap_top='q95', colorScale='viridis')
#' df2 = data.frame(a=rnorm(100, 0, 5), b=rnorm(100, 0, 5), d=rbinom(100, 50, 0.3), e = sample(c('A','B', 'C'), 100, replace=T))
#' jj_plot_features(reduction=df2, meta_features=c('d'), pt.size=4, facet_by = 'e', background_cells = T, order=T, custom_theme = theme_bw())
#' jj_plot_features(reduction=df2, meta_features='e', pt.size=4, facet_by = 'e', use_pointdensity = T, order=T, custom_theme = theme_bw())
#' jj_plot_features(reduction=df2, meta_features='e', pt.size=4, use_pointdensity = T, pointdensity_subset = c('B','C'), order=T, custom_theme = theme_bw())
#' jj_plot_features(reduction=df2, meta_features='e', pt.size=4, use_pointdensity = T, pointdensity_subset = c('B','C'), background_cells=T, order=T, custom_theme = theme_bw())
#'

jj_plot_features <- function(seurat_obj=NULL, reduction=NULL, features=NULL, meta_features=NULL,
                             assay=NULL, slot='data', 
                             colorScale=c('viridis', 'wbr', 'gbr', 'bry', 'seurat'),
                             facet_by=NULL, cap_top=NULL,  cap_bottom=NULL, 
                             custom_colors=NULL, custom_theme=theme_minimal(), shape = 16, alpha=1,
                             pt.size=0.5, return_gg_object=FALSE, my_title=NULL, 
                             no_legend=F, n_facet_rows=NULL, facet_subset=NULL,
                             cont_or_disc = 'a', use_pointdensity = FALSE,
                             pointdensity_subset=NULL, order=FALSE, 
                             background_cells=FALSE, label=FALSE, box_col=NULL, convert_factors=FALSE,
                             xlabel = 'UMAP 1', ylabel = 'UMAP 2'){
  
  if(is.null(reduction)){
    stop('reduction must be either string specifying the reduction to use from seurat object or a dr_df data.frame')
  }
  if(!is.null(seurat_obj)){
    stopifnot(class(seurat_obj)=='Seurat')
  }else{
    if(!is.null(features)){
      stop('`features` argument can only be used when seurat object is provided. Use `meta_features` instead.' )
    }
  }
  
  library(ggplot2)
  library(dplyr)
  
  gg_list <- list()
  message('getting reduction coordinates')
  colorScale <- match.arg(colorScale)
  if(is.data.frame(reduction) | is.matrix(reduction)){
    dr_df <- as.data.frame(reduction)
  }else{
    dr_df <- jj_get_reduction_coords(seurat_obj, reduction, exact_match = TRUE)
  }
  stopifnot(all(apply(dr_df[, c(1,2)], 2, is.numeric))) #first two columns should be pc/umap coordinates
  colnames(dr_df)[1:2] <- c('dim_1', 'dim_2')
  
  if(is.null(features) & is.null(meta_features)){
    meta_features = 'data'
    dr_df$data = ''
    goi = 'data'
  }
  
  cont_or_disc <- unlist(strsplit(cont_or_disc, split = ''))
  if(length(cont_or_disc) == 1){
    cont_or_disc <- rep(cont_or_disc, length(c(features, meta_features)))
  }else if(length(cont_or_disc) != length(c(features, meta_features))){
    stop('length of cont_or_disc must be either 1 or length(features, meta_features)')
  }else if(any(!cont_or_disc %in% c('c', 'd', 'a'))){
    stop('cont_or_disc can only be c, d, a: continuous, discrete, automatic')
  }
  
  if(!is.null(meta_features)){
    message('Getting features from metadata')
    goi <- meta_features
    if(!all(goi %in% colnames(dr_df))){
      warning(sprintf('%s not found in the assay meta.data.',
                      paste(goi[!goi %in% colnames(dr_df)], collapse = ', ')))
      goi = goi[goi %in% colnames(dr_df)]
      if(length(goi)==0){return(NULL)}
    }
    
    #convert_factors to numeric if possible and otherwise to character
    if(convert_factors){
      for(me in goi){
        if(is.factor(dr_df[, me])){
          dr_df[, me] = as.character(dr_df[, me])
          if(Hmisc::all.is.numeric(dr_df[, me])){
            dr_df[, me] <- as.numeric(dr_df[, me])
          }
        }
      }
    }

    if(!is.null(cap_top) | !is.null(cap_bottom)){
      message('Capping meta feature values')
      for(i in goi){
        if(is.numeric(dr_df[, i])){
          dr_df[, i] <- jj_cap_vals(dr_df[, i], cap_top=cap_top,  cap_bottom=cap_bottom)
        }
      }
    }
    #try to find continous/discrete status for each variable automatically
    if('a' %in% cont_or_disc){
      cont_or_disc = vector()
      #cont_or_disc = paste(ifelse(sapply(dr_df[, colnames(dr_df) %in% goi], class) == 'character', 'd', 'c'), collapse = '')
      for(i in goi){
        if(class(dr_df[, colnames(dr_df) %in% i]) %in% c('character','factor')){
          cont_or_disc = c(cont_or_disc, 'd')
        }else{
          if(n_distinct(dr_df[, colnames(dr_df) %in% i]) < 11){
            cont_or_disc = c(cont_or_disc, 'd')
          }else{
            cont_or_disc = c(cont_or_disc, 'c')
          }
        }
      }
      #cont_or_disc = paste(cont_or_disc, collapse = '')
      message(cont_or_disc)
    }
  }
  
  if(!is.null(features)){
    if(is.null(assay)){
      assay = DefaultAssay(seurat_obj)
    }
    message(sprintf('getting %s slot from %s assay of seurat object', slot, assay))
    if(!all(features %in% colnames(dr_df))){
      dr_df <- jj_bind_features_with_dr_df(seurat_obj, assay=assay, slot=slot, 
                                           features=features, dr_df=dr_df, cap_top=cap_top, 
                                           cap_bottom=cap_bottom, log10Transform=FALSE)
    }
    goi <- gsub('-', '_', features)
    #colnames(dr_df)[colnames(dr_df) %in% features] <- goi
    #replace :: with __ to avoid error in ggplot (naming from chromvar tfs)
    colnames(dr_df) <-  gsub('::', '___', colnames(dr_df))
    goi <-  gsub('::', '___', goi)
    colnames(dr_df) <-  gsub(':', '__', colnames(dr_df))
    goi <-  gsub(':', '__', goi)
    colnames(dr_df) <-  gsub('?', '', colnames(dr_df), fixed = T)
    goi <-  gsub('?', '', goi, fixed = T)
    
    if(!is.null(meta_features)) goi <- c(goi, meta_features)
  }
  #dr_df_molten <- reshape2::melt(dr_df, id=c('UMAP_1', 'UMAP_2'))
  range_x <- range(dr_df$dim_1)
  range_y <- range(dr_df$dim_2)
  
  message('plotting features')
  for(i in seq_along(goi)){
    message(sprintf('%i/%i: %s', i, length(goi), goi[i]))
    
    # if(is.factor(dr_df[, goi[i]])){
    #   dr_df[, goi[i]] <- as.character(dr_df[, goi[i]])
    #   if(Hmisc::all.is.numeric(dr_df[, goi[i]]) & n_distinct(dr_df[, goi[i]]) >=30){
    #     message("Converting factor to numeric")
    #     dr_df[, goi[i]] <- as.numeric(dr_df[, goi[i]])
    #   }else{
    #     message("Converting factor to character.")
    #   }
    # }
    
    if(cont_or_disc[i] == 'd'){
      #dr_df[, goi[i] ] = as.character(dr_df[, goi[i] ])
      if(!is.null(custom_colors)){
        if(is.null(names(custom_colors))){
          names(custom_colors) = levels(as.factor(dr_df[, goi[i] ]))
        }
        breaks_use = names(custom_colors)
        dr_df[, goi[i] ] <- factor(dr_df[, goi[i] ], levels= breaks_use)
      #}else if(suppressWarnings(sum(is.na(as.numeric(names(jj_get_jj_colours(levels(as.factor(dr_df[, goi[i] ]))))))) == 1 )){
      #  breaks_use <- suppressWarnings(as.character(sort(as.numeric(names(jj_get_jj_colours(dr_df[, goi[i]]))), decreasing = F)))
      #  dr_df[, goi[i] ] <- factor(dr_df[, goi[i] ], levels= breaks_use)
      }else{
        breaks_use <- names(jj_get_jj_colours(levels(as.factor(dr_df[, goi[i]]))))
      }
    }
    
    if(order){
      dr_df = dr_df[order(dr_df[, goi[i]]), ]
    }
    

    if(!is.null(facet_by) & background_cells){
      stopifnot(facet_by %in% colnames(dr_df))
      fac_levels = unique(dr_df[, facet_by])
      if(!is.null(facet_subset[1])){
        fac_levels = fac_levels[fac_levels %in% facet_subset]
      }
      dr_df_background_list = list()
      for(ii in seq_along(fac_levels)){
        dr_df_background_list[[ii]] = dr_df
        dr_df_background_list[[ii]][, facet_by] = fac_levels[ii]
      }
      dr_df_background = do.call(rbind, dr_df_background_list)
      rm(dr_df_background_list)
      gg <- ggplot() + 
        geom_point(data = dr_df_background, aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha, color='grey80')
    }else if(background_cells){
      gg <- ggplot() + 
        geom_point(data = dr_df, aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha, color='grey80')
    }else{
      gg = NULL
    }
    
    if(!is.null(facet_by) & !is.null(facet_subset[1])){
      dr_df = dr_df[dr_df[, facet_by] %in% facet_subset, ]
      if(!nrow(dr_df) > 0){
        stop('Groups in `facet_subset` are not available in `facet_by`')
      }
    }
    
    if(use_pointdensity){
      if(!is.null(gg)){
        if(any(pointdensity_subset %in% dr_df[, goi[i] ])){
          gg = gg + 
            ggpointdensity::geom_pointdensity(data = dr_df[dr_df[, goi[i]] %in% pointdensity_subset, ], mapping = aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }else{
          gg = gg + 
            ggpointdensity::geom_pointdensity(data = dr_df, mapping = aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }
      }else{
        if(any(pointdensity_subset %in% dr_df[, goi[i] ])){
          gg = ggplot() + 
            ggpointdensity::geom_pointdensity(data =  dr_df[dr_df[, goi[i]] %in% pointdensity_subset, ], mapping = aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }else{
          gg = ggplot() + 
            ggpointdensity::geom_pointdensity(data = dr_df, mapping = aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }
      }
    }else{
      if(is.character(shape)){
        if(!is.null(gg)){
          gg <- gg +  
            geom_point(data=dr_df, aes_string(x='dim_1', y='dim_2', colour=goi[i], shape=shape), size=pt.size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }else{
          gg = ggplot(dr_df, aes(x=dim_1, y=dim_2)) + 
            geom_point(aes_string(colour=goi[i], shape=shape), size=pt.size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }
      }else{
        if(!is.null(gg)){
          gg <- gg + 
            geom_point(data=dr_df, aes_string(x='dim_1', y='dim_2', colour=goi[i]), size=pt.size, shape=shape, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }else{
          gg <- ggplot(dr_df, aes(x=dim_1, y=dim_2)) + 
            geom_point(aes_string(colour=goi[i]), size=pt.size, shape=shape, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }
      }
    }
    
    if(use_pointdensity){
      gg = gg + viridis::scale_color_viridis()
    }else if(!all(is.null(custom_colors))){ #& !n_distinct(dr_df[, goi[i]]) < 30){
      if(length(custom_colors)==3 & cont_or_disc[i] == 'c'){
        #use gradient if 3 colours are supplied
        mean_acc <- (max(dr_df[, goi[i]], na.rm = T) + min(dr_df[, goi[i]], na.rm = T)) / 2 #mean(dr_df[, goi[i]])
        gg <- gg + scale_color_gradient2(low = custom_colors[1], mid = custom_colors[2], high = custom_colors[3], midpoint = mean_acc)
      }else if(!is.null(names(custom_colors))){
        #use discrete color mapping if named vector is supplied
        #custom_colors <- custom_colors[names(custom_colors) %in% dr_df[, goi[i]]]
        message('Using colors provided unter custom_colors.')
        gg <- gg +
          scale_colour_manual(values = custom_colors) +#, breaks = names(custom_colors)) + 
          guides(colour = guide_legend(override.aes = list(size=5)))
      }
    }else if(cont_or_disc[i] == 'd'){
      #TODO: can error if discrete number of values is > =30
      #plot with discrete colors if less than 30 distinct values in variable (or if d is set (check for <60 distinct groups as sanity check))
      gg <- gg +
        scale_colour_manual(values = jj_get_jj_colours(breaks_use)) + 
        guides(colour = guide_legend(override.aes = list(size=5)))
    }else if(colorScale=='viridis'){
      gg <- gg + viridis::scale_color_viridis()
    }else if(colorScale=='seurat'){
      gg <- gg + scale_color_gradient(low='#C3C3C3', high = '#1C0DFD')
    }else if(colorScale=='bry'){
      mean_acc <- (max(dr_df[, goi[i]],na.rm = T)+ min(dr_df[, goi[i]],na.rm = T)) / 2 #mean(dr_df[, goi[i]])
      gg <- gg +  scale_color_gradient2(low = "darkblue", mid = "red", high = "yellow", midpoint = mean_acc)
    }else if(colorScale=='wbr'){
      mean_acc <- (max(dr_df[, goi[i]], na.rm = T) + min(dr_df[, goi[i]], na.rm = T)) / 2 #mean(dr_df[, goi[i]])
      print(mean_acc)
      gg <- gg + scale_color_gradient2(low = "#ffffd9", mid = "blue", high = "red", midpoint = mean_acc)
    }else if(colorScale=='gbr'){
      mean_acc <- (max(dr_df[, goi[i]], na.rm = T) + min(dr_df[, goi[i]], na.rm = T)) / 2 #mean(dr_df[, goi[i]])
      print(mean_acc)
      gg <- gg + scale_color_gradient2(low = "grey80", mid = "blue", high = "red", midpoint = mean_acc)
    }
    
    gg = gg + xlab(xlabel) + ylab(ylabel)
    
    if(is.theme(custom_theme)){
      gg <- gg + custom_theme
    }
    if(no_legend){
      gg = gg + theme(legend.position='none')
    }
    if(!is.null(my_title)){
      if(length(my_title) == length(goi)){
        gg <- gg + ggtitle(my_title[i]) +  theme(plot.title = element_text(hjust = 0.5))
      }
      else if(length(my_title) == 1){
        gg <- gg + ggtitle(my_title) +  theme(plot.title = element_text(hjust = 0.5))
      }
    }
    
    if(!is.null(facet_by)){
      if(facet_by %in% colnames(dr_df) & n_distinct(dr_df[, facet_by]) < 100){
        message(sprintf('Facetting by %s.', facet_by))
        if(is.null(n_facet_rows)){
          n_facet_rows <- round(sqrt(n_distinct(dr_df[, facet_by])))
        }
        gg <- gg + facet_wrap(as.formula(sprintf('.~%s', facet_by)), nrow=n_facet_rows)
      }
    }
    # if(facet_by %in% c(TRUE, 1)){
    #   if(facet_by & n_distinct(dr_df[, goi[i]]) < 100){
    #     message(sprintf('Facetting by %s.', goi[i]))
    #     if(is.null(n_facet_rows)){
    #       n_facet_rows <- round(sqrt(n_distinct(dr_df[, goi[i]])))
    #     }
    #     if(!is.logical(facet_by)) print(gg) #print unfacetted + facetted if facet_by <- 1 is used
    #     gg <- gg + facet_wrap(as.formula(sprintf('.~%s', goi[i])), nrow=n_facet_rows)
    #   }
    # }
    
    if(label & cont_or_disc[i] == 'd'){
      gg$data[, goi[i]] = as.factor(gg$data[, goi[i]])
      gg = .LabelClusters(gg, goi[i], box = T, col_use = box_col)
    }
    
    if(return_gg_object){
      gg_list[[i]] <- gg
    }else{
      print(gg)
    }
  }
  if(return_gg_object){
    return(gg_list)
  }
}

#' @export
.LabelClusters = function(
  #function from Seurat to label clusters
  #additional:
  #1 check is done: if id is not factor, stop. In original function this results in another
  #error message difficult to understand
  #2 convert factor to character, otherwise custom labels results in error due to invalid factor levels
  plot,
  id,
  clusters = NULL,
  labels = NULL,
  split.by = NULL,
  repel = TRUE,
  box = FALSE,
  geom = 'GeomPoint',
  position = "median",
  col_use = NULL,
  ...
) {
  #dr_df <- get_reduction_coords(seurat_atac, 'umap')
  #dr_df$singler_label = as.factor(dr_df$singler_label)
  #gg = ggplot(dr_df[1:10, ], aes(x=dim_1, y=dim_2)) + geom_point(aes(colour=singler_label))
  #LabelClusters(gg, id='singler_label')
  if(repel) library(ggrepel)
  `%||%` <- function(lhs, rhs) {
    if (!is.null(x = lhs)) {
      return(lhs)
    } else {
      return(rhs)
    }
  }
  
  GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
    geoms <- sapply(
      X = plot$layers,
      FUN = function(layer) {
        return(class(x = layer$geom)[1])
      }
    )
    geoms <- which(x = geoms == geom)
    if (length(x = geoms) == 0) {
      stop("Cannot find a geom of class ", geom)
    }
    geoms <- min(geoms)
    if (plot.first) {
      x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
      y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
    } else {
      x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
      y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
    }
    return(list('x' = x, 'y' = y))
  }
  
  xynames <- unlist(x = GetXYAesthetics(plot = plot, geom = geom), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  
  #own condition:
  if(!is.factor(data[, id])){
    stop('id needs to be a factor')
  }
  
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% possible.clusters], collapse = ","))
  }
  pb <- ggplot_build(plot = plot)
  if (geom == 'GeomSpatial') {
    data[, xynames["y"]] = max(data[, xynames["y"]]) - data[, xynames["y"]] + min(data[, xynames["y"]])
    if (!pb$plot$plot_env$crop) {
      y.transform <- c(0, nrow(x = pb$plot$plot_env$image)) - pb$layout$panel_params[[1]]$y.range
      data[, xynames["y"]] <- data[, xynames["y"]] + sum(y.transform)
    }
  }
  data <- cbind(data, color = pb$data[[1]][[1]])
  labels.loc <- lapply(
    X = groups,
    FUN = function(group) {
      data.use <- data[data[, id] == group, , drop = FALSE]
      data.medians <- if (!is.null(x = split.by)) {
        do.call(
          what = 'rbind',
          args = lapply(
            X = unique(x = data.use[, split.by]),
            FUN = function(split) {
              medians <- apply(
                X = data.use[data.use[, split.by] == split, xynames, drop = FALSE],
                MARGIN = 2,
                FUN = median,
                na.rm = TRUE
              )
              medians <- as.data.frame(x = t(x = medians))
              medians[, split.by] <- split
              return(medians)
            }
          )
        )
      } else {
        as.data.frame(x = t(x = apply(
          X = data.use[, xynames, drop = FALSE],
          MARGIN = 2,
          FUN = median,
          na.rm = TRUE
        )))
      }
      data.medians[, id] <- group
      data.medians$color <- data.use$color[1]
      return(data.medians)
    }
  )
  if (position == "nearest") {
    labels.loc <- lapply(X = labels.loc, FUN = function(x) {
      group.data <- data[as.character(x = data[, id]) == as.character(x[3]), ]
      nearest.point <- nn2(data = group.data[, 1:2], query = as.matrix(x = x[c(1,2)]), k = 1)$nn.idx
      x[1:2] <- group.data[nearest.point, 1:2]
      return(x)
    })
  }
  labels.loc <- do.call(what = 'rbind', args = labels.loc)
  labels.loc[, id] <- factor(x = labels.loc[, id], levels = levels(data[, id]))
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels),  ") must be equal to the number of clusters being labeled (", length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  #convert factor to character, otherwise custom labels results in error due to invalid factor levels
  labels.loc[, id] = as.character(labels.loc[, id])
  
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  if (box) {
    geom.use <- ifelse(test = repel, yes = geom_label_repel, no = geom_label)
    if(!is.null(col_use)){
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id), fill = col_use,
        show.legend = FALSE,
        ...
      )
    }else{
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id, fill = id),
        show.legend = FALSE,
        ...
      ) + scale_fill_manual(values = labels.loc$color[order(labels.loc[, id])])
    }
  } else {
    geom.use <- ifelse(test = repel, yes = geom_text_repel, no = geom_text)
    if(!is.null(col_use)){
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id), colour=col_use,
        show.legend = FALSE,
        ...
      )
    }else{
      plot <- plot + geom.use(
        data = labels.loc,
        mapping = aes_string(x = xynames['x'], y = xynames['y'], label = id, colour=id),
        show.legend = FALSE,
        ...
      )
    }
  }
  return(plot)
}

#TODO replace repetitions of code with this helper function
cget_scale_midpoint = function(feat){
  (max(feat, na.rm = T) + min(feat, na.rm = T)) / 2 
}