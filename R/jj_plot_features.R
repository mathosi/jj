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
#' @param assay assay to use from the Seurat object
#' @param slot slot to use from the Seurat object
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
#' @param use_facets either TRUE/FALSE: facet by the meta_feature, or a string specifying another meta_feature to facet by
#' @param n_facet_rows number of rows for facetted plots
#' @param cont_or_disc string indicating whether the the meta_feature is continuous 'c' or discrete 'd'. Try to set this manually, when the function fails. Otherwise, set to 'u'
#' @param pointdensity colour by pointdensity instead of feature using the ggpointdensity package
#' @param order order points so that largest values are on top
#' @param background_cells when using facets, include the cells not part of the facet as grey background 
#' @param cont_or_disc_auto take the class from the df column to decide, and convert numeric to character if n_distinct < 11 visualize colours as ggplot heatmap, default: F
#' @keywords plot
#' @export
#' @examples
#' df = data.frame(umap1=c(1,5,3), umap2=c(2,5,2), fruit=c('Apple', 'Banana', 'Apple'), dish=c('Salad','Salad','Snack'))
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, shape=15, my_title='fruits')
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, use_facets=T, n_facet_rows=2)
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, use_facets='dish')
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, custom_colors=c(Apple='green', Banana='yellow'))
#' jj_plot_features(seurat_rna, features=c('CD4', 'CD8A'), cap_top='q95', colorScale='viridis')

jj_plot_features <- function(seurat_obj=NULL, reduction=NULL, features=NULL, meta_features=NULL,
                         assay='RNA', slot='counts', 
                         colorScale=c('wbr', 'bry', 'seurat', 'viridis'),
                         use_facets=FALSE, cap_top=NULL,  cap_bottom=NULL, 
                         custom_colors=NULL, custom_theme=theme_minimal(), shape = 16, alpha=1,
                         pt.size=0.1, return_gg_object=FALSE, my_title=NULL, 
                         no_legend=F, n_facet_rows=NULL,
                         cont_or_disc = 'u', cont_or_disc_auto = FALSE,
                         pointdensity=F, order=FALSE, background_cells=FALSE){

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
  
  cont_or_disc <- unlist(strsplit(cont_or_disc, split = ''))
  if(length(cont_or_disc) == 1){
    cont_or_disc <- rep(cont_or_disc, length(c(features, meta_features)))
  }else if(length(cont_or_disc) != length(c(features, meta_features))){
    stop('length of cont_or_disc must be either 1 or length(features, meta_features)')
  }else if(any(!cont_or_disc %in% c('c', 'd', 'u'))){
    stop('cont_or_disc can only be c, d, u: continuous, discrete, unknown')
  }
  
  gg_list <- list()
  message('getting reduction coordinates')
  colorScale <- match.arg(colorScale)
  if(is.data.frame(reduction)){
    dr_df <- reduction
  }else{
    dr_df <- jj_get_reduction_coords(seurat_obj, reduction, exact_match = TRUE)
  }
  stopifnot(all(apply(dr_df[, c(1,2)], 2, is.numeric))) #first two columns should be pc/umap coordinates
  colnames(dr_df)[1:2] <- c('dim_1', 'dim_2')
  
  if(!is.null(meta_features)){
    message('Getting features from metadata')
    goi <- meta_features
    if(!all(goi %in% colnames(dr_df))){
      warning(sprintf('%s not found in the assay meta.data.',
                      paste(goi[!goi %in% colnames(dr_df)], collapse = ', ')))
      goi = goi[goi %in% colnames(dr_df)]
      if(length(goi)==0){return(NULL)}
    }
    
    for(me in goi){
      if(is.factor(dr_df[, me])){
        dr_df[, me] <- as.numeric(as.character(dr_df[, me]))
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
    if(cont_or_disc_auto){
      cont_or_disc = vector()
      #cont_or_disc = paste(ifelse(sapply(dr_df[, colnames(dr_df) %in% goi], class) == 'character', 'd', 'c'), collapse = '')
      for(i in goi){
        if(class(dr_df[, colnames(dr_df) %in% i]) == 'character'){
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
      print(cont_or_disc)
    }
  }
  
  if(!is.null(features)){
    message('getting feature matrix')
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
    
    if(is.factor(dr_df[, goi[i]])){
      dr_df[, goi[i]] <- as.character(dr_df[, goi[i]])
      if(Hmisc::all.is.numeric(dr_df[, goi[i]]) & n_distinct(dr_df[, goi[i]]) >=30){
        message("Converting factor to numeric")
        dr_df[, goi[i]] <- as.numeric(dr_df[, goi[i]])
      }else{
        message("Converting factor to character.")
      }
    }
    
    if((n_distinct(dr_df[, goi[i]]) < 30 & cont_or_disc[i] == 'u') | cont_or_disc[i] == 'd'){
      if(!is.null(custom_colors)){
        if(is.null(names(custom_colors))){
          names(custom_colors) = levels(as.factor(dr_df[, goi[i] ]))
        }
        breaks_use = names(custom_colors)
        dr_df[, goi[i] ] <- factor(dr_df[, goi[i] ], levels= breaks_use)
      }else if(suppressWarnings(sum(is.na(as.numeric(names(jj_get_jj_colours(dr_df[, goi[i] ]))))) == 1 )){
        breaks_use <- suppressWarnings(as.character(sort(as.numeric(names(jj_get_jj_colours(dr_df[, goi[i]]))), decreasing = F)))
        dr_df[, goi[i] ] <- factor(dr_df[, goi[i] ], levels= breaks_use)
      }else{
        breaks_use <- names(jj_get_jj_colours(dr_df[, goi[i]]))
      }
    }
    
    if(order){
      dr_df = dr_df[order(dr_df[, goi[i]]), ]
    }
    
    if(use_facets %in% colnames(dr_df) & background_cells){
      fac_levels = unique(dr_df[, use_facets])
      dr_df_background_list = list()
      for(ii in seq_along(fac_levels)){
        dr_df_background_list[[ii]] = dr_df
        dr_df_background_list[[ii]][, use_facets] = fac_levels[ii]
      }
      dr_df_background = do.call(rbind, dr_df_background_list)
      rm(dr_df_background_list)
      gg <- ggplot() + 
        geom_point(data = dr_df_background, aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha, color='grey80')
    }else{
      gg = NULL
    }
    
    if(pointdensity){
      if(!is.null(gg)){
        gg = gg + 
          ggpointdensity::geom_pointdensity(data = dr_df, mapping = aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha) +
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
    
    if(pointdensity){
      gg = gg + viridis::scale_color_viridis()
    }else if(!all(is.null(custom_colors))){ #& !n_distinct(dr_df[, goi[i]]) < 30){
      if(length(custom_colors)==3 & cont_or_disc[i] %in% c('c', 'u')){
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
    }else if(n_distinct(dr_df[, goi[i]]) < 30 & cont_or_disc[i] == 'u' | cont_or_disc[i] == 'd'){
      #TODO: can error if discrete number of values is > =30
      #plot with discrete colors if less than 30 distinct values in variable (or if d is set (check for <60 distinct groups as sanity check))
      gg <- gg +
        scale_colour_manual(values = jj_get_jj_colours(dr_df[, goi[i]]), breaks = breaks_use) + 
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
    }
    
    if(is.theme(custom_theme)){
      gg <- gg + custom_theme
    }
    if(no_legend){
      gg = gg + theme(legend.position='none')
    }
    if(!is.null(my_title)){
      if(length(my_title) == length(goi)){
        gg <- gg + ggtitle(my_title[i])
      }
      else if(length(my_title) == 1){
        gg <- gg + ggtitle(my_title)
      }
    }
    
    if(!is.logical(use_facets)){
      if(use_facets %in% colnames(dr_df) & n_distinct(dr_df[, use_facets]) < 30){
        message(sprintf('Facetting by %s.', use_facets))
        if(is.null(n_facet_rows)){
          n_facet_rows <- round(sqrt(n_distinct(dr_df[, use_facets])))
        }
        gg <- gg + facet_wrap(as.formula(sprintf('.~%s', use_facets)), nrow=n_facet_rows)
      }
    }
    if(use_facets %in% c(TRUE, 1)){
      if(use_facets & n_distinct(dr_df[, goi[i]]) < 100){
        message(sprintf('Facetting by %s.', goi[i]))
        if(is.null(n_facet_rows)){
          n_facet_rows <- round(sqrt(n_distinct(dr_df[, goi[i]])))
        }
        if(!is.logical(use_facets)) print(gg) #print unfacetted + facetted if use_facets <- 1 is used
        gg <- gg + facet_wrap(as.formula(sprintf('.~%s', goi[i])), nrow=n_facet_rows)
      }
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
