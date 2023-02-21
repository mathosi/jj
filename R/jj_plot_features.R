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
#' @param use_no_legend omit legend, default: FALSE
#' @param facet_by string specifying a meta_feature to facet by
#' @param n_facet_rows number of rows for facetted plots
#' @param cont_or_disc string of length 1 or length n features indicating whether the features are continuous 'c' or discrete 'd'. Try to set this manually, when the function fails. Otherwise, set to 'a' to automatically determine c or d for each feature.
#' @param use_pointdensity colour by pointdensity using the ggpointdensity package. 
#' @param show_background_cells bool, works if foreground_subset_bool is specified. Also when using facets, include the cells not part of the facet as grey background 
#' @param foreground_subset_bool boolean vector with same length as number of cells. If show_background_cells=TRUE, show all cells with FALSE as grey background. Do also exlcude those cells for pointdensity calculation, if use_pointdensity=TRUE. Currently not working with do_label=T
#' @param facet_subset Only used when facet_by not NULL Only plot the facets for the groups supplied here. Background cells are still shown for the whole dataset 
#' @param do_order bool, order points so that largest values are on top
#' @param do_shuffle bool, randomly shuffle the plotting order. Supersedes do_order.
#' @param label_type one of geom_text, geom_text_repel, geom_label, geom_label_repel
#' @param label_col colour to fill boxes, if do_label = T. If NULL, use colours from the respective groups
#' @param xlabel x axis label, default: UMAP 1
#' @param ylabel y axis label, default: UMAP 2
#' @keywords plot
#' @export
#' @examples
#' df = data.frame(umap1=c(1,5,3), umap2=c(2,5,2), fruit=c('Apple', 'Banana', 'Apple'), dish=c('Salad','Salad','Snack'))
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, shape=15, my_title='fruits')
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, facet_by='fruit', n_facet_rows=2)
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, facet_by='dish')
#' jj_plot_features(reduction=df, meta_features=c('fruit'), pt.size=4, custom_colors=c(Apple='green', Banana='yellow'))
#' jj_plot_features(pbmc_small, reduction='tsne', features=c('MS4A1', 'CD79A'), cap_top='q95', colorScale='viridis', pt.size = 2, xlabel='tsne_1', ylabel='tsne_2',my_title=c('tnse: MS4A1','tsne: CD79A'))
#' df2 = data.frame(a=rnorm(100, 0, 5), b=rnorm(100, 0, 5), d=rbinom(100, 50, 0.3), e = sample(c('A','B', 'C'), 100, replace=T))
#' jj_plot_features(reduction=df2, meta_features=c('d'), pt.size=4, facet_by = 'e', show_background_cells = T, do_order=T, custom_theme = theme_bw())
#' jj_plot_features(reduction=df2, meta_features='e', pt.size=4, facet_by = 'e', use_pointdensity = T, do_order=T, custom_theme = theme_bw())
#' jj_plot_features(reduction=df2, meta_features='e', pt.size=4, use_pointdensity = T, foreground_subset_bool = df2$e %in% c('B','C'),show_background_cells=T, do_order=T, custom_theme = theme_bw())
#' jj_plot_features(reduction=df2, meta_features='e', pt.size=4, use_pointdensity = T, foreground_subset_bool = df2$e %in% c('B','C'), show_background_cells=T, facet_by = 'e', facet_subset = c('A','B'), custom_theme = theme_bw())
#' data_df = data.frame(x=c(rnorm(200, 5), rnorm(70,2), rnorm(50,2.5, 2)), y=c(rnorm(200, 5), rnorm(70,-2), rnorm(50,2.5,2)), clusterIdent = c(rep('gr1', 200), rep('gr2', 70), rep('gr3', 50)), groupVar = c(rep('AA', 180), rep('BB', 140)))
#' mycols = structure(c('red','green', 'pink'), names = c('gr1','gr2','gr3'))
#' jj_plot_features(reduction = data_df, meta_features = 'clusterIdent', return_gg_object = T, pt.size = 2,
#'                  facet_by = 'groupVar', label_type = 'geom_label_repel', label_subset = c('gr1','gr2'), custom_colors=mycols)[[1]]

jj_plot_features <- function(seurat_obj=NULL, reduction=NULL, features=NULL, meta_features=NULL,
                             assay=NULL, slot='data', 
                             colorScale=c('viridis', 'wbr', 'gbr', 'bry', 'seurat'),
                             facet_by=NULL, cap_top=NULL,  cap_bottom=NULL, 
                             custom_colors=NULL, custom_theme=theme_minimal(), shape = 16, alpha=1,
                             pt.size=0.5, return_gg_object=FALSE, my_title=NULL,
                             use_no_legend=F, n_facet_rows=NULL, facet_subset=NULL, foreground_subset_bool = NULL,
                             cont_or_disc = 'a', use_pointdensity = FALSE,
                             #pointdensity_subset=NULL, 
                             do_order=FALSE, do_shuffle=FALSE, 
                             show_background_cells=FALSE, label_type=NULL, label_subset = NULL, label_col=NULL, convert_factors=FALSE,
                             xlabel = 'UMAP 1', ylabel = 'UMAP 2'){
  #@param pointdensity_subset Only used if use_pointdensity=T. If NULL, use all cells. If set to groups within meta_features, only calculate density for these subgroups
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
  
  if(!is.null(foreground_subset_bool)[1]){
    if(!show_background_cells) warning('`foreground_subset_bool` will have no effect if `show_background_cells` = FALSE')
    dr_df$foreground_subset_bool = foreground_subset_bool
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
    
    if(do_order){
      dr_df = dr_df[order(dr_df[, goi[i]]), ]
    }
    if(do_shuffle){
      if(do_order) warning('Both `do_order` and `do_shuffle` are TRUE. Returning shuffled result.')
      dr_df = dr_df[sample(1:nrow(dr_df)), ]
    }
    

    if(!is.null(facet_by) & show_background_cells){
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
      if(!is.null(foreground_subset_bool)[1]){
        dr_df = dr_df[dr_df$foreground_subset_bool, ]
      }
    }else if(show_background_cells){
      gg <- ggplot() + 
        geom_point(data = dr_df, aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha, color='grey80')
      if(!is.null(foreground_subset_bool)[1]){
        dr_df = dr_df[dr_df$foreground_subset_bool, ]
      }
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
        #if(any(pointdensity_subset %in% dr_df[, goi[i] ])){
        #  gg = gg + 
        #    ggpointdensity::geom_pointdensity(data = dr_df[dr_df[, goi[i]] %in% pointdensity_subset, ], mapping = aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha) +
        #    coord_fixed() +  
        #    scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
        #    scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        #}else{
          gg = gg + 
            ggpointdensity::geom_pointdensity(data = dr_df, mapping = aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        #}
      }else{
        #if(any(pointdensity_subset %in% dr_df[, goi[i] ])){
        #  gg = ggplot() + 
        #    ggpointdensity::geom_pointdensity(data =  dr_df[dr_df[, goi[i]] %in% pointdensity_subset, ], mapping = aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha) +
        #    coord_fixed() +  
        #    scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
        #    scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        #}else{
          gg = ggplot() + 
            ggpointdensity::geom_pointdensity(data = dr_df, mapping = aes(x=dim_1, y=dim_2), size=pt.size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        #}
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
    if(use_no_legend){
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
    
    if(!is.null(label_type) & cont_or_disc[i] == 'd'){
      #gg$data[, goi[i]] = as.factor(gg$data[, goi[i]])
      #gg = .LabelClusters(gg, goi[i], box = T, col_use = box_col)
      #centroid_df = get_centroids(dr_df, id = goi[i], id_subset = label_subset, facet_by = facet_by)
      if(!is.null(label_col)){
        custom_colors = label_col
      }
      if(is.null(custom_colors)){
        custom_colors = 'grey50'
      }
      gg = add_labels(gg = gg, df = dr_df, id = goi[i], id_subset = label_subset,
                      facet_by = facet_by, col_vec = custom_colors, label_type = label_type, alpha = 0.9)    
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


# #TODO replace repetitions of code with this helper function
# cget_scale_midpoint = function(feat){
#   (max(feat, na.rm = T) + min(feat, na.rm = T)) / 2 
# }

#' @export
jj_add_labels = function(gg, df, id, id_subset=NULL, facet_by = NULL, col_vec='grey50', label_type='geom_label_repel', alpha=0.9){
  #gg ggplot object
  #df data.frame with columns in the order: x-coordinate, y-coordinate, id-column, (facet-column)
  #id column in df for which labels should be created
  #id_subset Values in the id column, for which labels should be created. If not specified, plot labels for all values in the id column
  #facet_by if plot has facets, the facet column name should be specified here
  #col_vec either one colour for all text/boxes, or a named vector of colors for each single value in id
  #label type geom_text, geom_label, geom_text_repel or geom_label_repel
  #alpha alpha level for the boxes/text
  get_centroids = function(df, id, id_subset=NULL, facet_by = NULL){
    stopifnot(id %in% colnames(df))
    stopifnot(facet_by %in% colnames(df))
    get_medians = function(df){df %>% dplyr::group_by(group) %>% dplyr::summarise(dim1=median(dim1), dim2=median(dim2)) }
    #id_subset: user defined subset of id values to highlight
    if(!is.null(facet_by)){
      df = as.data.frame(df[,c(colnames(df)[1:2], id, facet_by)])
      colnames(df) = c('dim1','dim2','group', 'splitvar')
      facet_list = list()
      df_list = split(df, df$splitvar)
      for(i in seq_along(df_list)){
        df_list[[i]] = get_medians(df_list[[i]])
        df_list[[i]][, facet_by] = names(df_list)[i]
      }
      centroid_df = do.call(rbind, df_list)
    }else{
      df = as.data.frame(df[,c(colnames(df)[1:2], id)])
      colnames(df) = c('dim1','dim2','group')
      centroid_df = get_medians(df)
    }
    
    if(!is.null(id_subset)){
      stopifnot(all(id_subset %in% centroid_df$group))
      centroid_df = centroid_df[centroid_df$group %in% id_subset, ]
    }
    
    centroid_df = dplyr::select(centroid_df, dim1, dim2, group, everything()) %>% as.data.frame
    centroid_df
  }
  label_type = match.arg(label_type, choices = c('geom_text', 'geom_text_repel', 'geom_label', 'geom_label_repel'))
  label_df = get_centroids(df, id = id, id_subset = id_subset, facet_by = facet_by)
  library(ggrepel)
  if(label_type == 'geom_label_repel'){
    if(length(col_vec) == 1){
      gg = gg + geom_label_repel(data = label_df, mapping = aes(x=dim1, y = dim2, label=group), fill = col_vec, alpha=alpha, show.legend = F)
    }else{
      gg = gg + geom_label_repel(data = label_df, mapping = aes(x=dim1, y = dim2, label=group, fill=group), alpha=alpha, show.legend = F) + 
        scale_fill_manual(values = col_vec)
    }
  }else if(label_type == 'geom_text_repel'){
    if(length(col_vec) == 1){
      gg = gg + geom_text_repel(data = label_df, mapping = aes(x=dim1, y = dim2, label=group), colour = col_vec, alpha=alpha, show.legend = F)
    }else{
      gg = gg + geom_text_repel(data = label_df, mapping = aes(x=dim1, y = dim2, label=group, colour=group), alpha=alpha, show.legend = F) + 
        scale_fill_manual(values = col_vec)
    }
  }else if(label_type == 'geom_label'){
    if(length(col_vec) == 1){
      gg = gg + geom_label(data = label_df, mapping = aes(x=dim1, y = dim2, label=group), fill = col_vec, alpha=alpha, show.legend = F)
    }else{
      gg = gg + geom_label(data = label_df, mapping = aes(x=dim1, y = dim2, label=group, fill=group), alpha=alpha, show.legend = F) + 
        scale_fill_manual(values = col_vec)
    }
  }else if(label_type == 'geom_text'){
    if(length(col_vec) == 1){
      gg = gg + geom_text(data = label_df, mapping = aes(x=dim1, y = dim2, label=group), colour = col_vec, alpha=alpha, show.legend = F)
    }else{
      gg = gg + geom_text(data = label_df, mapping = aes(x=dim1, y = dim2, label=group, colour=group), alpha=alpha, show.legend = F) + 
        scale_fill_manual(values = col_vec)
    }
  }
  gg
}
