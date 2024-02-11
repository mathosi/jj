#' plot features
#'
#' Extension of the functionalities provided by Seurats FeaturePlot/DimPlot functions
#' Can plot scatterplots of multiple features provided under `features` (from Seurat assay) or `meta_features` (from Seurat meta.data or data.frame)
#' Embeddings from Seurat need to be specified under `reduction`. Alternatively, a data.frame with the embedding coordinates in the first two columns 
#' can be passed to `reduction` and seurat_obj can be set to NULL.
#' Loads ggplot2 and dplyr libraries.
#' 
#' jj_plot_scatter is a wrapper function that calls jj_plot_features and allows 2D scatterplots with an optional linear regression line.
#' jj_add_labels can be used to add custom labels to a ggplot
#' 
#' @param seurat_obj optional Seurat object
#' @param reduction either a data.frame with embeddings in column 1+2 and further meta data or a string specifying a dimensionality reduction from the Seurat object
#' @param features If Seurat object is provided, extract features from the specified assay and slot
#' @param meta_features if Seurat object is provided, extract features from the meta.data. Otherwise extract features from the data.frame provided in `reduction`
#' @param assay assay to use from the Seurat object, default: DefaultAssay(seurat_obj)
#' @param slot slot to use from the Seurat object, default: data
#' @param color_scale color scale to use for continuous features, one of 'wbr', 'bry', 'seurat', 'viridis'
#' @param cap_top upper value threshold passed to jj_cap_vals
#' @param cap_bottom lower value threshold passed to jj_cap_vals
#' @param custom_colors named vector of colours used to colour categorical features
#' @param custom_theme custom theme to be used, default: theme_minimal()
#' @param shape point shape, default: 16
#' @param alpha alpha value, default: 1
#' @param pt_size size of points, default: 'auto'
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
#' @param do_shuffle bool with default = TRUE, randomly shuffle the plotting order. Supersedes do_order.
#' @param label_type one of geom_text, geom_text_repel, geom_label, geom_label_repel
#' @param fill_colors colour to fill boxes or areas, if label_type or area_by is not NULL. If NULL, use colours from the custom_colors argument or grey if this is NULL as well
#' @param xlabel x axis label, default: UMAP 1
#' @param ylabel y axis label, default: UMAP 2
#' @keywords plot
#' @export
#' @examples
#' set.seed(1)
#' df = df = data.frame(x=c(rnorm(380, 5), rnorm(70,2), rnorm(50, 2.5, 2)), y=c(rnorm(380, 5), rnorm(70,-2), rnorm(50,2.5,2)), mock_gene = sample(0:10, size = 500, replace = T, prob = 11:1), clusterIdent = c(rep('gr1', 380), rep('gr2', 70), rep('gr3', 50)), batch = c(rep('A', 380), rep('B', 120)))
#' jj_plot_features(obj=df, features=c('mock_gene', 'clusterIdent'))
#' jj_plot_features(obj=df, facet_by = 'batch', use_pointdensity = T, pt_size=4, custom_theme = theme_bw())
#' jj_plot_features(obj=df, facet_by = 'batch', use_pointdensity = T, foreground_subset_bool = df$clusterIdent %in% c('gr1','gr2'), show_background_cells=T, pt_size=4, custom_theme = theme_bw())
#' jj_plot_features(obj=df, features = 'batch', facet_by = 'clusterIdent', show_background_cells=T, pt_size=4)
#' jj_plot_features(obj=df, features = 'batch', facet_by = 'clusterIdent', facet_subset = c('gr1', 'gr3'), pt_size=4)
#' jj_plot_features(obj = df, features = 'clusterIdent', facet_by = 'batch', label_type = 'geom_label_repel', label_subset = c('gr1','gr2'), custom_colors=structure(c('red','green', 'darkblue'), names = c('gr1','gr2','gr3')))
#' jj_plot_features(obj = df, features = 'clusterIdent', label_type = 'geom_label_repel', fill_colors = 'white', return_gg_object = T)[[1]] + labs(colour='new\ncolourtitle')
#' jj_plot_features(obj = df, features = 'clusterIdent', area_by = 'batch', custom_colors = c(gr1='red', gr2='green', gr3='darkblue'))
#' #only colour specific subsets
#' jj_plot_features(obj = df, features = 'clusterIdent', area_by = 'batch', custom_colors = c(gr1='red', gr2='green'), fill_colors = c(A = 'yellow'))
#' #if both area_by and label_type are specified, fill_colors will overwrite the label colours and set them to grey
#' jj_plot_features(obj = df, features = 'clusterIdent', area_by = 'batch', label_type = 'geom_label', custom_colors = c(gr3='purple'), fill_colors = c(gr1 = 'yellow', A='orange'))
#' 
#' jj_plot_features(obj= pbmc_small, reduction ='tsne', features=c('MS4A1', 'CD79A'), cap_top='q95', color_scale='wbr', xlabel='tsne_1', ylabel='tsne_2',my_title=c('tnse: MS4A1','tsne: CD79A'))
#' jj_arrange_ggplots(jj_plot_features(obj= pbmc_small, reduction ='tsne', features=c('MS4A1', 'CD79A', 'CD8A'), cap_top='auto', return_gg_object = TRUE))
#' jj_plot_features(pbmc_small, reduction ='tsne', features= 'nCount_RNA', cap_top='q95', foreground_subset_bool = pbmc_small$groups == 'g1', show_background_cells = T)
#' jj_plot_features(pbmc_small, features= 'nCount_RNA', cap_top='q95', area_by = 'letter.idents', area_thres = 0.9)

jj_plot_features <- function(obj,
                             features = NULL,
                             reduction = NULL, 
                             assay_order = c('meta.data|cellColData','DefaultAssay', 'RNA', 'GeneScoreMatrix', 'PeakMatrix'), 
                             slot = 'data', 
                             color_scale = c('viridis', 'wbr', 'gbr', 'bry', 'seurat', 'spectral')[1],
                             facet_by = NULL, 
                             area_by = NULL,
                             area_thres = 0.98,
                             n_facet_rows = NULL, 
                             facet_subset = NULL, 
                             cap_top = NULL,  
                             cap_bottom = NULL, 
                             custom_colors = NULL, 
                             custom_theme = theme_minimal(), 
                             pt_size = 'a', 
                             return_gg_object = FALSE, 
                             use_pointdensity = FALSE,
                             show_background_cells = FALSE,
                             foreground_subset_bool = NULL,
                             cont_or_disc = 'a', 
                             order_type = c('shuffle', 'order', 'none')[1],
                             convert_factors=FALSE,
                             my_title = NULL,
                             label_type=NULL, 
                             label_subset = NULL,
                             fill_colors=NULL, 
                             use_no_legend = FALSE,
                             xlabel = NULL,
                             ylabel = NULL,
                             shape = 16, 
                             alpha = 1){
  library(ggplot2)
  library(dplyr)
  #@param pointdensity_subset Only used if use_pointdensity=T. If NULL, use all cells. If set to groups within meta_features, only calculate density for these subgroups
  # if(is.null(reduction)){
  #   stop('reduction must be either string specifying the reduction to use from seurat object or a dr_df data.frame')
  # }
  message('Validating object')
  obj_type = class(obj)
  if(any(c('matrix','tbl','data.frame') %in% obj_type)){
    dr_df = as.data.frame(obj)
  }else if('Seurat' %in% obj_type){
    dr_df <- jj_get_reduction_coords(obj, reduction, exact_match = TRUE)
  }else if('ArchRProject' %in% obj_type){
    dr_df <- jj_get_reduction_coords(obj, reduction, exact_match = TRUE)
  }else{
    stop('obj must be data.frame, Seurat object or ArchRProject')
  }
 
  
  color_scale <- match.arg(color_scale, choices = c('viridis', 'wbr', 'gbr', 'bry', 'seurat', 'spectral'))
  order_type = match.arg(order_type, choices = c('shuffle', 'order', 'none'))
  if(!all(apply(dr_df[, c(1,2)], 2, is.numeric))){
    stop('First two columns of the data.frame should contain numeric values used as x- and y-coordinates')
  } 
  
  if(is.null(xlabel)) xlabel = colnames(dr_df)[1]
  if(is.null(ylabel)) ylabel = colnames(dr_df)[2]
  
  colnames(dr_df)[1:2] <- c('dim_1', 'dim_2')
  na_coords = !complete.cases(dr_df[, 1:2])
  if(any(na_coords)){
    if(sum(na_coords) == nrow(dr_df)){
      stop('All observations have NA in their coordinates.')
    }
    warning(sprintf('Removing %i missing values from coordinates.', sum(na_coords)))
    dr_df = dr_df[!na_coords, ]
  }
  
  if(is.null(features)){
    dr_df$data = ''
    features = 'data'
  }
  
  if(!is.null(foreground_subset_bool)[1]){
    if(!show_background_cells) warning('`foreground_subset_bool` will have no effect if `show_background_cells` = FALSE')
    dr_df$foreground_subset_bool = foreground_subset_bool
  }
  
  cont_or_disc <- unlist(strsplit(cont_or_disc, split = ''))
  if(length(cont_or_disc) == 1){
    cont_or_disc <- rep(cont_or_disc, length(features))
  }else if(length(cont_or_disc) != length(features)){
    stop('length of cont_or_disc must be either 1 or length(features)')
  }else if(any(!cont_or_disc %in% c('c', 'd', 'a'))){
    stop('cont_or_disc can only be c, d, a: continuous, discrete, automatic')
  }
  
  if(pt_size == 'a'){
    pt_size = ifelse(nrow(dr_df) > 20000, 0.5, ifelse(nrow(dr_df) > 10000, 1, ifelse(nrow(dr_df) > 1000, 1.5, 2)))
  }
  
  # if(assay_order[1] == 'meta.data|cellColData' | any(c('matrix','tbl','data.frame') %in% obj_type)){
  #   message('Getting features from metadata')
  #   if(!all(features %in% colnames(dr_df))){
  #     warning(sprintf('%s not found in the assay meta.data.',
  #                     paste(features[!features %in% colnames(dr_df)], collapse = ', ')))
  #     #features = features[features %in% colnames(dr_df)]
  #   }
  # }
    
  if('Seurat' %in% obj_type){
    ass_avail = Seurat::Assays(obj)
    assay_order = c(assay_order, ass_avail)
    for(i in seq_along(assay_order)){
      if(assay_order[i] == 'meta.data|cellColData'){
        if(any(features %in% colnames(dr_df))){
          message('Getting features from meta.data')
          break
        }else{
          next
        }
      }else if(assay_order[i] == 'DefaultAssay'){
        assay_use = DefaultAssay(obj)
      }else if(assay_order[i] %in% ass_avail){
        assay_use = assay_order[i]
      }else{
        next
      }
      ass_m = GetAssayData(obj, assay = assay_use, slot = slot)
      if(any(features %in% rownames(ass_m))){
        message(sprintf('Getting features from %s slot of the %s assay in the Seurat object', slot, assay_use))
        dr_df <- jj_bind_features_with_dr_df(ass_m, features=features[features %in% rownames(ass_m)], dr_df=dr_df[, 1:2], log10Transform=FALSE)
        break
      }
    }
  }else if('ArchRProject' %in% obj_type){
    ass_avail = getAvailableMatrices(obj)
    assay_order = c(assay_order, ass_avail)
    for(i in seq_along(assay_order)){
      if(assay_order[i] == 'meta.data|cellColData'){
        if(any(features %in% colnames(dr_df))){
          message('Getting features from meta.data')
          break
        }else{
          next
        }
      }else if(assay_order[i] %in% ass_avail){
        assay_use = assay_order[i]
      }else{
        next
      }
      if(assay_use == 'PeakMatrix'){
        favail = getPeakSet(obj)
        favail = try(paste(seqnames(favail), ranges(favail), sep = "-"))
      }else{
        favail = getFeatures(obj, useMatrix = assay_use)
      }
     
      if(any(features %in% favail)){
        message(sprintf('Getting features from %s assay in the ArchR object', assay_use))
        ass_m = get_archr_mat(obj, mat = assay_use, assay_use = 1)
        dr_df <- jj_bind_features_with_dr_df(ass_m, features=features[features %in% rownames(ass_m)], dr_df=dr_df[, 1:2], log10Transform=FALSE)
        break
      }
    }
  }else{
    assay_order = 'meta.data|cellColData'
    i = 1
  }
  
  #replace characters that can lead to problems in ggplot2
  features <- gsub('-', '_', features)
  colnames(dr_df) <- gsub('-', '_', colnames(dr_df))
  features <-  gsub('::', '___', features)
  colnames(dr_df) <-  gsub('::', '___', colnames(dr_df))
  features <-  gsub(':', '__', features)
  colnames(dr_df) <-  gsub(':', '__', colnames(dr_df))
  features <-  gsub('?', '', features, fixed = T)
  colnames(dr_df) <-  gsub('?', '', colnames(dr_df), fixed = T)
  
  if(!any(features %in% colnames(dr_df))){
    message('None of the features was found in the provided assay_order')
    return(NULL)
  }else if(!all(features %in% colnames(dr_df))){
    message(sprintf('Features missing in the assay %s:\n%s', assay_order[i], paste(features[!features %in% colnames(dr_df)], collapse = ', ')))
    features = features[features %in% colnames(dr_df)]
  }
  
  if(convert_factors & assay_order[i] == 'meta.data|cellColData'){
    #convert_factors to numeric if possible and otherwise to character
    for(me in features){
      if(is.factor(dr_df[, me])){
        dr_df[, me] = as.character(dr_df[, me])
        if(Hmisc::all.is.numeric(dr_df[, me])){
          dr_df[, me] <- as.numeric(dr_df[, me])
        }
      }
    }
  }
  
  #try to find continous/discrete status for each variable automatically
  if(any(c('a','auto') %in% cont_or_disc) & assay_order[i] == 'meta.data|cellColData'){
    cont_or_disc = vector()
    #cont_or_disc = paste(ifelse(sapply(dr_df[, colnames(dr_df) %in% features], class) == 'character', 'd', 'c'), collapse = '')
    for(i in features){
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
    message('Continous or discrete feature order: ', cont_or_disc)
  }
  
  if(!is.null(cap_top) | !is.null(cap_bottom)){
    message('Capping meta feature values')
    for(i in features){
      if(is.numeric(dr_df[, i])){
        dr_df[, i] <- jj_cap_vals(dr_df[, i], cap_top=cap_top,  cap_bottom=cap_bottom)
      }
    }
  }


  range_x <- range(dr_df$dim_1)
  range_y <- range(dr_df$dim_2)
  
  gg_list <- list()
  message('plotting features')
  for(i in seq_along(features)){
    message(sprintf('%i/%i: %s', i, length(features), features[i]))
    
    # if(is.factor(dr_df[, features[i]])){
    #   dr_df[, features[i]] <- as.character(dr_df[, features[i]])
    #   if(Hmisc::all.is.numeric(dr_df[, features[i]]) & n_distinct(dr_df[, features[i]]) >=30){
    #     message("Converting factor to numeric")
    #     dr_df[, features[i]] <- as.numeric(dr_df[, features[i]])
    #   }else{
    #     message("Converting factor to character.")
    #   }
    # }
    
    if(cont_or_disc[i] == 'd'){
      #dr_df[, features[i] ] = as.character(dr_df[, features[i] ])
      if(!is.null(custom_colors)){
        if(is.null(names(custom_colors))){
          names(custom_colors) = levels(as.factor(dr_df[, features[i] ]))
        }
        breaks_use = names(custom_colors)
        dr_df[, features[i] ] <- factor(dr_df[, features[i] ], levels= breaks_use)
      #}else if(suppressWarnings(sum(is.na(as.numeric(names(jj_get_jj_colours(levels(as.factor(dr_df[, features[i] ]))))))) == 1 )){
      #  breaks_use <- suppressWarnings(as.character(sort(as.numeric(names(jj_get_jj_colours(dr_df[, features[i]]))), decreasing = F)))
      #  dr_df[, features[i] ] <- factor(dr_df[, features[i] ], levels= breaks_use)
      }else{
        breaks_use <- names(jj_get_jj_colours(levels(as.factor(dr_df[, features[i]]))))
      }
    }
    
    if(order_type == 'order'){
      dr_df = dr_df[order(dr_df[, features[i]]), ]
    }else if(order_type == 'shuffle'){
      #if(do_order) warning('Both `do_order` and `do_shuffle` are TRUE. Returning shuffled result.')
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
        geom_point(data = dr_df_background, aes(x=dim_1, y=dim_2), size=pt_size, alpha=alpha, color='grey80')
      if(!is.null(foreground_subset_bool)[1]){
        dr_df = dr_df[dr_df$foreground_subset_bool, ]
      }
    }else if(show_background_cells){
      gg <- ggplot() + 
        geom_point(data = dr_df, aes(x=dim_1, y=dim_2), size=pt_size, alpha=alpha, color='grey80')
      if(!is.null(foreground_subset_bool)[1]){
        dr_df = dr_df[dr_df$foreground_subset_bool, ]
      }
    }else{
      gg = NULL
    }
    
    if(!is.null(area_by)){
      stopifnot(area_by %in% colnames(dr_df))
      cont_df = jj_get_contour_lines(dr_df, area_by, area_thres)
      if(!is.null(gg)){
        gg = gg + geom_polygon(data = cont_df, aes_string(x='x', y='y', fill = area_by, 
                                                          group = 'cont_group'), alpha = 0.5, size=1)
        #gg = gg + geom_path(data = cont_df, aes_string(x='x', y='y', colour = area_by, group = 'cont_group'), size=1)
      }else{
        gg = ggplot() + geom_polygon(data = cont_df, aes_string(x='x', y='y', fill = area_by, 
                                                                group = 'cont_group'), alpha = 0.5, size=1)
        #gg = ggplot() + geom_path(data = cont_df, aes_string(x='x', y='y', colour = area_by, group = 'cont_group'), size=1)
      }
    }
    
    if(!is.null(facet_by) & !is.null(facet_subset[1])){
      dr_df = dr_df[dr_df[, facet_by] %in% facet_subset, ]
      if(!nrow(dr_df) > 0){
        stop('Groups in `facet_subset` are not available in `facet_by`')
      }
    }
    
    if(use_pointdensity){
      if(!is.null(gg)){
        #if(any(pointdensity_subset %in% dr_df[, features[i] ])){
        #  gg = gg + 
        #    ggpointdensity::geom_pointdensity(data = dr_df[dr_df[, features[i]] %in% pointdensity_subset, ], mapping = aes(x=dim_1, y=dim_2), size=pt_size, alpha=alpha) +
        #    coord_fixed() +  
        #    scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
        #    scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        #}else{
          gg = gg + 
            ggpointdensity::geom_pointdensity(data = dr_df, mapping = aes(x=dim_1, y=dim_2), size=pt_size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        #}
      }else{
        #if(any(pointdensity_subset %in% dr_df[, features[i] ])){
        #  gg = ggplot() + 
        #    ggpointdensity::geom_pointdensity(data =  dr_df[dr_df[, features[i]] %in% pointdensity_subset, ], mapping = aes(x=dim_1, y=dim_2), size=pt_size, alpha=alpha) +
        #    coord_fixed() +  
        #    scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
        #    scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        #}else{
          gg = ggplot() + 
            ggpointdensity::geom_pointdensity(data = dr_df, mapping = aes(x=dim_1, y=dim_2), size=pt_size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        #}
      }
    }else{
      if(is.character(shape)){
        if(!is.null(gg)){
          gg <- gg +  
            geom_point(data=dr_df, aes_string(x='dim_1', y='dim_2', colour=features[i], shape=shape), size=pt_size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }else{
          gg = ggplot(dr_df, aes(x=dim_1, y=dim_2)) + 
            geom_point(aes_string(colour=features[i], shape=shape), size=pt_size, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }
      }else{
        if(!is.null(gg)){
          gg <- gg + 
            geom_point(data=dr_df, aes_string(x='dim_1', y='dim_2', colour=features[i]), size=pt_size, shape=shape, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }else{
          gg <- ggplot(dr_df, aes(x=dim_1, y=dim_2)) + 
            geom_point(aes_string(colour=features[i]), size=pt_size, shape=shape, alpha=alpha) +
            coord_fixed() +  
            scale_x_continuous(breaks=seq(floor(range_x[1]),ceiling(range_x[2]),2)) + 
            scale_y_continuous(breaks=seq(floor(range_y[1]),ceiling(range_y[2]),2))
        }
      }
    }
    
    if(use_pointdensity){
      gg = gg + viridis::scale_color_viridis()
    }else if(!all(is.null(custom_colors))){ #& !n_distinct(dr_df[, features[i]]) < 30){
      if(length(custom_colors)==3 & cont_or_disc[i] == 'c'){
        #use gradient if 3 colours are supplied
        mean_acc <- (max(dr_df[, features[i]], na.rm = T) + min(dr_df[, features[i]], na.rm = T)) / 2 #mean(dr_df[, features[i]])
        gg <- gg + scale_color_gradient2(low = custom_colors[1], mid = custom_colors[2], high = custom_colors[3], midpoint = mean_acc)
      }else if(!is.null(names(custom_colors))){
        #use discrete color mapping if named vector is supplied
        #custom_colors <- custom_colors[names(custom_colors) %in% dr_df[, features[i]]]
        message('Using colors provided unter custom_colors.')
        gg <- gg +
          scale_colour_manual(values = custom_colors[names(custom_colors) %in% dr_df[, features[i]]]) +#, breaks = names(custom_colors)) + 
          guides(colour = guide_legend(override.aes = list(size=5)))
      }
    }else if(cont_or_disc[i] == 'd'){
      #TODO: can error if discrete number of values is > =30
      #plot with discrete colors if less than 30 distinct values in variable (or if d is set (check for <60 distinct groups as sanity check))
      gg <- gg +
        scale_colour_manual(values = jj_get_jj_colours(breaks_use)) + 
        guides(colour = guide_legend(override.aes = list(size=5)))
    }else if(color_scale=='viridis'){
      gg <- gg + viridis::scale_color_viridis()
    }else if(color_scale=='seurat'){
      gg <- gg + scale_color_gradient(low='#C3C3C3', high = '#1C0DFD')
    }else if(color_scale=='bry'){
      mean_acc <- (max(dr_df[, features[i]],na.rm = T)+ min(dr_df[, features[i]],na.rm = T)) / 2 #mean(dr_df[, features[i]])
      gg <- gg +  scale_color_gradient2(low = "darkblue", mid = "red", high = "yellow", midpoint = mean_acc)
    }else if(color_scale=='wbr'){
      mean_acc <- (max(dr_df[, features[i]], na.rm = T) + min(dr_df[, features[i]], na.rm = T)) / 2 #mean(dr_df[, features[i]])
      #print(mean_acc)
      gg <- gg + scale_color_gradient2(low = "#ffffd9", mid = "blue", high = "red", midpoint = mean_acc)
    }else if(color_scale=='gbr'){
      mean_acc <- (max(dr_df[, features[i]], na.rm = T) + min(dr_df[, features[i]], na.rm = T)) / 2 #mean(dr_df[, features[i]])
      #print(mean_acc)
      gg <- gg + scale_color_gradient2(low = "grey80", mid = "blue", high = "red", midpoint = mean_acc)
    }else if(color_scale=='spectral'){
      gg <- gg + scale_color_gradientn(colours = c('#5E4FA2', '#3288BD', '#66C2A5', '#ABDDA4', '#E6F598', '#FFFFBF', 
                                                   '#FEE08B', '#FDAE61', '#F46D43', '#D53E4F', '#9E0142'))
    }
    
    gg = gg + xlab(xlabel) + ylab(ylabel)
    
    if(is.theme(custom_theme)){
      gg <- gg + custom_theme
    }
    if(use_no_legend){
      gg = gg + theme(legend.position='none')
    }
    if(!is.null(my_title)){
      if(length(my_title) == length(features)){
        gg <- gg + ggtitle(my_title[i]) +  theme(plot.title = element_text(hjust = 0.5))
      }
      else if(length(my_title) == 1){
        gg <- gg + ggtitle(my_title) +  theme(plot.title = element_text(hjust = 0.5))
      }
    }
    
    if(!is.null(label_type) & cont_or_disc[i] == 'd'){
      #gg$data[, features[i]] = as.factor(gg$data[, features[i]])
      #gg = .LabelClusters(gg, features[i], box = T, col_use = box_col)
      #centroid_df = get_centroids(dr_df, id = features[i], id_subset = label_subset, facet_by = facet_by)
      if(!is.null(fill_colors)){
        custom_colors = fill_colors
      }
      if(is.null(custom_colors)){
        custom_colors = 'grey50'
      }
      gg = jj_add_labels(gg = gg, df = dr_df, id = features[i], id_subset = label_subset,
                      facet_by = facet_by, col_vec = custom_colors, label_type = label_type, alpha = 0.9)    
    }
    if(!is.null(area_by) & !is.null(fill_colors)){
      gg = gg + scale_fill_manual(values = fill_colors[names(fill_colors) %in% dr_df[, area_by]]) + 
        guides(fill = guide_legend(override.aes = list(size=5)))
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
    #   if(facet_by & n_distinct(dr_df[, features[i]]) < 100){
    #     message(sprintf('Facetting by %s.', features[i]))
    #     if(is.null(n_facet_rows)){
    #       n_facet_rows <- round(sqrt(n_distinct(dr_df[, features[i]])))
    #     }
    #     if(!is.logical(facet_by)) print(gg) #print unfacetted + facetted if facet_by <- 1 is used
    #     gg <- gg + facet_wrap(as.formula(sprintf('.~%s', features[i])), nrow=n_facet_rows)
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


#' @export
jj_plot_scatter = function(df, x, y, col=NULL, add_regression_line=FALSE, ...){
  stopifnot(all(c(x,y) %in% colnames(df)))
  if(!is.null(col)){
    stopifnot(col %in% colnames(df))
  }else{
    df$group = " "
    col = 'group'
  }
  df = df[, c(x, y, col)]
  colnames(df)[1:2] = c('x', 'y')
  gg = jj_plot_features(obj=df, 
                        features = col,
                        xlabel = x, 
                        ylabel = y,
                        return_gg_object = T,
                        ...)[[1]]
  
  #overwrite axes system and labeling
  gg = suppressMessages(gg + coord_cartesian() + 
    scale_x_continuous(breaks = waiver()) + 
    scale_y_continuous(breaks = waiver()))
  
  if(add_regression_line){
    # GET EQUATION AND R-SQUARED AS STRING
    # SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
    lm_eqn <- function(df){
      m <- lm(y ~ x, df);
      eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                       list(a = format(unname(coef(m)[1]), digits = 2),
                            b = format(unname(coef(m)[2]), digits = 2),
                            r2 = format(summary(m)$r.squared, digits = 3)))
      as.character(as.expression(eq));
    }
    gg <- gg + 
      annotate('text',  x = -Inf, y = Inf, hjust = -0.5, vjust = 1, label = lm_eqn(df), parse = TRUE) + 
      geom_smooth(method = "lm", se=T, color="black", formula = y ~ x)
  }
  gg
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