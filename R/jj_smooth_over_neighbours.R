#' smooth values using nearest neighbours
#'
#' smooth values in a data.frame by taking the mean (numeric columns) or mode (character columns) among the k nearest neighbours for each cell (row)
#'
#' @param dr_df data.frame to smooth values over 
#' @param reduction_df data.frame with embeddings from dimensionality reduction (same nrow as dr_df)
#' @param k number of nearest neighbours to smooth over for each cell
#' @export
#' @examples
#' emb_df = Embeddings(pbmc_small, 'tsne')
#' meta_df = pbmc_small[[]]
#' res = jj_summarize_nearest_neighbour_data(dr_df = meta_df, reduction_df = emb_df, k = 5)
#' pbmc_small$nCount_RNA_smoothed = res$nCount_RNA
#' pbmc_small$groups_smoothed = res$groups
#' jj_arrange_ggplots(jj_plot_features(pbmc_small, reduction = 'tsne', 
#'                                    meta_features = c('nCount_RNA','nCount_RNA_smoothed', 'groups', 'groups_smoothed'), 
#'                                    pt.size = 2, return_gg_object = T))

jj_summarize_nearest_neighbour_data <- function(dr_df, reduction_df, k=30){
  #dr_df: metadata to summarize: indices in nn_mat must match (dr_df and nn_mat generated from same seurat object!)
  #nn_mat: nearest neighbour matrix with selected cells as rownames and indices of nearest neighbours in the row values
  nn_mat = get_nn_mat(reduction_df, k = k)
  
  stopifnot(!anyNA(dr_df))
  #stopifnot('cells' %in% colnames(dr_df))
  #add indices of selected cells to the neighbour matrix in the first row
  neighbour_mat <- t(cbind(match(rownames(nn_mat), rownames(dr_df)), nn_mat))
  stopifnot(!anyNA(neighbour_mat))
  #summarize numeric and non-numeric values separately
  numeric_or_not <- sapply(dr_df,  is.numeric)
  new_cols <- c(colnames(dr_df)[!numeric_or_not], colnames(dr_df)[numeric_or_not]) 
  aggr_df <- jj_initialize_df(ncol = ncol(dr_df), nrow = nrow(nn_mat),
                            init = NA, row.names = rownames(nn_mat),
                            col.names = new_cols)
  
  for(i in 1:ncol(neighbour_mat)){
    if(i %% 500 == 0) message(sprintf('%.2f%%', i/ncol(neighbour_mat)*100))
    #get subset with selected neighbours for a given cell
    dr_df_temp <- dr_df[neighbour_mat[, i],  ]
    #calculate means for all numeric variables
    numeric_aggr <- apply(dr_df_temp[, numeric_or_not], 2, mean)
    #select most abundant value for non-numeric variables
    not_numeric_aggr <- apply(dr_df_temp[, !numeric_or_not], 2, function(x) names(sort(table(as.character(x)),decreasing=TRUE)[1]))
    #add to the aggregation matrix
    aggr_df[i, ] <- c(not_numeric_aggr, numeric_aggr)
  }
  if(length(numeric_aggr)>0){
    aggr_df[, (length(not_numeric_aggr)+1):ncol(aggr_df)] <- apply(aggr_df[, (length(not_numeric_aggr)+1):ncol(aggr_df)], 2, as.numeric) 
  }
  stopifnot(!anyNA(aggr_df))
  return(aggr_df)
}

get_nn_mat = function(reduction_df, k){
  nn_map <- FNN::knn.index(reduction_df, k = k) #data.frame of cell embeddings
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  row.names(nn_map) <- row.names(reduction_df)
  return(nn_map)
}


get_mean_expression = function(seurat_obj, gene, red_use='umap', k=5){
  #get mean expression values from `k` neighbours of each cell based on `red_use` reduction
  DefaultAssay(seurat_obj) = 'RNA'
  stopifnot(gene %in% rownames(seurat_obj))
  reduced_coordinates = Embeddings(seurat_obj, red_use)
  nn_map <- FNN::knn.index(reduced_coordinates, k = k)
  nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
  row.names(nn_map) <- row.names(reduced_coordinates)
  goi = GetAssayData(seurat_obj, assay = 'RNA', slot = 'data')[gene, ]
  mean_goi = sapply(1:nrow(nn_map), function(x) mean(goi[nn_map[x, ]]))
  summary(mean_goi)
  return(mean_goi)
}

# annotate_from_nearest_neighbours = function(dr_df, annot_col, disc_or_cont = 'd', only_na = TRUE, n_neighbours=30){
#   #fill in annotation in `annot_col` in dr_df based on the annotations of nearest neighbours
#   #if `only_na` = TRUE, only fill in the cells that have a NA value (otherwise all annotations are smoothed)
#   nn_map <- FNN::knn.index(dr_df[, 1:2], k = n_neighbours) #columns 1 and 2 must contain cell embeddings
#   nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
#   row.names(nn_map) <- row.names(dr_df)
#   annot_col_vals = dr_df[, annot_col]
#   if(disc_or_cont=='d'){
#     annot_val = lapply(1:nrow(nn_map), function(x) names(sort(table(annot_col_vals[nn_map[x, ]]), decreasing = T))[1])
#     #if all neighbours are NA, then Null is returned and must be replaced by NA before using unlist to keep vector length
#     annot_val = unlist(lapply(annot_val, function(x) ifelse(is.null(x), NA, x))  )
#   }else if(disc_or_cont=='c'){
#     annot_val = sapply(1:nrow(nn_map), function(x) mean(annot_col_vals[nn_map[x, ]]))
#   }else{
#     stop("disc_or_cont must be either d or c")
#   }
#   if(only_na){
#     annot_col_vals[is.na(annot_col_vals)] = annot_val[is.na(annot_col_vals)]
#     annot_val = annot_col_vals
#   }
#   #if there are still na values, rerun a second round only getting the indices of cells with known annotations
#   all_na_neighbours = is.na(annot_val) #which(sapply(1:nrow(nn_map), function(x) all(is.na(annot_col_vals[nn_map[x, ]]))))
#   if(sum(all_na_neighbours) > 0){
#     #na_positions = dr_df[all_na_neighbours, 1:2]
#     nn_map2 <- FNN::knnx.index(data=dr_df[!all_na_neighbours, 1:2],
#                                query=dr_df[all_na_neighbours, 1:2], k = n_neighbours)
#     nn_map2 <- cbind(nn_map2, seq_len(nrow(nn_map2)))
#     row.names(nn_map2) <- row.names(dr_df[all_na_neighbours, ])
#     annot_col_vals = annot_val[!all_na_neighbours]
#     if(disc_or_cont=='d'){
#       annot_val2 = sapply(1:nrow(nn_map2), function(x) names(sort(table(annot_col_vals[nn_map2[x, ]]), decreasing = T))[1])
#     }else if(disc_or_cont=='c'){
#       annot_val2 = sapply(1:nrow(nn_map2), function(x) mean(annot_col_vals[nn_map2[x, ]]))
#     }
#     annot_val[is.na(annot_val)] = annot_val2
#   }   
#   return(annot_val)
# }



# make_cicero_cds_return_nn <- function (cds, reduced_coordinates, k = 50, summary_stats = NULL, 
#                                        size_factor_normalize = TRUE, silent = FALSE, iterations = 5000) 
# {
#   #function from cicero, that was modified to also return the k-nearest neighbour matrix and be more verbose
#   assertthat::assert_that(is(cds, "CellDataSet"))
#   assertthat::assert_that(is.data.frame(reduced_coordinates) | 
#                             is.matrix(reduced_coordinates))
#   assertthat::assert_that(assertthat::are_equal(nrow(reduced_coordinates), 
#                                                 nrow(pData(cds))))
#   assertthat::assert_that(setequal(row.names(reduced_coordinates), 
#                                    colnames(cds)))
#   assertthat::assert_that(assertthat::is.count(k) & k > 1)
#   assertthat::assert_that(is.character(summary_stats) | is.null(summary_stats))
#   if (!is.null(summary_stats)) {
#     assertthat::assert_that(all(summary_stats %in% names(pData(cds))), 
#                             msg = paste("One of your summary_stats is missing", 
#                                         "from your pData table. Either add a", "column with the name in", 
#                                         "summary_stats, or remove the name", "from the summary_stats parameter.", 
#                                         collapse = " "))
#     assertthat::assert_that(sum(vapply(summary_stats, function(x) {
#       !(is(pData(cds)[, x], "numeric") | is(pData(cds)[, 
#                                                        x], "integer"))
#     }, 1)) == 0, msg = paste("All columns in summary_stats must be", 
#                              "of class numeric or integer.", collapse = " "))
#   }
#   assertthat::assert_that(is.logical(size_factor_normalize))
#   assertthat::assert_that(is.logical(silent))
#   reduced_coordinates <- as.data.frame(reduced_coordinates)
#   reduced_coordinates <- reduced_coordinates[colnames(cds), 
#   ]
#   #return the indices of the k nearest neigbours (each row is one barcode, each column the index of the neigherst neighbour)
#   #with k=50, matrix of n cells * k neaerest neighbours
#   nn_map <- FNN::knn.index(reduced_coordinates, k = (k - 1))
#   row.names(nn_map) <- row.names(reduced_coordinates)
#   #add 1:nrow(cds) to matrix to include the index of the original cell
#   nn_map <- cbind(nn_map, seq_len(nrow(nn_map)))
#   #initialize chosen cells
#   good_choices <- seq_len(nrow(nn_map))
#   #draw index of 1 cell randomly
#   choice <- sample(seq_len(length(good_choices)), size = 1, 
#                    replace = FALSE)
#   #add it to the vector of chosen cells
#   chosen <- good_choices[choice]
#   #remove the cell from the remaining cells to choose from subsequently
#   good_choices <- good_choices[good_choices != good_choices[choice]]
#   it <- 0
#   k2 <- k * 2
#   #get the number of common nearest neigbours for two cells
#   get_shared <- function(other, this_choice) {
#     k2 - length(union(cell_sample[other, ], this_choice))
#   }
#   
#   #in 5000 iterations (and as long as there are cells to choose from randomly), 
#   #add a randomly drawn cell to the vector of chosen cells if the number of common neighbours
#   #is < 0.9 * k in a comparison with any of the chosen cells so far
#   message('sampling cells that are not too close to each other')
#   while (length(good_choices) > 0 & it < iterations) {
#     if(it %% 500 == 0) message(sprintf('%i/%i iterations', it, iterations))
#     it <- it + 1
#     choice <- sample(seq_len(length(good_choices)), size = 1, 
#                      replace = FALSE)
#     new_chosen <- c(chosen, good_choices[choice])
#     good_choices <- good_choices[good_choices != good_choices[choice]]
#     cell_sample <- nn_map[new_chosen, ]
#     others <- seq_len(nrow(cell_sample) - 1)
#     this_choice <- cell_sample[nrow(cell_sample), ]
#     shared <- sapply(others, get_shared, this_choice = this_choice)
#     if (max(shared) < 0.9 * k) {
#       chosen <- new_chosen
#     }
#   }
#   
#   #get nearest neighbour matrix for all chosen cells
#   cell_sample <- nn_map[chosen, ]
#   if (!silent) {
#     message('calculating statistics')
#     #Generate all combinations of the elements of x taken m at a time.
#     #eg combn(1:3, 2) -> 
#     #       [,1] [,2] [,3]
#     # [1,]    1    1    2
#     # [2,]    2    3    3
#     combs <- combn(nrow(cell_sample), 2)
#     #for every pair of cells, get the number of common nearest neighbours 
#     shared <- apply(combs, 2, function(x) {
#       k2 - length(unique(as.vector(cell_sample[x, ])))
#     })
#     message(paste0("Overlap QC metrics:\nCells per bin: ", 
#                    k, "\nMaximum shared cells bin-bin: ", max(shared), 
#                    "\nMean shared cells bin-bin: ", mean(shared), "\nMedian shared cells bin-bin: ", 
#                    median(shared)))
#     if (mean(shared)/k > 0.1) 
#       warning("On average, more than 10% of cells are shared between paired bins.")
#   }
#   message('calculate new expression matrix')
#   exprs_old <- exprs(cds)
#   mask <- sapply(seq_len(nrow(cell_sample)), function(x) seq_len(ncol(exprs_old)) %in% 
#                    cell_sample[x, , drop = FALSE])
#   mask <- Matrix::Matrix(mask)
#   new_exprs <- exprs_old %*% mask
#   new_exprs <- Matrix::t(new_exprs)
#   new_exprs <- as.matrix(new_exprs)
#   pdata <- pData(cds)
#   new_pcols <- "agg_cell"
#   if (!is.null(summary_stats)) {
#     new_pcols <- c(new_pcols, paste0("mean_", summary_stats))
#   }
#   message('Calculate means for metadata')
#   new_pdata <- plyr::adply(cell_sample, 1, function(x) {
#     sub <- pdata[x, ]
#     df_l <- list()
#     df_l["temp"] <- 1
#     for (att in summary_stats) {
#       df_l[paste0("mean_", att)] <- mean(sub[, att])
#     }
#     data.frame(df_l)
#   })
#   message('Construct CellDataSet')
#   new_pdata$agg_cell <- paste("agg", chosen, sep = "")
#   new_pdata <- new_pdata[, new_pcols, drop = FALSE]
#   row.names(new_pdata) <- new_pdata$agg_cell
#   row.names(new_exprs) <- new_pdata$agg_cell
#   new_exprs <- as.matrix(t(new_exprs))
#   fdf <- fData(cds)
#   new_pdata$temp <- NULL
#   fd <- new("AnnotatedDataFrame", data = fdf)
#   pd <- new("AnnotatedDataFrame", data = new_pdata)
#   cicero_cds <- suppressWarnings(newCellDataSet(new_exprs, 
#                                                 phenoData = pd, featureData = fd, expressionFamily = negbinomial.size(), 
#                                                 lowerDetectionLimit = 0))
#   cicero_cds <- monocle::detectGenes(cicero_cds, min_expr = 0.1)
#   cicero_cds <- BiocGenerics::estimateSizeFactors(cicero_cds)
#   if (any(!c("chr", "bp1", "bp2") %in% names(fData(cicero_cds)))) {
#     fData(cicero_cds)$chr <- NULL
#     fData(cicero_cds)$bp1 <- NULL
#     fData(cicero_cds)$bp2 <- NULL
#     fData(cicero_cds) <- cbind(fData(cicero_cds), row.names(fData(cicero_cds)))
#   }
#   if (size_factor_normalize) {
#     Biobase::exprs(cicero_cds) <- t(t(Biobase::exprs(cicero_cds))/Biobase::pData(cicero_cds)$Size_Factor)
#   }
#   return(list(cicero_cds = cicero_cds, nearest_neighour_mat = cell_sample))
# }


