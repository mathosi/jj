% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jj_plot_feature_by_group.R
\name{plot_feature_by_group}
\alias{plot_feature_by_group}
\alias{jj_plot_sparse_by_group_seurat}
\alias{jj_plot_sparse_by_group}
\alias{jj_plot_numeric_by_group}
\alias{jj_plot_categorical_by_group}
\title{violin plot of feature by group}
\usage{
jj_plot_sparse_by_group_seurat(
  seurat_obj,
  gene_plot,
  group_column,
  assay = NULL,
  slot = "data",
  ...
)

jj_plot_sparse_by_group(
  rna_mat,
  gene_plot,
  group_vec,
  x_lab = "Group",
  theme_use = theme_minimal(),
  plot_cell_sample = FALSE,
  plot_zero_fraction = TRUE,
  plot_mean = TRUE,
  plot_group_size = FALSE,
  type = "violin",
  custom_colors = NULL,
  order = FALSE,
  flip_coordinates = FALSE,
  cap_top = NULL,
  ...
)

jj_plot_numeric_by_group(
  df,
  feature_column,
  group_column,
  custom_colors = NULL,
  plot_cell_sample = FALSE,
  plot_mean = TRUE,
  show_annotation = NULL,
  annotation_at = NULL,
  theme_use = theme_minimal(),
  type = "violin",
  order = FALSE,
  flip_coordinates = FALSE,
  ...
)

jj_plot_categorical_by_group(
  df,
  feature_column,
  group_column,
  custom_colors = NULL,
  absolute_numbers = FALSE,
  return_df = FALSE,
  flip_coordinates = FALSE,
  add_text = FALSE,
  text_size = 5,
  theme_use = theme_minimal(),
  ...
)
}
\arguments{
\item{seurat_obj}{Seurat object that contains the \code{assay} with the feature \code{gene_plot}}

\item{gene_plot}{Gene from the RNA matrix to plot}

\item{group_column}{string of the column with the variable used to group the feature}

\item{...}{further arguments passed to the main geom of ggplot}

\item{group_vec}{Vector of strings with the grouping information}

\item{x_lab}{label for the groups}

\item{theme_use}{theme to use, default: theme_minimal()}

\item{plot_cell_sample}{if TRUE, plot a sample of cells for each group (equal number)}

\item{plot_zero_fraction}{for sparse data, plot the fraction of zero counts per group as pie}

\item{plot_mean}{Plot the mean value per group as horizontal line}

\item{plot_group_size}{Plot number of cells per group, in jj_plot_numeric_by_group replaced by show_annotation with choices c('n','mean','sd')}

\item{type}{Type of plot, options are 'violin' or 'boxplot'}

\item{custom_colors}{named vector of colors to use to fill the violins/boxplots}

\item{order}{if TRUE, order the groups  by their mean value}

\item{flip_coordinates}{flip coordinate system}

\item{cap_top}{cap the values at a quantile or fixed value}

\item{df}{data.frame containing the columns passed in feature_column and group_column}

\item{feature_column}{string of the column with feature that should be quantified}

\item{absolute_numbers}{if TRUE, plot absolute counts per category instead of relative fractions per group}

\item{return_df}{if TRUE, instead of plotting, return the data.frame with the data}

\item{text_size}{size of numbers in plots}
}
\description{
Plot a gene/numeric/categorical feature by group as violin or boxplot
}
\examples{
#plot as boxplot with additional mean, number of cells per group and cell sample (requires ggbeeswarm)
jj_plot_numeric_by_group(pbmc_small@meta.data, feature_column = 'nFeature_RNA', group_column = 'groups',
                         plot_mean = T, plot_cell_sample = T, type = 'boxplot')
#plot as violin with custom colours
jj_plot_numeric_by_group(pbmc_small@meta.data, feature_column = 'nFeature_RNA', group_column = 'groups',
                         custom_colors = c(g1='green', g2='blue'), type = 'violin')
jj_plot_numeric_by_group(pbmc_small@meta.data, feature_column = 'nFeature_RNA', group_column = 'groups',
                         show_annotation = c('n', 'mean', 'sd'), annotation_at = 100, type = 'violin')
#plot a sparse feature directly from Seurat
jj_plot_sparse_by_group_seurat(pbmc_small, 'CD79A', 'groups', assay='RNA', slot='data')
#or from a sparse matrix
sp_mat = GetAssayData(pbmc_small)
jj_plot_sparse_by_group(sp_mat, gene_plot = 'MS4A1', group_vec = pbmc_small$groups, order=T)
#barplot of fractions by group, eg fractions of cluster annotations per group
jj_plot_categorical_by_group(pbmc_small[[]], feature_column = 'RNA_snn_res.1', group_column =  'groups')
#or using absolute counts
jj_plot_categorical_by_group(pbmc_small[[]], feature_column = 'RNA_snn_res.1', group_column =  'groups', absolute_numbers = T)
}
