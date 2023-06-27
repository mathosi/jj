# jj

jj is an R package with functions revolving around single-cell analysis. 
It allows an easy generation of scatterplots, heatmaps and other visualizations based on the packages ggplot2 or ComplexHeatmap.
The most important functions are:

- `jj_plot_features`: plot a 2D scatterplot of a dimensionality reduction (e.g. UMAP). Accepts a Seurat object, ArchR object or data.frame as input. Can plot a feature from any of the available data matrices in the object such as the meta data, RNA count matrix, TF deviations scores etc.
- `j`: Get a quick overview on any R object. Prints the object's dimensions and its first few entries 
- `jj_plot_sparse_by_group`: Plot a sparse numeric feature summarized by group as violins or boxplots
- `jj_plot_numeric_by_group`: Plot a numeric feature summarized by group as violins or boxplots
- `jj_plot_categorical_by_group`: Plot a categorical feature summarized by group as violins or boxplots
- `jj_plot_heatmap` / `jj_plot_dotplot` : Plot a heatmap/dotplot for a selection of genes by group
- `jj_plot_hist` : Plot a histogram and optionally count values above and below a threshold
- `jj_plot_volcano` : Plot a volcano plot given a data.frame with fold change and p-value column
- `jj_plot_fc_fc` : Plot two fold change columns against each other to compare differential analysis results
