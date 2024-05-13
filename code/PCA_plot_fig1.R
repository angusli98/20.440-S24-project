## The analyses performed in this script were not used in the final project.
## It is included to provide context for development of the codebase.
library(here)
library(ggplot2)
library(TMExplorer)
library(Seurat)
here::i_am("code/PCA_plot_fig1.R")
res = queryTME(geo_accession = "GSE131907") # Download dataset
res.seurat=as.Seurat(res[[1]], data = NULL) # Convert to Seurat
res.seurat <- NormalizeData(res.seurat) # Normalize data by total expression in cell
res.seurat <- FindVariableFeatures(res.seurat) # Take 2000 most variable features
res.seurat <- ScaleData(res.seurat) # Scale and center data
res.seurat <- RunPCA(res.seurat) # PCA on variable features
dp <- DimPlot(res.seurat, reduction = "pca") + NoLegend() # Plot of first two PCA dimensions
ggsave(here("figures", "Figure1.png"), plot = dp) # Export figure to figures folder