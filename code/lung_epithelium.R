## The analyses performed in this script were not used in the final project.
## It is included to provide context for development of the codebase.
library(here)
library(ggplot2)
library(TMExplorer)
library(Seurat)
here::i_am("code/lung_epithelium.R")
res = queryTME(geo_accession = "GSE131907") # Download dataset
res.seurat=as.Seurat(res[[1]], data = NULL) # Convert to Seurat
rm(res)
nLung <- res.seurat[,grepl('LUNG_N',names(Idents(res.seurat)))] # Extract normal lung samples
tLung <- res.seurat[,grepl('LUNG_T',names(Idents(res.seurat)))] # Extract tumor lung samples
nLung$origin <- "Normal lung" # Assign label for tissue type of origin
tLung$origin <- "Early stage tumor lung"
lung <- merge(nLung, tLung) # Create merged dataset of samples
rm(nLung,tLung)
lungEpi <- lung[,grepl("Epithelial",lung$label,fixed = TRUE)] # Extract epithelial cells
lungEpiv5 <- as(lungEpi[["originalexp"]], Class = "Assay5") # Conversion to Seurat v5
lungEpi5 <- CreateSeuratObject(lungEpiv5, meta.data=lungEpi@meta.data)
lst = as.list(as.data.frame(matrix(unlist(strsplit(names(Idents(lungEpi5)),"_")),ncol=3,byrow=TRUE)[,2:3])) # Extract patient of origin from sample name
lungEpi5$experiment <- paste(lst[[1]],lst[[2]],sep ="_") # Create new metadata field for patient of origin
lungEpi5[["RNA"]] <- split(lungEpi5[["RNA"]], f = lungEpi5$experiment) # Separate dataset into layers by patient (i.e. sample)
lungEpi5 <- NormalizeData(lungEpi5) # Normalization by expression per cell
lungEpi5 <- FindVariableFeatures(lungEpi5) # Finding 2000 most variable features
lungEpi5 <- ScaleData(lungEpi5) # Scaling and centering data
lungEpi5 <- RunPCA(lungEpi5) # PCA on variable features
lungEpi5 <- IntegrateLayers(lungEpi5, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca") # Integration of data from different patients
lungEpi5[["RNA"]] <- JoinLayers(lungEpi5[["RNA"]]) # Reunification into single dataset
lungEpi5 <- RunUMAP(lungEpi5, dims = 1:30, reduction = "integrated.rpca") # UMAP on integrated data
DimPlot(lungEpi5, reduction = "umap", group.by=c("origin","label")) # UMAP plot
