library(here)
library(ggplot2)
library(TMExplorer)
library(Seurat)
here::i_am("code/lung_epithelium.R")
res = queryTME(geo_accession = "GSE131907")
res.seurat=as.Seurat(res[[1]], data = NULL)
rm(res)
nLung <- res.seurat[,grepl('LUNG_N',names(Idents(res.seurat)))]
tLung <- res.seurat[,grepl('LUNG_T',names(Idents(res.seurat)))]
nLung$origin <- "Normal lung"
tLung$origin <- "Early stage tumor lung"
lung <- merge(nLung, tLung)
rm(nLung,tLung)
lungEpi <- lung[,grepl("Epithelial",lung$label,fixed = TRUE)]
lungEpiv5 <- as(lungEpi[["originalexp"]], Class = "Assay5")
lungEpi5 <- CreateSeuratObject(lungEpiv5, meta.data=lungEpi@meta.data)
lst = as.list(as.data.frame(matrix(unlist(strsplit(names(Idents(lungEpi5)),"_")),ncol=3,byrow=TRUE)[,2:3]))
lungEpi5$experiment <- paste(lst[[1]],lst[[2]],sep ="_")
lungEpi5[["RNA"]] <- split(lungEpi5[["RNA"]], f = lungEpi5$experiment)
lungEpi5 <- NormalizeData(lungEpi5)
lungEpi5 <- FindVariableFeatures(lungEpi5)
lungEpi5 <- ScaleData(lungEpi5)
lungEpi5 <- RunPCA(lungEpi5)
lungEpi5 <- IntegrateLayers(lungEpi5, method = RPCAIntegration, orig.reduction = "pca", new.reduction = "integrated.rpca")
lungEpi5[["RNA"]] <- JoinLayers(lungEpi5[["RNA"]])
lungEpi5 <- RunUMAP(lungEpi5, dims = 1:30, reduction = "integrated.rpca")
DimPlot(lungEpi5, reduction = "umap", group.by=c("origin","label"))
