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
res.seurat[,grepl('.+LUNG_N.+',names(Idents(res.seurat)))]$origin <- "nLung"
nLung$origin <- "Normal lung"
tLung$origin <- "Early stage tumor lung"
lung <- merge(nLung, tLung)
rm(nLung,tLung)
lungEpi <- lung[,grepl("Epithelial",lung$label,fixed = TRUE)]
lungEpi <- NormalizeData(lungEpi)
lungEpi <- FindVariableFeatures(lungEpi)
lungEpi <- ScaleData(lungEpi)
lungEpi <- RunPAC(lungEpi)
lungEpi <- RunPCA(lungEpi)
DimPlot(lungEpi, reduction = "pca", group.by=c("origin","label"))