library(here)
library(ggplot2)
library(ggrepel)
library(TMExplorer)
library(Seurat)
source("processing_functions.R")
here::i_am("code/lymph_nodes.R")
res = queryTME(geo_accession = "GSE131907")
res.seurat=as.Seurat(res[[1]], data = NULL)
rm(res)
res.seurat$origin <- "none"
res.seurat$origin[grepl('LUNG_N',names(Idents(res.seurat)))] <- "Normal lung"
res.seurat$origin[grepl('LUNG_T',names(Idents(res.seurat)))] <- "Early stage lung tumor"
res.seurat$origin[grepl('LN_',names(Idents(res.seurat)))] <- "Normal lymph node"
res.seurat$origin[grepl('NS_',names(Idents(res.seurat)))] <- "Brain metastasis"
res.seurat$origin[grepl('EFFUSION_',names(Idents(res.seurat)))] <- "Pleural effusion"
res.seurat$origin[grepl('BRONCHO_58',names(Idents(res.seurat)))] <- "Advanced stage lung tumor"
res.seurat$origin[grepl('BRONCHO_11',names(Idents(res.seurat)))] <- "Metastatic lymph node"
res.seurat$origin[grepl('EBUS',names(Idents(res.seurat)))] <- "Metastatic lymph node"
res.seurat$origin[grepl('EBUS_06',names(Idents(res.seurat)))] <- "Advanced stage lung tumor"
res.seurat$origin[grepl('EBUS_28',names(Idents(res.seurat)))] <- "Advanced stage lung tumor"
res.seurat$origin[grepl('EBUS_49',names(Idents(res.seurat)))] <- "Advanced stage lung tumor"
RNAv5 <- as(res.seurat[["originalexp"]], Class = "Assay5")
resv5 <- CreateSeuratObject(RNAv5, meta.data=res.seurat@meta.data)
rm(res.seurat,RNAv5)
lymphFB <- PCAIntegrate(resv5, "lymph", "B lymphocytes")
UMAPDp(lymphFB, "LN B.png")
de_markers <- DEget(lymphFB, "Normal lymph node", "Metastatic lymph node")
DEplot(de_markers, "DE lymph node B")
