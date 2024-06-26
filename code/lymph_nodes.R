library(here)
library(ggplot2)
library(ggrepel)
library(TMExplorer)
library(Seurat)
library(rbioapi)
source("processing_functions.R")
here::i_am("code/lymph_nodes.R")
res = queryTME(geo_accession = "GSE131907") # Download dataset
res.seurat=as.Seurat(res[[1]], data = NULL) # Conversion to Seurat
rm(res)
res.seurat$origin <- "none"
res.seurat$origin[grepl('LUNG_N',names(Idents(res.seurat)))] <- "Normal lung" # Labeling samples by tissue type of origin
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
RNAv5 <- as(res.seurat[["originalexp"]], Class = "Assay5") # Conversion to Seurat V5
resv5 <- CreateSeuratObject(RNAv5, meta.data=res.seurat@meta.data)
rm(res.seurat,RNAv5)
for (tissue in c("Myeloid cells", "T/NK cells", "Epithelial cells", "Endothelial cells", "Fibroblasts", "B lymphocytes")) {
  lymphFB <- PCAIntegrate(resv5, "lymph", tissue) # PCA, integration, and UMAP pipeline of specified tissue and cell type
  UMAPDp(lymphFB, tissue) # Plots UMAP
  de_markers <- DEget(lymphFB, "Metastatic lymph node", "Normal lymph node") # Pseudobulk analysis for differentially expressed genes between conditions
  DEanalyze(de_markers, tissue) # Volcano plot, extraction of significantly changed genes, GO analysis
}
