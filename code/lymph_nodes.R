library(here)
library(ggplot2)
library(ggrepel)
library(TMExplorer)
library(Seurat)
here::i_am("code/lymph_nodes.R")
res = queryTME(geo_accession = "GSE131907")
res.seurat=as.Seurat(res[[1]], data = NULL)
rm(res)
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
lymph = subset(resv5, subset = (origin == "Metastatic lymph node" | origin == "Normal lymph node"))
lymphFB <- lymph[,grepl("Fibroblasts",lymph$label,fixed = TRUE)]
lst = as.list(as.data.frame(matrix(unlist(strsplit(names(Idents(lymphFB)),"_")),ncol=3,byrow=TRUE)[,2:3]))
lymphFB$experiment <- paste(lst[[1]],lst[[2]],sep ="_")
lymphFB[["RNA"]] <- split(lymphFB[["RNA"]], f = lymphFB$experiment)
lymphFB <- NormalizeData(lymphFB)
lymphFB <- FindVariableFeatures(lymphFB)
lymphFB <- ScaleData(lymphFB)
lymphFB <- RunPCA(lymphFB)
lymphFB <- IntegrateLayers(lymphFB, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony")
lymphFB[["RNA"]] <- JoinLayers(lymphFB[["RNA"]])
lymphFB <- RunUMAP(lymphFB, dims = 1:30, reduction = "integrated.harmony")
DimPlot(lymphFB, reduction = "umap", group.by=c("origin","label"))
dp <- DimPlot(lymphFB, reduction = "umap", group.by=c("origin","label"))
ggsave(here("figures", "LN fibroblasts.png"), plot = dp)
bulk <- AggregateExpression(lymphFB, return.seurat = TRUE, assays = "RNA", group.by = c("label","experiment","origin"))
Idents(bulk) <- "origin"
de_markers = FindMarkers(bulk, ident.1 = "Normal lymph node", ident.2 = "Metastatic lymph node", slot = "counts", test.use = "wilcox")
ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
  ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,
                                                                          "")), colour = "red", size = 3)