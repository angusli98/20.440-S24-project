PCAIntegrate <- function(resv5,sample,celltype) {
  lymph <- resv5[,grepl(sample,resv5$origin,fixed=TRUE)]
  lymphFB <- lymph[,grepl(celltype,lymph$label,fixed = TRUE)]
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
  return(lymphFB)
}
UMAPDp <- function(lymphFB, fn) {
  dp <- DimPlot(lymphFB, reduction = "umap", group.by=c("origin","label"))
  ggsave(here("figures", fn), plot = dp)
  return(dp)
}
DEget <- function(lymphFB, ident1, ident2) {
  bulk <- AggregateExpression(lymphFB, return.seurat = TRUE, assays = "RNA", group.by = c("label","experiment","origin"))
  Idents(bulk) <- "origin"
  de_markers = FindMarkers(bulk, ident.1 = ident1, ident.2 = ident2, slot = "counts", test.use = "wilcox")
  de_markers$gene <- rownames(de_markers)
  return(de_markers)
}
DEplot <- function(de_markers, fn) {
  dep <- ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
    ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,
                                                                            "")), colour = "red", size = 3)
  ggsave(here("figures", sprintf("%s.png",fn)), plot = dep)
  write.csv(de_markers, here("figures", sprintf("%s.csv",fn)))
  return(dep)
}