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
  dp <- DimPlot(lymphFB, reduction = "umap", group.by=c("origin","label"), combine = FALSE)
  dp[[1]] <- dp[[1]] + ggtitle(sprintf("%s grouped by disease state", fn))
  dp[[2]] <- dp[[2]] + ggtitle(sprintf("%s grouped by cell subtype", fn))
  dp <- dp[[1]] + dp[[2]]
  ggsave(here("figures", sprintf("LN %s.png", gsub("/","+",fn))), plot = dp, width = 13, height = 3.5)
  return(dp)
}
DEget <- function(lymphFB, ident1, ident2) {
  bulk <- AggregateExpression(lymphFB, return.seurat = TRUE, assays = "RNA", group.by = c("label","experiment","origin"))
  Idents(bulk) <- "origin"
  de_markers = FindMarkers(bulk, ident.1 = ident1, ident.2 = ident2, slot = "counts", test.use = "wilcox")
  de_markers$gene <- rownames(de_markers)
  return(de_markers)
}
DEanalyze <- function(de_markers, fn) {
  dep <- ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
    ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,
                                                                            "")), colour = "red", size = 3)
  significant <- de_markers[de_markers$p_val_adj < 0.01,]
  upregulated <- significant[significant$pct.2 > significant$pct.1,]
  downregulated <- significant[significant$pct.2 < significant$pct.1,]
  ggsave(here("figures", sprintf("DE lymph node %s.png", gsub("/","+",fn))), plot = dep, width = 6, height = 4)
  write.csv(upregulated, here("figures", sprintf("%s upregulated in tumor LN.csv",gsub("/","+",fn))))
  write.csv(downregulated, here("figures", sprintf("%s downregulated in tumor LN.csv",gsub("/","+",fn))))
  pantherup <- rba_panther_enrich(genes = rownames(upregulated), organism = 9606, annot_dataset = "GO:0008150")
  pantherdown <- rba_panther_enrich(genes = rownames(downregulated), organism = 9606, annot_dataset = "GO:0008150")
  write.csv(pantherup$result, here("figures", sprintf("%s upregulated pathways in tumor LN.csv",gsub("/","+",fn))))
  write.csv(downregulated$result, here("figures", sprintf("%s downregulated pathways in tumor LN.csv",gsub("/","+",fn))))
  return(dep)
}