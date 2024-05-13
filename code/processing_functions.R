# This file contains functions used in the lymph node analysis pipeline.
PCAIntegrate <- function(resv5,sample,celltype) { # PCA, integration, and UMAP pipeline of specified tissue and cell type
  lymph <- resv5[,grepl(sample,resv5$origin,fixed=TRUE)] # Extraction of samples from desired tissue (e.g. lymph node)
  lymphFB <- lymph[,grepl(celltype,lymph$label,fixed = TRUE)] # Extraction of desired cell type from sample
  lst = as.list(as.data.frame(matrix(unlist(strsplit(names(Idents(lymphFB)),"_")),ncol=3,byrow=TRUE)[,2:3])) # Extract patient of origin from sample name
  lymphFB$experiment <- paste(lst[[1]],lst[[2]],sep ="_") # Create new metadata field for patient of origin
  lymphFB[["RNA"]] <- split(lymphFB[["RNA"]], f = lymphFB$experiment) # Separate dataset into layers by patient (i.e. sample)
  lymphFB <- NormalizeData(lymphFB) # Normalization by expression per cell
  lymphFB <- FindVariableFeatures(lymphFB) # Finding 2000 most variable features
  lymphFB <- ScaleData(lymphFB) # Scaling and centering data
  lymphFB <- RunPCA(lymphFB) # PCA on variable features
  lymphFB <- IntegrateLayers(lymphFB, method = HarmonyIntegration, orig.reduction = "pca", new.reduction = "integrated.harmony") # Integration of data from different patients
  lymphFB[["RNA"]] <- JoinLayers(lymphFB[["RNA"]]) # Reunification into single dataset
  lymphFB <- RunUMAP(lymphFB, dims = 1:30, reduction = "integrated.harmony") # UMAP on integrated data
  return(lymphFB)
}
UMAPDp <- function(lymphFB, fn) { # Plots UMAP
  dp <- DimPlot(lymphFB, reduction = "umap", group.by=c("origin","label"), combine = FALSE) # Plots UMAP colored by disease state and cell subtype
  dp[[1]] <- dp[[1]] + ggtitle(sprintf("%s grouped by disease state", fn))
  dp[[2]] <- dp[[2]] + ggtitle(sprintf("%s grouped by cell subtype", fn))
  dp <- dp[[1]] + dp[[2]]
  ggsave(here("figures", sprintf("LN %s.png", gsub("/","+",fn))), plot = dp, width = 13, height = 3.5) # Saves UMAP plot to figures folder
  return(dp)
}
DEget <- function(lymphFB, ident1, ident2) { # Pseudobulk analysis for differentially expressed genes between conditions
  bulk <- AggregateExpression(lymphFB, return.seurat = TRUE, assays = "RNA", group.by = c("label","experiment","origin")) # Pseudobulk data
  Idents(bulk) <- "origin" # Specify to perform analysis by disease state
  de_markers = FindMarkers(bulk, ident.1 = ident1, ident.2 = ident2, slot = "counts", test.use = "wilcox") # Find differentially expressed genes by Wilcoxon test
  de_markers$gene <- rownames(de_markers) # Assign gene names as row labels
  return(de_markers)
}
DEanalyze <- function(de_markers, fn) { # Volcano plot, extraction of significantly changed genes, GO analysis
  dep <- ggplot(de_markers, aes(avg_log2FC, -log10(p_val))) + geom_point(size = 0.5, alpha = 0.5) + theme_bw() +
    ylab("-log10(unadjusted p-value)") + geom_text_repel(aes(label = ifelse(p_val_adj < 0.01, gene,"")),
      colour = "red", size = 3) + ggtitle(sprintf("Differentially expressed genes in metastatic lymph node %s", fn)) # Volcano plot
  significant <- de_markers[de_markers$p_val_adj < 0.01,] # # Extract significantly changed genes
  upregulated <- significant[significant$pct.1 > significant$pct.2,] # Extract significantly upregulated genes in diseased state
  downregulated <- significant[significant$pct.1 < significant$pct.2,] # Extract significantly downregulated genes in diseased state
  ggsave(here("figures", sprintf("DE lymph node %s.png", gsub("/","+",fn))), plot = dep, width = 7.5, height = 5) # Save volcano plot to figures folder
  write.csv(upregulated, here("figures", sprintf("%s upregulated in tumor LN.csv",gsub("/","+",fn)))) # Save upregulated genes to file
  write.csv(downregulated, here("figures", sprintf("%s downregulated in tumor LN.csv",gsub("/","+",fn)))) # Save downregulated genes to file
  pantherup <- rba_panther_enrich(genes = rownames(upregulated), organism = 9606, annot_dataset = "GO:0008150") # Gene Ontology enrichment analysis by PANTHER web service on upregulated genes, using Homo sapiens biological process annotations
  pantherdown <- rba_panther_enrich(genes = rownames(downregulated), organism = 9606, annot_dataset = "GO:0008150") # Gene Ontology enrichment analysis by PANTHER web service on downregulated genes, using Homo sapiens biological process annotations
  write.csv(pantherup$result, here("figures", sprintf("%s upregulated pathways in tumor LN.csv",gsub("/","+",fn)))) # Save enriched process annotations among upregulated genes to file
  write.csv(downregulated$result, here("figures", sprintf("%s downregulated pathways in tumor LN.csv",gsub("/","+",fn)))) # Save enriched process annotations among downregulated genes to file
  return(dep)
}