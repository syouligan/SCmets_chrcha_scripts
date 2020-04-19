#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells using Seurat.
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data")
  place <- "wolfpack"
}

# Libraries
library('Seurat')
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('readr')
library('Matrix')
library('phateR')
library('cowplot')
library('factoextra')
library('gtools')
library('doParallel')
library('foreach')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')

# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp_sce <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/practice_all_data/Prefiltered_QC_experiment_practice.rds") # uses practice data if local
} else {
  filtered_exp_sce <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/Prefiltered_QC_experiment_all.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Set up prefiltered object, add cell cycle difference info and split by sample
# --------------------------------------------------------------------------
filtered_exp <- as.Seurat(filtered_exp_sce, counts = "counts", data = NULL) # convert to Seurat

# Add Entrez IDs for KEGG and REACTOME analyses
rowData(filtered_exp_sce)$EntrezID <- mapIds(org.Hs.eg.db, keys=rowData(filtered_exp_sce)$Ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
rowMetaData <- data.frame(rowData(filtered_exp_sce))
filtered_exp@assays$RNA@meta.features <- merge(filtered_exp@assays$RNA@meta.features, rowMetaData, by.x = 0, by.y = 0)

# Add cell cycle and stress scores
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
filtered_exp <- CellCycleScoring(filtered_exp, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
filtered_exp$CC.Difference <- filtered_exp$S.Score - filtered_exp$G2M.Score
digest_stress <- list(c("FOS", "CXCL2", "ZFP36", "FOSB", "DUSP1", "ATF3", "CXCL8", "NR4A1", "CXCL3", "PPP1R15A", "JUNB", "EGR1", "HSPA1A", "HSPA1B", "SOCS3", "KLF6", "JUN", "IER2", "CXCL1", "NKFBIA", "HSPA6", "DNAJB1", "IER3", "CCNL1", "MTRNR2L2", "IER5", "ID1", "CEBPD", "KRT6A", "CYR61", "DEPP1", "CLDN4", "IRF1", "DUSP2", "BTG2", "PLAUR", "MAFF", "KLF4", "PHLDA2", "TNFAIP3"))
filtered_exp <- AddModuleScore(object = filtered_exp, features = digest_stress, name = 'digest_stress')
Idents(object = filtered_exp) <- "Tissue"

# Find markers for 1:1 comparison between tissues
# --------------------------------------------------------------------------

# Make table of comparisons across tissue
comparisons <- combinations(4, 2, unique(as.character(Idents(object = filtered_exp))))
colnames(comparisons) <- c("Tissue_1", "Tissue_2")
rownames(comparisons) <- 1:nrow(comparisons)

cores <- detectCores()
cl <- makeCluster(round(0.5*cores[1]))
registerDoParallel(cl)

# Find markers for each comparison
markers.filtered_exp <- foreach(x = 1:nrow(comparisons), .packages='Seurat') %dopar% {  

# Make vectors for each tissue in the comparison
tissue1 <- comparisons[x, 1]
tissue2 <- comparisons[x, 2]

# Find markers for each comparison
features <- as.character(unique(filtered_exp@assays$RNA@meta.features[filtered_exp@assays$RNA@meta.features$Any_Active, "Row.names"]))
markers <- FindMarkers(filtered_exp, ident.1 = tissue1, ident.2 = tissue2, assay = "RNA", test.use = "LR", latent.vars = c("Replicate", "Lib_size", "S.Score", "G2M.Score", "digest_stress1"), features = features)
markers$Upregulated <- markers$avg_logFC > 0
markers$Downregulated <- markers$avg_logFC < 0

# Add ENSEMBL gene ids to Marker dataframe
idx <- match(rownames(markers), filtered_exp@assays$RNA@meta.features$Row.names)
markers$Ensembl <- filtered_exp@assays$RNA@meta.features$Ensembl [idx]
markers$EntrezID <- filtered_exp@assays$RNA@meta.features$EntrezID [idx]

return(markers)
}

names(markers.filtered_exp) <- paste0(comparisons[, 1], "_", comparisons[, 2])
filtered_exp@misc$tissue_markers <- markers.filtered_exp

stopCluster(cl)

# Perfrom GO analysis of each marker list
# --------------------------------------------------------------------------

# Make gene universe(s)
universe <- as.character(unique(filtered_exp@assays$RNA@meta.features[filtered_exp@assays$RNA@meta.features$Any_Active, "Ensembl"]))
universe_entrez <- as.character(unique(filtered_exp@assays$RNA@meta.features[filtered_exp@assays$RNA@meta.features$Any_Active, "EntrezID"]))

# All-regulated markers in each cluster
for(i in names(markers.filtered_exp)){
  dir.create(paste0("markers/Tissue_", i))
  interesting <- markers.filtered_exp[[i]]
  write.csv(interesting, paste0("markers/Tissue_", i,"/Tissue_", i, "_allregulated_markers.csv"))
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/Tissue_", i,"/GOBP_markers_all_Tissue_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/Tissue_", i,"/GOMF_markers_all_Tissue_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/Tissue_", i,"/GOCC_markers_all_Tissue_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/Tissue_", i,"/KEGG_markers_all_Tissue_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/Tissue_", i,"/Reactome_markers_all_Tissue_", i, ".csv"))
}

# Up-regulated markers in each cluster
for(i in names(markers.filtered_exp)){
  interesting <- markers.filtered_exp[[i]]
  interesting <- interesting[interesting$Upregulated, ]
  write.csv(interesting, paste0("markers/Tissue_", i,"/Tissue_", i, "_upregulated_markers.csv"))
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/Tissue_", i,"/GOBP_markers_up_Tissue_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/Tissue_", i,"/GOMF_markers_up_Tissue_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/Tissue_", i,"/GOCC_markers_up_Tissue_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/Tissue_", i,"/KEGG_markers_up_Tissue_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/Tissue_", i,"/Reactome_markers_up_Tissue_", i, ".csv"))
}

# Down-regulated markers in each cluster
for(i in names(markers.filtered_exp)){
  interesting <- markers.filtered_exp[[i]]
  interesting <- interesting[interesting$Downregulated, ]
  write.csv(interesting, paste0("markers/Tissue_", i,"/Tissue_", i, "_downregulated_markers.csv"))
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/Tissue_", i,"/GOBP_markers_down_Tissue_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/Tissue_", i,"/GOMF_markers_down_Tissue_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting[interesting$p_val_adj < 0.05, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/Tissue_", i,"/GOCC_markers_down_Tissue_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/Tissue_", i,"/KEGG_markers_down_Tissue_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/Tissue_", i,"/Reactome_markers_down_Tissue_", i, ".csv"))
}

if(place == "local") {
  saveRDS(filtered_exp, "Prefiltered_experiment_practice_seurat_tissueDGE.rds")
} else {
  saveRDS(filtered_exp, "Prefiltered_experiment_all_seurat_tissueDGE.rds")
}