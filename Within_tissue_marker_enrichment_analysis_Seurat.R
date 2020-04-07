#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Find markers of each cluster and perform enrichment analysis within each tissue
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_by_organ/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/by_organ")
  place <- "wolfpack"
}

# Libraries
library('Seurat')
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')

# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp_sce <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/practice_all_data/Prefiltered_QC_experiment_practice.rds") # uses practice data if local
} else {
  filtered_exp_sce <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/Prefiltered_QC_experiment_all.rds") # uses whole dataset if wolfpack
}

for(tissue in c("Liver", "LN", "Lung", "Primary")){
# Load prefiltered and clustered Seurat Object
if(place == "local") {
  tissue_exp <- readRDS(paste0(tissue, "/Prefiltered_experiment_practice_seurat_integrated_", tissue,".rds")) # uses practice data if local
} else {
  tissue_exp <- readRDS(paste0(tissue, "/Prefiltered_experiment_all_seurat_integrated_", tissue,".rds")) # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Find markers for all clusters and perform GO enrichment
# --------------------------------------------------------------------------

all_markers <- FindAllMarkers(tissue_exp, assay = "RNA")
tissue_exp@misc$all_markers <- all_markers

# Add Entrez IDs for KEGG and REACTOME analyses
rowData(filtered_exp_sce)$EntrezID <- mapIds(org.Hs.eg.db, keys=rowData(filtered_exp_sce)$Ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
rowMetaData <- data.frame(rowData(filtered_exp_sce))
tissue_exp@assays$RNA@meta.features <- merge(tissue_exp@assays$RNA@meta.features, rowMetaData, by.x = 0, by.y = 0)

# Make gene universe(s)
universe <- as.character(unique(tissue_exp@assays$RNA@meta.features[tissue_exp@assays$RNA@meta.features$Any_Active, "Ensembl"]))
universe_entrez <- as.character(unique(tissue_exp@assays$RNA@meta.features[tissue_exp@assays$RNA@meta.features$Any_Active, "EntrezID"]))

# Add ENSEMBL gene ids to Marker dataframe
idx <- match(all_markers$gene, tissue_exp@assays$RNA@meta.features$Row.names)
all_markers$Ensembl <- tissue_exp@assays$RNA@meta.features$Ensembl [idx]
all_markers$EntrezID <- tissue_exp@assays$RNA@meta.features$EntrezID [idx]

# Split into clusters and find GO annotations all, up or down regulated
all_markers$Upregulated <- all_markers$avg_logFC > 0
all_markers$Downregulated <- all_markers$avg_logFC < 0
markers.tissue_exp <- split(all_markers, all_markers$cluster)

# All-regulated markers in each cluster
for(i in names(markers.tissue_exp)){
  dir.create(paste0(tissue, "/markers"))
  dir.create(paste0(tissue, "/markers/Cluster_", i))
  interesting <- markers.tissue_exp[[i]]
  write.csv(interesting, paste0(tissue, "/markers/Cluster_", i,"/Cluster_", i, "_allregulated_markers.csv"))
  
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
  write.csv(GOBP_GOI, paste0(tissue, "/markers/Cluster_", i,"/GOBP_markers_all_Cluster_", i, ".csv"))
  
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
  write.csv(GOMF_GOI, paste0(tissue, "/markers/Cluster_", i,"/GOMF_markers_all_Cluster_", i, ".csv"))
  
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
  write.csv(GOCC_GOI, paste0(tissue, "/markers/Cluster_", i,"/GOCC_markers_all_Cluster_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0(tissue, "/markers/Cluster_", i,"/KEGG_markers_all_Cluster_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0(tissue, "/markers/Cluster_", i,"/Reactome_markers_all_Cluster_", i, ".csv"))
}

# Up-regulated markers in each cluster
for(i in names(markers.tissue_exp)){
  interesting <- markers.tissue_exp[[i]]
  interesting <- interesting[interesting$Upregulated, ]
  write.csv(interesting, paste0(tissue, "/markers/Cluster_", i,"/Cluster_", i, "_upregulated_markers.csv"))
  
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
  write.csv(GOBP_GOI, paste0(tissue, "/markers/Cluster_", i,"/GOBP_markers_up_Cluster_", i, ".csv"))
  
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
  write.csv(GOMF_GOI, paste0(tissue, "/markers/Cluster_", i,"/GOMF_markers_up_Cluster_", i, ".csv"))
  
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
  write.csv(GOCC_GOI, paste0(tissue, "/markers/Cluster_", i,"/GOCC_markers_up_Cluster_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0(tissue, "/markers/Cluster_", i,"/KEGG_markers_up_Cluster_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0(tissue, "/markers/Cluster_", i,"/Reactome_markers_up_Cluster_", i, ".csv"))
}

# Down-regulated markers in each cluster
for(i in names(markers.tissue_exp)){
  interesting <- markers.tissue_exp[[i]]
  interesting <- interesting[interesting$Downregulated, ]
  write.csv(interesting, paste0(tissue, "/markers/Cluster_", i,"/Cluster_", i, "_downregulated_markers.csv"))
  
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
  write.csv(GOBP_GOI, paste0(tissue, "/markers/Cluster_", i,"/GOBP_markers_down_Cluster_", i, ".csv"))
  
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
  write.csv(GOMF_GOI, paste0(tissue, "/markers/Cluster_", i,"/GOMF_markers_down_Cluster_", i, ".csv"))
  
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
  write.csv(GOCC_GOI, paste0(tissue, "/markers/Cluster_", i,"/GOCC_markers_down_Cluster_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0(tissue, "/markers/Cluster_", i,"/KEGG_markers_down_Cluster_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p_val_adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0(tissue, "/markers/Cluster_", i,"/Reactome_markers_down_Cluster_", i, ".csv"))
}

if(place == "local") {
  saveRDS(tissue_exp, paste0(tissue, "/Prefiltered_experiment_practice_seurat_integrated_markers_", tissue,".rds"))
} else {
  saveRDS(tissue_exp, paste0(tissue, "/Prefiltered_experiment_all_seurat_integrated_markers_", tissue,".rds"))
}
}

