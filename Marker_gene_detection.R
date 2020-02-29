#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Find and save markers for each cluster
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data")
  place <- "wolfpack"
}

# Libraries
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('readr')
library('Matrix')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')

set.seed(100)
# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_experiment_practice_merge_cluster.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_experiment_all_merge_cluster.rds") # uses whole dataset if wolfpack
}

rowData(filtered_exp)$EntrezID <- mapIds(org.Hs.eg.db, keys=rowData(filtered_exp)$Ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Find conservative markers between each cluster and save a csv
markers.filtered_exp <- findMarkers(filtered_exp, test="wilcox", filtered_exp$cluster, direction="up", block = filtered_exp$Sample, pval.type = "all")
for(i in names(markers.filtered_exp)){
  dir.create(paste0("markers/Cluster_", i))
  interesting <- markers.filtered_exp[[i]]
  write.csv(interesting, paste0("markers/Cluster_", i,"/Cluster_", i, "_conservative_markers.csv"))
}

# Find markers between each cluster and save a csv
# Genes up regulated in each cluster
markers.filtered_exp <- findMarkers(filtered_exp, test="wilcox", filtered_exp$cluster, direction="up", block = filtered_exp$Sample, pval.type = "some", row.data = rowData(filtered_exp))
for(i in names(markers.filtered_exp)){
  interesting <- markers.filtered_exp[[i]]
  write.csv(interesting, paste0("markers/Cluster_", i,"/Cluster_", i, "_permissive_markers_up.csv"))
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$FDR < 0.001, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = interesting[, "Ensembl"])
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/Cluster_", i,"/GOBP_markers_up_Cluster_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting[interesting$FDR < 0.001, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = interesting[, "Ensembl"])
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/Cluster_", i,"/GOMF_markers_up_Cluster_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting[interesting$FDR < 0.001, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = interesting[, "Ensembl"])
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/Cluster_", i,"/GOCC_markers_up_Cluster_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$FDR < 0.001, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = interesting[, "EntrezID"])
  write.csv(KEGGenrichsig, paste0("markers/Cluster_", i,"/KEGG_markers_up_Cluster_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$FDR < 0.001, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = interesting[, "EntrezID"])
  write.csv(REACTOMEenrichsig, paste0("markers/Cluster_", i,"/Reactome_markers_up_Cluster_", i, ".csv"))
  }

markers.filtered_exp <- findMarkers(filtered_exp, test="wilcox", filtered_exp$cluster, direction="down", block = filtered_exp$Sample, pval.type = "some", row.data = rowData(filtered_exp))

# Genes down regulated in each cluster
for(i in names(markers.filtered_exp)){
  interesting <- markers.filtered_exp[[i]]
  write.csv(interesting, paste0("markers/Cluster_", i,"/Cluster_", i, "_permissive_markers_down.csv"))

  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$FDR < 0.001, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = interesting[, "Ensembl"])
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/Cluster_", i,"/GOBP_markers_down_Cluster_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting[interesting$FDR < 0.001, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = interesting[, "Ensembl"])
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/Cluster_", i,"/GOMF_markers_down_Cluster_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting[interesting$FDR < 0.001, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = interesting[, "Ensembl"])
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/Cluster_", i,"/GOCC_markers_down_Cluster_", i, ".csv"))

  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$FDR < 0.001, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = interesting[, "EntrezID"])
  write.csv(KEGGenrichsig, paste0("markers/Cluster_", i,"/KEGG_markers_down_Cluster_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$FDR < 0.001, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = interesting[, "EntrezID"])
  write.csv(REACTOMEenrichsig, paste0("markers/Cluster_", i,"/Reactome_markers_down_Cluster_", i, ".csv"))
}

