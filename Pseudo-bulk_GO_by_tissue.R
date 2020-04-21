#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Look at GO terms enriched among DEGs between tissues (pseudobulk)
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data")
  place <- "wolfpack"
}

# Libraries
library('dplyr')
library('tidyr')
library('ggplot2')
library('readr')
library('Matrix')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')

# Make list of differentially expressed genes
liver_primary <- read.csv("pseudo-bulk_DGE/Liver_Primary_DEG_0LFC.csv", header = TRUE)
liver_primary$EntrezID <- mapIds(org.Hs.eg.db, keys=as.character(liver_primary$X), column="ENTREZID", keytype="SYMBOL", multiVals="first")

ln_primary <- read.csv("pseudo-bulk_DGE/LN_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary$EntrezID <- mapIds(org.Hs.eg.db, keys=as.character(ln_primary$X), column="ENTREZID", keytype="SYMBOL", multiVals="first")

lung_primary <- read.csv("pseudo-bulk_DGE/Lung_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary$EntrezID <- mapIds(org.Hs.eg.db, keys=as.character(lung_primary$X), column="ENTREZID", keytype="SYMBOL", multiVals="first")

liver_lung <- read.csv("pseudo-bulk_DGE/Liver_Lung_DEG_0LFC.csv", header = TRUE)
liver_lung$EntrezID <- mapIds(org.Hs.eg.db, keys=as.character(liver_lung$X), column="ENTREZID", keytype="SYMBOL", multiVals="first")

liver_ln <- read.csv("pseudo-bulk_DGE/Liver_LN_DEG_0LFC.csv", header = TRUE)
liver_ln$EntrezID <- mapIds(org.Hs.eg.db, keys=as.character(liver_ln$X), column="ENTREZID", keytype="SYMBOL", multiVals="first")

ln_lung <- read.csv("pseudo-bulk_DGE/LN_Lung_DEG_0LFC.csv", header = TRUE)
ln_lung$EntrezID <- mapIds(org.Hs.eg.db, keys=as.character(ln_lung$X), column="ENTREZID", keytype="SYMBOL", multiVals="first")

DGE_list <- list("Liver" = liver_primary, "LN" = ln_primary, "Lung" = lung_primary, "Liver_Lung" = liver_lung, "Liver_LN" = liver_ln, "LN_Lung" = ln_lung)

# Perform enrichment analysis on each element of list
for(i in names(DGE_list)){
  interesting <- DGE_list[[i]]
  interesting$X <- as.character(interesting$X)
  interesting$EntrezID <- as.character(interesting$EntrezID)
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "X"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "X"])
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("pseudo-bulk_DGE/", i,"_GOBP_pseudo-bulk_up.csv"))
  
  MFenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "X"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "X"])
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("pseudo-bulk_DGE/", i,"_GOMF_pseudo-bulk_up.csv"))
  
  CCenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "X"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "X"])
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("pseudo-bulk_DGE/", i,"_GOCC_pseudo-bulk_up.csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = interesting[, "EntrezID"])
  write.csv(KEGGenrichsig, paste0("pseudo-bulk_DGE/", i,"_KEGG_pseudo-bulk_up.csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     # readable = TRUE,
                                     universe = interesting[, "EntrezID"])
  write.csv(REACTOMEenrichsig, paste0("pseudo-bulk_DGE/", i,"_Reactome_pseudo-bulk_up.csv"))
}

  
# Perform enrichment analysis on each element of list
for(i in names(DGE_list)){
  interesting <- DGE_list[[i]]
  interesting$X <- as.character(interesting$X)
  interesting$EntrezID <- as.character(interesting$EntrezID)
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "X"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "X"])
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("pseudo-bulk_DGE/", i,"_GOBP_pseudo-bulk_down.csv"))
  
  MFenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "X"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "X"])
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("pseudo-bulk_DGE/", i,"_GOMF_pseudo-bulk_down.csv"))
  
  CCenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "X"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "X"])
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("pseudo-bulk_DGE/", i,"_GOCC_pseudo-bulk_down.csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = interesting[, "EntrezID"])
  write.csv(KEGGenrichsig, paste0("pseudo-bulk_DGE/", i,"_KEGG_pseudo-bulk_down.csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     # readable = TRUE,
                                     universe = interesting[, "EntrezID"])
  write.csv(REACTOMEenrichsig, paste0("pseudo-bulk_DGE/", i,"_Reactome_pseudo-bulk_down.csv"))
}


  