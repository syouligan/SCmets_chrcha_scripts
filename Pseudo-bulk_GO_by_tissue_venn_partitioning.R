#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Look at GO terms enriched among DEGs between tissues (pseudobulk) with venn partitioning
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
library('gplots')

# Make list of differentially expressed genes
liver_primary <- read.csv("pseudo-bulk_DGE/Liver_Primary_DEG_0LFC.csv", header = TRUE)
liver_primary$X <- as.character(liver_primary$X)
liver_primary_DGE <- liver_primary[liver_primary$adj.P.Val < 0.05 & liver_primary$logFC > 0, ]
ln_primary <- read.csv("pseudo-bulk_DGE/LN_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary_DGE <- ln_primary[ln_primary$adj.P.Val < 0.05 & ln_primary$logFC > 0, ]
lung_primary <- read.csv("pseudo-bulk_DGE/Lung_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary_DGE <- lung_primary[lung_primary$adj.P.Val < 0.05 & lung_primary$logFC > 0, ]

DGE_list <- list("Liver" = liver_primary_DGE$X, "LN" = ln_primary_DGE$X, "Lung" = lung_primary_DGE$X)

pdf("pseudo-bulk_DGE/DGE_psuedo-bulk_venn.pdf")
ItemsList <- venn(DGE_list, show.plot = TRUE)
dev.off()

lists <- attributes(ItemsList)$intersections

DGE_list <- list("Liver&Lung" = lists[["Liver:Lung"]], "LiverOnly" = lists[["Liver"]], "LungOnly" = lists[["Lung"]], "All" = lists[["Liver:LN:Lung"]])

# Perform enrichment analysis on each element of list
for(i in names(DGE_list)){
  interesting <- DGE_list[[i]]

  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = liver_primary[, "X"])
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("pseudo-bulk_DGE/", i,"_GOBP_pseudo-bulk.csv"))
  
  MFenrich <- enrichGO(interesting,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = liver_primary[, "X"])
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("pseudo-bulk_DGE/", i,"_GOMF_pseudo-bulk.csv"))
  
  CCenrich <- enrichGO(interesting,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = liver_primary[, "X"])
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("pseudo-bulk_DGE/", i,"_GOCC_pseudo-bulk.csv"))
  
  # # Perform KEGG enrichment analysis
  # KEGGenrichsig <- enrichKEGG(interesting[interesting$adj.P.Val < 0.05, "EntrezID"],
  #                             organism = "hsa",
  #                             keyType = "kegg",
  #                             pvalueCutoff = 0.05,
  #                             pAdjustMethod = "bonferroni",
  #                             universe = interesting[, "EntrezID"])
  # write.csv(KEGGenrichsig, paste0("pseudo-bulk_DGE/", i,"_KEGG_pseudo-bulk.csv"))
  # 
  # # Perform REACTOME enrichment analysis
  # REACTOMEenrichsig <- enrichPathway(interesting[interesting$adj.P.Val < 0.05, "EntrezID"],
  #                                    organism = "human",
  #                                    pvalueCutoff = 0.05,
  #                                    pAdjustMethod = "bonferroni",
  #                                    # readable = TRUE,
  #                                    universe = interesting[, "EntrezID"])
  # write.csv(REACTOMEenrichsig, paste0("pseudo-bulk_DGE/", i,"_Reactome_pseudo-bulk.csv"))
}

