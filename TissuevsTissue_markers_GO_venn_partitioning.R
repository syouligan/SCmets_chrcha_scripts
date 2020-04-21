#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Look at GO terms enriched among DEGs between tissues (pseudobulk) with venn partitioning
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data")
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
liver_primary <- read.csv("markers/Tissue_Liver_Primary/Tissue_Liver_Primary_downregulated_markers.csv", header = TRUE)
liver_primary_DGE <- liver_primary[liver_primary$p_val_adj < 0.05, ]
ln_primary <- read.csv("markers/Tissue_LN_Primary/Tissue_LN_Primary_downregulated_markers.csv", header = TRUE)
ln_primary_DGE <- ln_primary[ln_primary$p_val_adj < 0.05, ]
lung_primary <- read.csv("markers/Tissue_Lung_Primary/Tissue_Lung_Primary_downregulated_markers.csv", header = TRUE)
lung_primary_DGE <- lung_primary[lung_primary$p_val_adj < 0.05, ]

DGE_list <- list("Liver" = liver_primary_DGE$X, "LN" = ln_primary_DGE$X, "Lung" = lung_primary_DGE$X)

pdf("markers/Up_markers_venn.pdf")
ItemsList <- venn(DGE_list, show.plot = TRUE)
dev.off()

lists <- attributes(ItemsList)$intersections

DGE_list <- list("Liver&Lung" = lists[["Liver:Lung"]], "Liver&LN" = lists[["Liver:LN"]], "LN&Lung" = lists[["LN:Lung"]], "LiverOnly" = lists[["Liver"]], "LungOnly" = lists[["Lung"]], "LNOnly" = lists[["LN"]], "All" = lists[["Liver:LN:Lung"]])

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
                       universe = as.character(liver_primary[, "X"]))
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/", i,"_GOBP_up_markers.csv"))
  
  MFenrich <- enrichGO(interesting,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = as.character(liver_primary[, "X"]))
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/", i,"_GOMF_up_markers.csv"))
  
  CCenrich <- enrichGO(interesting,
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = as.character(liver_primary[, "X"]))
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/", i,"_GOCC_up_markers.csv"))
  
  # Perform KEGG enrichment analysis
  interesting_entrez <- liver_primary[is.element(liver_primary$X, interesting), "EntrezID"]
  interesting_entrez <- interesting_entrez[!is.na(interesting_entrez)]
  KEGGenrichsig <- enrichKEGG(interesting_entrez,
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = as.character(liver_primary[! is.na(liver_primary$EntrezID), "EntrezID"]))
  write.csv(KEGGenrichsig, paste0("markers/", i,"_KEGG_up_markers.csv"))

  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting_entrez,
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     # readable = TRUE,
                                     universe = as.character(liver_primary[! is.na(liver_primary$EntrezID), "EntrezID"]))
  write.csv(REACTOMEenrichsig, paste0("markers/", i,"_Reactome_up_markers.csv"))
}

