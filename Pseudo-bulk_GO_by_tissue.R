#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Look at GO terms enriched among DEGs between tissues (pseudobulk)
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
library('ggplot2')
library('readr')
library('Matrix')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')
library('gplots')
library('msigdbr')

# Make list of differentially expressed genes
liver_primary <- read.csv("pseudo-bulk_DGE/Liver_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary <- read.csv("pseudo-bulk_DGE/LN_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary <- read.csv("pseudo-bulk_DGE/Lung_Primary_DEG_0LFC.csv", header = TRUE)
liver_lung <- read.csv("pseudo-bulk_DGE/Liver_Lung_DEG_0LFC.csv", header = TRUE)
liver_ln <- read.csv("pseudo-bulk_DGE/Liver_LN_DEG_0LFC.csv", header = TRUE)
ln_lung <- read.csv("pseudo-bulk_DGE/LN_Lung_DEG_0LFC.csv", header = TRUE)

DGE_list <- list("Liver_Primary" = liver_primary, "LN_Primary" = ln_primary, "Lung_Primary" = lung_primary, "Liver_Lung" = liver_lung, "Liver_LN" = liver_ln, "LN_Lung" = ln_lung)

dir.create("pseudo-bulk_DGE/GO")

# Makes hallmark geneset for enrichment testing
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_t2g <- h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()

# Perform enrichment analysis on each element of list
for(i in names(DGE_list)){
  interesting <- DGE_list[[i]]
  interesting$GeneSymbol <- as.character(interesting$GeneSymbol)
  interesting$EntrezID <- as.character(interesting$EntrezID)
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "GeneSymbol"])
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_GOBP_pseudo-bulk_up.csv"))
  
  MFenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "GeneSymbol"])
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_GOMF_pseudo-bulk_up.csv"))
  
  CCenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "GeneSymbol"])
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_GOCC_pseudo-bulk_up.csv"))
  
  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = interesting[, "GeneSymbol"])
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_HALLMARK_pseudo-bulk_up.csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = interesting[, "EntrezID"])
  write.csv(KEGGenrichsig, paste0("pseudo-bulk_DGE/GO/", i,"_KEGG_pseudo-bulk_up.csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC > 0, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     # readable = TRUE,
                                     universe = interesting[, "EntrezID"])
  write.csv(REACTOMEenrichsig, paste0("pseudo-bulk_DGE/GO/", i,"_Reactome_pseudo-bulk_up.csv"))
}

  
# Perform enrichment analysis on each element of list
for(i in names(DGE_list)){
  interesting <- DGE_list[[i]]
  interesting$GeneSymbol <- as.character(interesting$GeneSymbol)
  interesting$EntrezID <- as.character(interesting$EntrezID)
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "GeneSymbol"])
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_GOBP_pseudo-bulk_down.csv"))
  
  MFenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "GeneSymbol"])
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_GOMF_pseudo-bulk_down.csv"))
  
  CCenrich <- enrichGO(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = interesting[, "GeneSymbol"])
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_GOCC_pseudo-bulk_down.csv"))
  
  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = interesting[, "GeneSymbol"])
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_HALLMARK_pseudo-bulk_down.csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = interesting[, "EntrezID"])
  write.csv(KEGGenrichsig, paste0("pseudo-bulk_DGE/GO/", i,"_KEGG_pseudo-bulk_down.csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$adj.P.Val < 0.05 & interesting$logFC < 0, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     # readable = TRUE,
                                     universe = interesting[, "EntrezID"])
  write.csv(REACTOMEenrichsig, paste0("pseudo-bulk_DGE/GO/", i,"_Reactome_pseudo-bulk_down.csv"))
}
  