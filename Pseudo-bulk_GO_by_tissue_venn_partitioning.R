#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Look at GO terms enriched among DEGs between tissues (pseudobulk) with venn partitioning
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
liver_primary$GeneSymbol <- as.character(liver_primary$GeneSymbol)
liver_primary_DGE <- liver_primary[liver_primary$adj.P.Val < 0.05 & liver_primary$logFC > 0, ]
ln_primary <- read.csv("pseudo-bulk_DGE/LN_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary_DGE <- ln_primary[ln_primary$adj.P.Val < 0.05 & ln_primary$logFC > 0, ]
lung_primary <- read.csv("pseudo-bulk_DGE/Lung_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary_DGE <- lung_primary[lung_primary$adj.P.Val < 0.05 & lung_primary$logFC > 0, ]
DGE_list <- list("Liver" = liver_primary_DGE$GeneSymbol, "LN" = ln_primary_DGE$GeneSymbol, "Lung" = lung_primary_DGE$GeneSymbol)

# Venn diagram overlap DEGs between metastatic sites and tissues
pdf("pseudo-bulk_DGE/DGE_psuedo-bulk_venn.pdf")
ItemsList <- venn(DGE_list, show.plot = TRUE)
dev.off()

# List of DEGs partitioned by their venn segment
lists <- attributes(ItemsList)$intersections
DGE_list <- list("Liver&Lung" = liver_primary[is.element(liver_primary$GeneSymbol, lists[["Liver:Lung"]]), ], "Liver&LN" = liver_primary[is.element(liver_primary$GeneSymbol, lists[["Liver:LN"]]), ], "LN&Lung" = ln_primary[is.element(ln_primary$GeneSymbol, lists[["LN:Lung"]]), ], "LiverOnly" = liver_primary[is.element(liver_primary$GeneSymbol, lists[["Liver"]]), ], "LungOnly" = lung_primary[is.element(lung_primary$GeneSymbol, lists[["Lung"]]), ], "LNOnly" = ln_primary[is.element(ln_primary$GeneSymbol, lists[["LN"]]), ], "All" = ln_primary[is.element(ln_primary$GeneSymbol, lists[["Liver:LN:Lung"]]), ])
  
# Perform enrichment analysis on each element of list
universe <- as.character(unique(c(lung_primary$GeneSymbol, ln_primary$GeneSymbol, liver_primary$GeneSymbol)))
universe_entrez <- as.character(unique(c(lung_primary$EntrezID, ln_primary$EntrezID, liver_primary$EntrezID)))

# Makes hallmark geneset for enrichment testing
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_t2g <- h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()

for(i in names(DGE_list)){
  interesting <- DGE_list[[i]]

  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_GOBP_upVSprimary_pseudo-bulk.csv"))
  
  MFenrich <- enrichGO(interesting[, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_GOMF_upVSprimary_pseudo-bulk.csv"))
  
  CCenrich <- enrichGO(interesting[, "GeneSymbol"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       # readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_GOCC_upVSprimary_pseudo-bulk.csv"))
  
  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = universe)
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("pseudo-bulk_DGE/GO/", i,"_HALLMARK_upVSprimary_pseudo-bulk.csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("pseudo-bulk_DGE/GO/", i,"_KEGG_upVSprimary_pseudo-bulk.csv"))

  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     # readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("pseudo-bulk_DGE/GO/", i,"_Reactome_upVSprimary_pseudo-bulk.csv"))
}

