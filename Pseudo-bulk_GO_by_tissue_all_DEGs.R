#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Plot heatmaps of DEGs between tissues
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/")
  place <- "wolfpack"
}

# Libraries
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('Matrix')
library('limma')
library('edgeR')
library('gtools')
library('gplots')
library('msigdbr')
library('org.Hs.eg.db')
library('fgsea')
library('data.table')
library('sva')
library('viridis')
library('dendextend')
library('clusterProfiler')
library('ReactomePA')
library('Matrix')
library('phateR')

set.seed(100)

# Load Pseudo-bulk counts and row data SingleCellExperiment
exprs <- readRDS("Pseudo-bulk_whole_experiment_all.rds") # uses whole dataset if wolfpack
rowData_summed <- read.csv("Pseudo-bulk_rowData_all.csv", header = TRUE, row.names = 1)

# Load DGE data
# --------------------------------------------------------------------------
# Make list of top 100 differentially expressed genes
liver_primary <- read.csv("Liver_Primary_DEG_0LFC.csv", header = TRUE)
liver_primary <- liver_primary[liver_primary$adj.P.Val < 0.05, ]
liver_primary <- liver_primary[order(liver_primary$t),]
liver_primary <- as.character(liver_primary[, "X"])

ln_primary <- read.csv("LN_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary <- ln_primary[ln_primary$adj.P.Val < 0.05, ]
ln_primary <- ln_primary[order(ln_primary$t),]
ln_primary <- as.character(ln_primary[, "X"])

lung_primary <- read.csv("Lung_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary <- lung_primary[lung_primary$adj.P.Val < 0.05, ]
lung_primary <- lung_primary[order(lung_primary$t),]
lung_primary <- as.character(lung_primary[, "X"])

liver_lung <- read.csv("Liver_Lung_DEG_0LFC.csv", header = TRUE)
liver_lung <- liver_lung[liver_lung$adj.P.Val < 0.05, ]
liver_lung <- liver_lung[order(liver_lung$t),]
liver_lung <- as.character(liver_lung[, "X"])

liver_ln <- read.csv("Liver_LN_DEG_0LFC.csv", header = TRUE)
liver_ln <- liver_ln[liver_ln$adj.P.Val < 0.05, ]
liver_ln <- liver_ln[order(liver_ln$t),]
liver_ln <- as.character(liver_ln[, "X"])

ln_lung <- read.csv("LN_Lung_DEG_0LFC.csv", header = TRUE)
ln_lung <- ln_lung[ln_lung$adj.P.Val < 0.05, ]
ln_lung <- ln_lung[order(ln_lung$t),]
ln_lung <- as.character(ln_lung[, "X"])

total_DGE <- as.character(unique(c(liver_primary, ln_primary, lung_primary, liver_lung, ln_lung, liver_ln)))

DGE_list <- list("Total" = total_DGE)

# Make gene universe(s)
universe <- as.character(unique(rowData_summed[rowData_summed$Any_Active, "Ensembl"]))
universe_entrez <- as.character(unique(rowData_summed[rowData_summed$Any_Active, "EntrezID"]))
universe_genesymbol <- as.character(unique(rowData_summed[rowData_summed$Any_Active, "GeneSymbol"]))

# Makes hallmark geneset for enrichment testing
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_t2g <- h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()

# Makes hallmark geneset for enrichment testing
metabolic_pathways <- read.csv("~/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Metabolic_pathways_genes.csv", header = TRUE)
metabolic_pathways$Gene <- as.character(metabolic_pathways$Gene)
metabolic_pathways$Metabolic_pathway <- as.character(metabolic_pathways$Metabolic_pathway)

# Perform enrichment analysis on markers
# --------------------------------------------------------------------------
dir.create("markers")

for(i in names(DGE_list)){
  dir.create(paste0("markers/HC_", i))
  interesting <- data.frame(DGE_list[[i]])
  idx <- match(interesting[,1], rownames(rowData_summed))
  interesting$Ensembl <- rowData_summed$Ensembl [idx]
  interesting$EntrezID <- rowData_summed$EntrezID [idx]
  interesting$GeneSymbol <- rowData_summed$GeneSymbol [idx]
 
  # All-regulated markers in each tissue
  # --------------------------------------------------------------------------
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[, "Ensembl"],
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
  write.csv(GOBP_GOI, paste0("markers/HC_", i,"/GOBP_markers_all_HC_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting[, "Ensembl"],
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
  write.csv(GOMF_GOI, paste0("markers/HC_", i,"/GOMF_markers_all_HC_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting[, "Ensembl"],
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
  write.csv(GOCC_GOI, paste0("markers/HC_", i,"/GOCC_markers_all_HC_", i, ".csv"))
  
  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = universe_genesymbol)
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("markers/HC_", i,"/HALLMARK_markers_all_HC_", i, ".csv"))

  # Perform Metabolic_pathways enrichment analysis
  Metabolic_pathwaysEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                                       TERM2GENE = metabolic_pathways,
                                       pvalueCutoff = 1,
                                       qvalueCutoff = 1,
                                       pAdjustMethod = "bonferroni",
                                       minGSSize = 10,
                                       maxGSSize = 500,
                                       universe = universe_genesymbol)
  Metabolic_pathways_GOI <- as.data.frame(Metabolic_pathwaysEnrich@result)
  write.csv(Metabolic_pathways_GOI, paste0("markers/HC_", i,"/Metabolic_pathways_all_HC_", i, ".csv"))
    
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/HC_", i,"/KEGG_markers_all_HC_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/HC_", i,"/Reactome_markers_all_HC_", i, ".csv"))
}
