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
library('NMF')

set.seed(100)

# Load Pseudo-bulk counts and row data SingleCellExperiment
exprs <- readRDS("Pseudo-bulk_whole_experiment_all.rds") # uses whole dataset if wolfpack
rowData_summed <- read.csv("Pseudo-bulk_rowData_all.csv", header = TRUE, row.names = 1)

# Make heatmap of DEGs between tissues and primary tumors
# --------------------------------------------------------------------------

# Find genes detected in all replicates of a given tissue type
active_liver <- rowSums(exprs[,1:4] >= 1) == 4
active_LN <- rowSums(exprs[,5:8] >= 1) == 4
active_lung <- rowSums(exprs[,9:12] >= 1) == 4
active_primary <- rowSums(exprs[,13:16] >= 1) == 4
geneActivity <- data.frame("Liver" = active_liver, "LN" = active_LN, "Lung" = active_lung, "Primary" = active_primary)

# Make sample info table
sampleName <- colnames(exprs)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTable <- data.frame(sampleName, sampleType, sampleRep)

# Make DGEList object and normalise using TMM
allDGEList <- DGEList(counts = exprs, group = sampleType)
allDGEList <- calcNormFactors(allDGEList, method = "TMM")

# Make design matrix, fit linear models, list comparisons
type <- factor(sampleTable$sampleType)
rep <- factor(sampleTable$sampleRep, levels = c(1:4))
design <- model.matrix(~0 + type)
rownames(design) <- sampleName
allDGEList <- voom(allDGEList, design, plot = FALSE)

# Removes batch effect and calculates the inverse of the log2 transformation
exprs_no_batch <- 2^(removeBatchEffect(allDGEList$E, sampleRep))

preprocessed.nmf <- nmf(exprs_no_batch, 2:10, nrun = 50, seed = 123456, .opt="vp8")
rand <- randomize(exprs_no_batch)
random.nmf <- nmf(rand, 2:10, nrun = 50, seed = 123456, .opt="vp8")

saveRDS(preprocessed.nmf, "Pseudo-bulk_NMF_organised_object.rds")
saveRDS(preprocessed.nmf, "Pseudo-bulk_NMF_random_object.rds")

plot(preprocessed.nmf, random.nmf)
consensusmap(preprocessed.nmf)

k4.nmf <- nmf(exprs_no_batch, 4, seed = 123456, .opt="vp8")

summary(k4.nmf, class=as.factor(sampleType))
k4.nmf.features <- extractFeatures(k4.nmf)
top5 <- extractFeatures(k4.nmf, 10L)

names(k4.nmf.features) <- c("cluster1", "cluster2", "cluster3", "cluster4")

# Scale gene expression values
zscores <- t(apply(exprs_no_batch, 1, scale))
colnames(zscores) <- colnames(allDGEList)
zscores <- zscores[! is.na(zscores[ ,1]), ]

# Make heatmaps of each cluster
for(i in names(k4.nmf.features)) {
  zscores_GOI <- zscores[rownames(exprs_no_batch)[k4.nmf.features[[i]]], ]
  pdf(paste0("NMF_", i,"_expression_heatmap.pdf"))
  heatmap.2(zscores_GOI, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="none", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-2.5, 2.5, 1), col = magma(5, direction = -1))
  dev.off()
}

# Make gene universe(s)
universe <- as.character(unique(rowData_summed[rowData_summed$Any_Active, "Ensembl"]))
universe_entrez <- as.character(unique(rowData_summed[rowData_summed$Any_Active, "EntrezID"]))
universe_genesymbol <- as.character(unique(rowData_summed[rowData_summed$Any_Active, "GeneSymbol"]))

# Makes hallmark geneset for enrichment testing
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_t2g <- h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()

for(i in names(k4.nmf.features)){
  dir.create(paste0("markers/NMF_", i))
  interesting <- data.frame(rownames(exprs_no_batch)[k4.nmf.features[[i]]])
  idx <- match(interesting[,1], rownames(rowData_summed))
  interesting$Ensembl <- rowData_summed$Ensembl [idx]
  interesting$EntrezID <- rowData_summed$EntrezID [idx]
  interesting$GeneSymbol <- rowData_summed$GeneSymbol [idx]
  write.csv(interesting, paste0("markers/NMF_", i,"/NMF_", i, "_markers.csv"))
  
  # All-regulated markers in each tissue
  # --------------------------------------------------------------------------
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "fdr",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/NMF_", i,"/GOBP_markers_all_NMF_", i, ".csv"))
  
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
  write.csv(GOMF_GOI, paste0("markers/NMF_", i,"/GOMF_markers_all_NMF_", i, ".csv"))
  
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
  write.csv(GOCC_GOI, paste0("markers/NMF_", i,"/GOCC_markers_all_NMF_", i, ".csv"))
  
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
  write.csv(Hallmark_GOI, paste0("markers/NMF_", i,"/HALLMARK_markers_all_NMF_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/NMF_", i,"/KEGG_markers_all_NMF_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/NMF_", i,"/Reactome_markers_all_NMF_", i, ".csv"))
}




