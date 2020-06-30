#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Barplots of HALLMARK pathways.
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
library('Seurat')
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('viridis')
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
library('Rmagic')
library('msigdbr')
library('colorblindr')
library('ggridges')

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

exprs_no_batch <- removeBatchEffect(allDGEList$E, sampleRep)

# Scale gene expression values
zscores <- t(apply(exprs_no_batch, 1, scale))
colnames(zscores) <- colnames(allDGEList)
zscores <- zscores[! is.na(zscores[ ,1]), ]

# Import of significant HALLMARK pathways
# --------------------------------------------------------------------------

# Load Camera hallmark enrichment results for each tissue
if(place == "local") {
  liver_hallmark_up <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Liver/HALLMARK_markers_all_HC_cluster1.csv", header = TRUE)
  liver_hallmark_down <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Liver/HALLMARK_markers_all_HC_cluster2.csv", header = TRUE)
  } else {
  liver_hallmark_up <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Liver_Primary_Hallmark_camera.csv", header = TRUE)
}

liver_hallmark_up <- liver_hallmark_up[liver_hallmark_up$p.adjust < 0.05, ]
liver_hallmark_up$log10padj <- -log10(liver_hallmark_up$p.adjust)
liver_hallmark_up$Direction <- "Up"
liver_hallmark_down <- liver_hallmark_down[liver_hallmark_down$p.adjust < 0.05, ]
liver_hallmark_down$log10padj <- log10(liver_hallmark_down$p.adjust)
liver_hallmark_down$Direction <- "Down"

liver_hallmark <- rbind(liver_hallmark_up, liver_hallmark_down)
liver_hallmark <- liver_hallmark[order(liver_hallmark$log10padj, decreasing = FALSE), ]
liver_hallmark$Description <- factor(liver_hallmark$Description, levels = liver_hallmark$Description)


ggplot(liver_hallmark, aes(x = Description, y = log10padj)) + 
  geom_bar(aes(fill = Direction), stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(c(-15, 10)) +
  coord_flip() +
  theme_classic() +
  ggsave("markers/HC_Total_Liver/Liver_hallmark_pathways.pdf")


# Heatmaps of significant HALLMARK pathways
# --------------------------------------------------------------------------

# mSigDB lists
h_df <- msigdbr(species = "Homo sapiens", category = "H")

# Subset to genes active in any tissue
h_df <- h_df[is.element(h_df$gene_symbol, rowData_summed$GeneSymbol), ]

# Add rownames to hallmark lists
idx <- match(h_df$gene_symbol, rowData_summed$GeneSymbol)
h_df$Names <- rownames(rowData_summed) [idx]

# Split hallmark dataframe into lists of Names by geneset
h_list <- h_df %>% split(x = .$Names, f = .$gs_name)

# Make a list of only differentially regulated hallmark datasets 
h_sig_list <- h_list[names(h_list) %in% liver_hallmark$Description]

# Make heatmaps for each DEG list
for(pw in names(h_sig_list)) {
  zscores_GOI <- zscores[h_sig_list[[pw]], ]
  pdf(paste0("HALLMARK_", pw,"_expression_heatmap.pdf"))
  heatmap.2(zscores_GOI, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="row", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-2.5, 2.5, 1), col = hcl.colors(n = 5, palette = "Blue-Red 2"))
  dev.off()
}

# Import of significant HALLMARK pathways
# --------------------------------------------------------------------------

# Load Camera hallmark enrichment results for each tissue
if(place == "local") {
  lung_hallmark_up <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Lung/HALLMARK_markers_all_HC_cluster1.csv", header = TRUE)
  lung_hallmark_down <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Lung/HALLMARK_markers_all_HC_cluster2.csv", header = TRUE)
} else {
  lung_hallmark_up <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Lung_Primary_Hallmark_camera.csv", header = TRUE)
}

lung_hallmark_up <- lung_hallmark_up[lung_hallmark_up$p.adjust < 0.05, ]
lung_hallmark_up$log10padj <- -log10(lung_hallmark_up$p.adjust)
lung_hallmark_up$Direction <- "Up"
lung_hallmark_down <- lung_hallmark_down[lung_hallmark_down$p.adjust < 0.05, ]
lung_hallmark_down$log10padj <- log10(lung_hallmark_down$p.adjust)
lung_hallmark_down$Direction <- "Down"

lung_hallmark <- rbind(lung_hallmark_up, lung_hallmark_down)
lung_hallmark <- lung_hallmark[order(lung_hallmark$log10padj, decreasing = FALSE), ]
lung_hallmark$Description <- factor(lung_hallmark$Description, levels = lung_hallmark$Description)


ggplot(lung_hallmark, aes(x = Description, y = log10padj)) + 
  geom_bar(aes(fill = Direction), stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(c(-15, 10)) +
  coord_flip() +
  theme_classic() +
  ggsave("markers/HC_Total_Lung/Lung_hallmark_pathways.pdf")


# Heatmaps of significant HALLMARK pathways
# --------------------------------------------------------------------------

# mSigDB lists
h_df <- msigdbr(species = "Homo sapiens", category = "H")

# Subset to genes active in any tissue
h_df <- h_df[is.element(h_df$gene_symbol, rowData_summed$GeneSymbol), ]

# Add rownames to hallmark lists
idx <- match(h_df$gene_symbol, rowData_summed$GeneSymbol)
h_df$Names <- rownames(rowData_summed) [idx]

# Split hallmark dataframe into lists of Names by geneset
h_list <- h_df %>% split(x = .$Names, f = .$gs_name)

# Make a list of only differentially regulated hallmark datasets 
h_sig_list <- h_list[names(h_list) %in% lung_hallmark$Description]

# Make heatmaps for each DEG list
for(pw in names(h_sig_list)) {
  zscores_GOI <- zscores[h_sig_list[[pw]], ]
  pdf(paste0("HALLMARK_", pw,"_expression_heatmap.pdf"))
  heatmap.2(zscores_GOI, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="row", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-2.5, 2.5, 1), col = hcl.colors(n = 5, palette = "Blue-Red 2"))
  dev.off()
}

# Import of significant HALLMARK pathways
# --------------------------------------------------------------------------

# Load Camera hallmark enrichment results for each tissue
if(place == "local") {
  primary_hallmark_up <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Primary/HALLMARK_markers_all_HC_cluster2.csv", header = TRUE)
  primary_hallmark_down <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Primary/HALLMARK_markers_all_HC_cluster1.csv", header = TRUE)
} else {
  primary_hallmark_up <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Primary_Primary_Hallmark_camera.csv", header = TRUE)
}

primary_hallmark_up <- primary_hallmark_up[primary_hallmark_up$p.adjust < 0.05, ]
primary_hallmark_up$log10padj <- -log10(primary_hallmark_up$p.adjust)
primary_hallmark_up$Direction <- "Up"
primary_hallmark_down <- primary_hallmark_down[primary_hallmark_down$p.adjust < 0.05, ]
primary_hallmark_down$log10padj <- log10(primary_hallmark_down$p.adjust)
primary_hallmark_down$Direction <- "Down"

primary_hallmark <- rbind(primary_hallmark_up, primary_hallmark_down)
primary_hallmark <- primary_hallmark[order(primary_hallmark$log10padj, decreasing = FALSE), ]
primary_hallmark$Description <- factor(primary_hallmark$Description, levels = primary_hallmark$Description)


ggplot(primary_hallmark, aes(x = Description, y = log10padj)) + 
  geom_bar(aes(fill = Direction), stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(c(-15, 10)) +
  coord_flip() +
  theme_classic() +
  ggsave("markers/HC_Total_Primary/Primary_hallmark_pathways.pdf")


# Heatmaps of significant HALLMARK pathways
# --------------------------------------------------------------------------

# mSigDB lists
h_df <- msigdbr(species = "Homo sapiens", category = "H")

# Subset to genes active in any tissue
h_df <- h_df[is.element(h_df$gene_symbol, rowData_summed$GeneSymbol), ]

# Add rownames to hallmark lists
idx <- match(h_df$gene_symbol, rowData_summed$GeneSymbol)
h_df$Names <- rownames(rowData_summed) [idx]

# Split hallmark dataframe into lists of Names by geneset
h_list <- h_df %>% split(x = .$Names, f = .$gs_name)

# Make a list of only differentially regulated hallmark datasets 
h_sig_list <- h_list[names(h_list) %in% primary_hallmark$Description]

# Make heatmaps for each DEG list
for(pw in names(h_sig_list)) {
  zscores_GOI <- zscores[h_sig_list[[pw]], ]
  pdf(paste0("HALLMARK_", pw,"_expression_heatmap.pdf"))
  heatmap.2(zscores_GOI, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="row", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-2.5, 2.5, 1), col = hcl.colors(n = 5, palette = "Blue-Red 2"))
  dev.off()
}


