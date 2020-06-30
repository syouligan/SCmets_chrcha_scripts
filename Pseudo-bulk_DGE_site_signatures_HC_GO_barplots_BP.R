#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Barplots of GOBP pathways.
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
active_liver <- rowSums(exprs[,1:5] >= 1) == 4
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
rep <- factor(sampleTable$sampleRep, levels = c(1:5))
design <- model.matrix(~0 + type)
rownames(design) <- sampleName
allDGEList <- voom(allDGEList, design, plot = FALSE)

exprs_no_batch <- removeBatchEffect(allDGEList$E, sampleRep)

# Scale gene expression values
zscores <- t(apply(exprs_no_batch, 1, scale))
colnames(zscores) <- colnames(allDGEList)
zscores <- zscores[! is.na(zscores[ ,1]), ]

# Import of significant GOBP pathways
# --------------------------------------------------------------------------

# Load Camera GOBP enrichment results for each tissue
if(place == "local") {
  liver_GOBP_up <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Liver/GOBP_markers_all_HC_cluster1.csv", header = TRUE)
  liver_GOBP_down <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Liver/GOBP_markers_all_HC_cluster2.csv", header = TRUE)
} else {
  liver_GOBP_up <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Liver_Primary_GOBP_camera.csv", header = TRUE)
}

liver_GOBP_up <- liver_GOBP_up[1:5, ]
liver_GOBP_up$log10padj <- -log10(liver_GOBP_up$p.adjust)
liver_GOBP_up$Direction <- "Up"
liver_GOBP_down <- liver_GOBP_down[1:5, ]
liver_GOBP_down$log10padj <- log10(liver_GOBP_down$p.adjust)
liver_GOBP_down$Direction <- "Down"

liver_GOBP <- rbind(liver_GOBP_up, liver_GOBP_down)
liver_GOBP <- liver_GOBP[order(liver_GOBP$log10padj, decreasing = FALSE), ]
liver_GOBP$Description <- factor(liver_GOBP$Description, levels = liver_GOBP$Description)


ggplot(liver_GOBP, aes(x = Description, y = log10padj)) + 
  geom_bar(aes(fill = Direction), stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(c(-8, 9)) +
  coord_flip() +
  theme_classic() +
  ggsave("markers/HC_Total_Liver/Liver_GOBP_pathways.pdf")

# Import of significant GOBP pathways
# --------------------------------------------------------------------------

# Load Camera GOBP enrichment results for each tissue
if(place == "local") {
  lung_GOBP_up <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Lung/GOBP_markers_all_HC_cluster1.csv", header = TRUE)
  lung_GOBP_down <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Lung/GOBP_markers_all_HC_cluster2.csv", header = TRUE)
} else {
  lung_GOBP_up <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Lung_Primary_GOBP_camera.csv", header = TRUE)
}

lung_GOBP_up <- lung_GOBP_up[1:5, ]
lung_GOBP_up$log10padj <- -log10(lung_GOBP_up$p.adjust)
lung_GOBP_up$Direction <- "Up"
lung_GOBP_down <- lung_GOBP_down[1:5, ]
lung_GOBP_down$log10padj <- log10(lung_GOBP_down$p.adjust)
lung_GOBP_down$Direction <- "Down"

lung_GOBP <- rbind(lung_GOBP_up, lung_GOBP_down)
lung_GOBP <- lung_GOBP[order(lung_GOBP$log10padj, decreasing = FALSE), ]
lung_GOBP$Description <- factor(lung_GOBP$Description, levels = lung_GOBP$Description)


ggplot(lung_GOBP, aes(x = Description, y = log10padj)) + 
  geom_bar(aes(fill = Direction), stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(c(-8, 9)) +
  coord_flip() +
  theme_classic() +
  ggsave("markers/HC_Total_Lung/Lung_GOBP_pathways.pdf")

# Import of significant GOBP pathways
# --------------------------------------------------------------------------

# Load Camera GOBP enrichment results for each tissue
if(place == "local") {
  primary_GOBP_up <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Primary/GOBP_markers_all_HC_cluster2.csv", header = TRUE)
  primary_GOBP_down <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/markers/HC_Total_Primary/GOBP_markers_all_HC_cluster1.csv", header = TRUE)
} else {
  primary_GOBP_up <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Primary_Primary_GOBP_camera.csv", header = TRUE)
}

primary_GOBP_up <- primary_GOBP_up[1:5, ]
primary_GOBP_up$log10padj <- -log10(primary_GOBP_up$p.adjust)
primary_GOBP_up$Direction <- "Up"
primary_GOBP_down <- primary_GOBP_down[1:5, ]
primary_GOBP_down$log10padj <- log10(primary_GOBP_down$p.adjust)
primary_GOBP_down$Direction <- "Down"

primary_GOBP <- rbind(primary_GOBP_up, primary_GOBP_down)
primary_GOBP <- primary_GOBP[order(primary_GOBP$log10padj, decreasing = FALSE), ]
primary_GOBP$Description <- factor(primary_GOBP$Description, levels = primary_GOBP$Description)


ggplot(primary_GOBP, aes(x = Description, y = log10padj)) + 
  geom_bar(aes(fill = Direction), stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ylim(c(-8, 9)) +
  coord_flip() +
  theme_classic() +
  ggsave("markers/HC_Total_Primary/Primary_GOBP_pathways.pdf")


