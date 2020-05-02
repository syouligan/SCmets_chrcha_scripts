#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Makes heatmaps of Hallmark pathways across tissues
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
library('Seurat')
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('org.Hs.eg.db')
library('msigdbr')
library('gtools')
library('cowplot')


# Load prefiltered and clustered Seurat Object
if(place == "local") {
  filtered_exp <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/Prefiltered_experiment_all_seurat_integrated.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Normalise counts for plotting and save object
filtered_exp <- SCTransform(filtered_exp, verbose = TRUE, vars.to.regress = c("Replicate", "Lib_size", "S.Score", "G2M.Score", "digest_stress1"), variable.features.n = 5000, return.only.var.genes = FALSE, new.assay.name = "SCT_plotting")

if(place == "local") {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_practice_seurat_integrated_plotting.rds")
} else if(place == "wolfpack") {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_all_seurat_integrated_plotting.rds")
} else {
  print("Not overwritten")
}

# Make heatmap of DEGs by Pseudobulk
# --------------------------------------------------------------------------

# Make list of top 100 differentially expressed genes
liver_primary <- read.csv("pseudo-bulk_DGE/Liver_Primary_DEG_0LFC.csv", header = TRUE)
liver_primary <- liver_primary[liver_primary$adj.P.Val < 0.05, ]
liver_primary <- liver_primary[order(liver_primary$t),]
liver_primary <- as.character(liver_primary[, "X"])

ln_primary <- read.csv("pseudo-bulk_DGE/LN_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary <- ln_primary[ln_primary$adj.P.Val < 0.05, ]
ln_primary <- ln_primary[order(ln_primary$t),]
ln_primary <- as.character(ln_primary[, "X"])

lung_primary <- read.csv("pseudo-bulk_DGE/Lung_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary <- lung_primary[lung_primary$adj.P.Val < 0.05, ]
lung_primary <- lung_primary[order(lung_primary$t),]
lung_primary <- as.character(lung_primary[, "X"])

liver_lung <- read.csv("pseudo-bulk_DGE/Liver_Lung_DEG_0LFC.csv", header = TRUE)
liver_lung <- liver_lung[liver_lung$adj.P.Val < 0.05, ]
liver_lung <- liver_lung[order(liver_lung$t),]
liver_lung <- as.character(liver_lung[, "X"])

liver_ln <- read.csv("pseudo-bulk_DGE/Liver_LN_DEG_0LFC.csv", header = TRUE)
liver_ln <- liver_ln[liver_ln$adj.P.Val < 0.05, ]
liver_ln <- liver_ln[order(liver_ln$t),]
liver_ln <- as.character(liver_ln[, "X"])

ln_lung <- read.csv("pseudo-bulk_DGE/LN_Lung_DEG_0LFC.csv", header = TRUE)
ln_lung <- ln_lung[ln_lung$adj.P.Val < 0.05, ]
ln_lung <- ln_lung[order(ln_lung$t),]
ln_lung <- as.character(ln_lung[, "X"])

DGE_list <- list("Liver_Primary" = liver_primary, "LN_Primary" = ln_primary, "Lung_Primary" = lung_primary, "Liver_Lung" = liver_lung, "Liver_LN" = liver_ln, "LN_Lung" = ln_lung)

# Plot heatmaps for top 100 DEGs
for(pw in names(DGE_list)) {
  hm <- DoHeatmap(object = filtered_exp, features = as.character(DGE_list[[pw]]), group.by = "Tissue", slot = "scale.data", assay = "SCT_plotting")
  ggsave(paste0("pseudo-bulk_DGE/", pw, "_DEGs_heatmap.png"), plot = hm, device = "png")
}