#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells using Seurat.
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data")
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

# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp.list <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_QC_experiment_practice_filtered_exp_list.rds") # uses practice data if local
  all_features <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/All_features.rds") # uses practice data if local
} else {
  filtered_exp.list <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/Prefiltered_QC_experiment_all_filtered_exp_list.rds") # uses whole dataset if wolfpack
  all_features <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/All_features.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 1000000*1024^2)
  
}

# Normalise transform counts within each experiment. Note: will not overwrite if already exists.
# --------------------------------------------------------------------------

all_features <- as.character(all_features[,,drop = TRUE])

# Integrate datasets based on highly correlated features
filtered_exp.features <- SelectIntegrationFeatures(object.list = filtered_exp.list, nfeatures = 5000)
filtered_exp.list <- PrepSCTIntegration(object.list = filtered_exp.list, verbose = TRUE, anchor.features = filtered_exp.features)
filtered_exp.list <- lapply(X = filtered_exp.list, FUN = RunPCA, verbose = TRUE, features = filtered_exp.features) # Perform PCA on each object individually (needed for rpca)
reference_datasets <- which(names(filtered_exp.list) == "3")
filtered_exp.anchors <- FindIntegrationAnchors(object.list = filtered_exp.list, normalization.method = "SCT", anchor.features = filtered_exp.features, verbose = TRUE, reduction = "rpca", reference = reference_datasets)
filtered_exp.integrated <- IntegrateData(anchorset = filtered_exp.anchors, normalization.method = "SCT", verbose = TRUE)

if(place == "local" & exists("phate.out")) {
  write.csv(t(filtered_exp.integrated@assays$integrated@scale.data), "Prefiltered_experiment_practice_seurat_integrated_SCT.csv")
  write.csv(filtered_exp.integrated@meta.data, "Prefiltered_experiment_practice_seurat_integrated_colData.csv")
} else {
  write.csv(t(filtered_exp.integrated@assays$integrated@scale.data), "Prefiltered_experiment_all_seurat_integrated_SCT.csv")
  write.csv(filtered_exp.integrated@meta.data, "Prefiltered_experiment_all_seurat_integrated_colData.csv")
}

# Save seurat objects
# --------------------------------------------------------------------------

if(place == "local") {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_practice_seurat_integrated_for_MELD.rds")
} else if(place == "wolfpack") {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_all_seurat_integrated_for_MELD.rds")
} else {
  print("Not overwritten")
}

# Make and save log normalised integrated object for imputation
# --------------------------------------------------------------------------
rm(filtered_exp.integrated)
rm()
filtered_exp.integrated.ln <- IntegrateData(anchorset = filtered_exp.anchors, normalization.method = "LogNormalize", verbose = TRUE)
if(place == "local" & exists("phate.out")) {
  write.csv(t(filtered_exp.integrated.ln@assays$integrated@data), "Prefiltered_experiment_practice_seurat_integrated_LogNorm.csv")
} else {
  write.csv(t(filtered_exp.integrated.ln@assays$integrated@data), "Prefiltered_experiment_all_seurat_integrated_LogNorm.csv")
}



