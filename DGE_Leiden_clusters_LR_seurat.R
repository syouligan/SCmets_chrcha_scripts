#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Perform differential expression (likelihood ratio) between clusters identified in MELD
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/meld/all_data")
  place <- "wolfpack"
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Libraries
library('Seurat')
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('batchelor')
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
library('msigdbr')
library('fgsea')
library('clusterProfiler')

if (place == 'local'){
  obs <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/practice_all_data/Prefiltered_experiment_practice_normalised_AnnData_PHATE_clusterings_obs.csv", header = TRUE, row.names = 1)
  var <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/practice_all_data/Prefiltered_experiment_practice_normalised_AnnData_PHATE_clusterings_var.csv", header = TRUE, row.names = 1)
  obsm <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/practice_all_data/Prefiltered_experiment_practice_normalised_AnnData_PHATE_clusterings_embeddings.csv", header = TRUE, row.names = 1)
  counts <- as(readMM("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/practice_all_data/Prefiltered_experiment_practice_normalised_AnnData_PHATE_clusterings_normalised_counts.mtx"), "dgCMatrix")
} else {
  obs <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/meld/all_data/Prefiltered_experiment_all_normalised_AnnData_PHATE_clusterings_obs.csv", header = TRUE, row.names = 1)
  var <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/meld/all_data/Prefiltered_experiment_all_normalised_AnnData_PHATE_clusterings_var.csv", header = TRUE, row.names = 1)
  obsm <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/meld/all_data/Prefiltered_experiment_all_normalised_AnnData_PHATE_clusterings_embeddings.csv", header = TRUE, row.names = 1)
  counts <- as(readMM("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/meld/all_data/Prefiltered_experiment_all_normalised_AnnData_PHATE_clusterings_normalised_counts.mtx"), "dgCMatrix")
}

# Perform differential gene expression
filtered_exp_sce <- SingleCellExperiment(assays = list('logcounts' = t(counts)), colData = obs, rowData = var, reducedDims = SimpleList('PHATE' = as.matrix(obsm)))
saveRDS(filtered_exp_sce, 'Prefiltered_experiment_all_normalised_AnnData_PHATE_clusterings_sce.rds')
filtered_exp_seurat <- as.Seurat(filtered_exp_sce, data = 'logcounts', counts = NULL)
saveRDS(filtered_exp_sce, 'Prefiltered_experiment_all_normalised_AnnData_PHATE_clusterings_seurat.rds')


Idents(filtered_exp_seurat) <- "Leiden_start_cluster_0.7_opt"
markers <- FindAllMarkers(filtered_exp_seurat, test.use = 'LR', latent.vars = 'Replicate')
write.csv(markers, 'Prefiltered_experiment_all_normalised_AnnData_PHATE_clusterings_seurat_markers.csv')
