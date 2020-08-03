#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells using Seurat.
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
  filtered_exp_sce <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/practice_all_data/Prefiltered_QC_experiment_practice.rds") # uses practice data if local
} else {
  filtered_exp_sce <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/Prefiltered_QC_experiment_all.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Save counts, rowData and colData for filtered counts
if(place == "local") {
  writeMM(counts(filtered_exp_sce), "Prefiltered_experiment_practice_filtered_counts.mtx")
  write.csv(colData(filtered_exp_sce), "Prefiltered_experiment_practice_colData.csv")
  write.csv(rowData(filtered_exp_sce), "Prefiltered_experiment_practice_rowData.csv")
} else {
  writeMM(counts(filtered_exp_sce), "Prefiltered_experiment_all_filtered_counts.mtx")
  write.csv(colData(filtered_exp_sce), "Prefiltered_experiment_all_colData.csv")
  write.csv(rowData(filtered_exp_sce), "Prefiltered_experiment_all_rowData.csv")
}
