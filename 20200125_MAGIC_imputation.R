#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Generate and save imputed counts object using MAGIC. To be used for plotting
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCmets_chrcha/project_results/prefiltered/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCmets_chrcha/project_results/prefiltered/all_data")
  place <- "wolfpack"
}

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('phateR')
library('ggplot2')
library('readr')
library('Rmagic')
library('Matrix')

set.seed(100)
# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_experiment_Practice_merge_cluster_wCC.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_experiment_All_merge_cluster_wCC.rds") # uses whole dataset if wolfpack
}

# Perform imputation on logcounts
counts_MAGIC <- magic(Matrix::t(assay(filtered_exp, "logcounts")), genes = "all_genes")
assay(filtered_exp, "logcounts_MAGIC") <- Matrix(Matrix::t(counts_MAGIC$result), sparse = TRUE)

# Perform imputation on fastMNN corrected counts
counts_MAGIC <- magic(Matrix::t(assay(filtered_exp, "reconstructed_fastMNN")), genes = "all_genes")
assay(filtered_exp, "reconstructed_fastMNN_MAGIC") <- Matrix(Matrix::t(counts_MAGIC$result), sparse = TRUE)

# Save SCE object with imputed counts
if(place == "local") {
  saveRDS(filtered_exp, "Prefiltered_experiment_Practice_merge_cluster_wCC_MAGIC.rds") # saves practice data if local
} else {
  saveRDS(filtered_exp, "Prefiltered_experiment_All_merge_cluster_wCC_MAGIC.rds") # saves whole dataset if wolfpack
}
