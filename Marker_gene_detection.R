#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Find and save markers for each cluster
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
library('scater')
library('scran')
library('ggplot2')
library('readr')
library('Matrix')

set.seed(100)
# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_experiment_practice_merge_cluster.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_experiment_all_merge_cluster.rds") # uses whole dataset if wolfpack
}

# Find markers between each cluster and save a csv
markers.filtered_exp <- findMarkers(filtered_exp, test="wilcox", filtered_exp$cluster, direction="up", block = filtered_exp$Sample, BPPARAM  = MulticoreParam(), pval.type = "some")
for(i in names(markers.filtered_exp)){
  interesting <- markers.filtered_exp[[i]]
  write.csv(interesting, paste0("markers/Cluster_", i, "_markers.csv"))
}



