#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Make summed experiment for processing locally
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
library('Matrix')
library('limma')
library('edgeR')
library('gtools')
library('msigdbr')
library('org.Hs.eg.db')
library('fgsea')
library('data.table')

set.seed(100)

# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_QC_experiment_practice.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_QC_experiment_all.rds") # uses whole dataset if wolfpack
}

rowData(filtered_exp)$EntrezID <- mapIds(org.Hs.eg.db, keys=as.character(rowData(filtered_exp)$Ensembl), column="ENTREZID", keytype="ENSEMBL", multiVals="first")

# Sum counts across cells from each tissue within each replicate
summed <- aggregateAcrossCells(filtered_exp, ids=DataFrame(label=filtered_exp$Tissue, sample=filtered_exp$Replicate))
exprs <- counts(summed)
colnames(exprs) <- unique(paste0(filtered_exp$Tissue, "_", filtered_exp$Replicate))
exprs <- exprs[ ,order(colnames(exprs))]

# Save Pseudo-bulk dataset for working locally
if(place == "local") {
  saveRDS(exprs, "pseudo-bulk_DGE/Pseudo-bulk_whole_experiment_practice.rds")
} else if(place == "wolfpack") {
  saveRDS(exprs, "pseudo-bulk_DGE/Pseudo-bulk_whole_experiment_all.rds")
} else {
  print("Not overwritten")
}

# Save rowData for working locally
rowDataSummed <- rowData(summed)
if(place == "local") {
  write.csv(rowDataSummed, "pseudo-bulk_DGE/Pseudo-bulk_rowData_practice.csv")
} else if(place == "wolfpack") {
  write.csv(rowDataSummed, "pseudo-bulk_DGE/Pseudo-bulk_rowData_all.csv")
} else {
  print("Not overwritten")
}

