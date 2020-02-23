#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Impute cell expression values using SAVER
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  location <- "/Users/mac/cloudstor/"
  place <- "local"
  folder <- "practice_all_data"
  
} else {
  location <- "/share/ScratchGeneral/scoyou/"
  place <- "wolfpack"
  folder <- "all_data"
}

setwd(paste0(location, "sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/", folder))

# Libraries
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('readr')
library('Matrix')
library('SAVER')

set.seed(100)
# Load prefiltered SingleCellExperiment
if(place == "local") {
  cores <- 3
  filtered_exp <- readRDS("Prefiltered_QC_experiment_practice.rds") # uses practice data if local
} else {
  cores <- 32
  filtered_exp <- readRDS("Prefiltered_QC_experiment_all.rds") # uses whole dataset if wolfpack
}

# Calculate SAVER standard error estimates and plot relative to library size.
if(file.exists("Prefiltered_QC_experiment_SAVER_object.rds")) {
  saver.out <- readRDS("Prefiltered_QC_experiment_SAVER_object.rds")
  } else {
  row_names <- split(rownames(filtered_exp), ceiling(seq_along(rownames(filtered_exp))/2500))
  for(i in names(row_names)){
    saver.new <- saver(counts(filtered_exp), pred.genes = which(rownames(filtered_exp) %in% row_names[[i]]), pred.genes.only = TRUE, do.fast = FALSE, ncores = cores)
    saver.out <- combine.saver(list(saver.out, saver.new))
  }
  saveRDS(saver.out, "Prefiltered_QC_experiment_SAVER_object.rds")
  }

# Plot SAVER estimated values to filtered experiment object
pdf("Filtered_SAVER_SEvsLibSize.pdf")
plot(colSums(counts(filtered_exp)), colMeans(saver.out$se), log="x")
dev.off()

# Add SAVER estimated values to filtered experiment object
assays(filtered_exp, "logestimate") <- log(saver.out$estimate)

# Save total filtered dataset
if (place == "wolfpack") {
  saveRDS(filtered_exp, "Prefiltered_experiment_SAVER_all.rds")
} else {
  saveRDS(filtered_exp, "Prefiltered_experiment_SAVER_practice.rds")
}
