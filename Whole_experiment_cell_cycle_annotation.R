#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Annotate cells with cell cycle information using Cyclone ('scran')
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

# Assign cell cycle phases using Cyclone
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
cyc.values <- cyclone(filtered_exp, hs.pairs, gene.names = rowData(filtered_exp)$Ensembl)
filtered_exp$CC.Phase <- cyc.values$phases
filtered_exp$CC.G1 <- cyc.values$normalized.scores$G1
filtered_exp$CC.S <- cyc.values$normalized.scores$S
filtered_exp$CC.G2M <- cyc.values$normalized.scores$G2M

# Plot clusters annotated with CC phase
plotReducedDim(filtered_exp, dimred="corrected_logestimate", colour_by = "CC.Phase", text_by = "cluster") +
  scale_fill_viridis_d(option = "D") +
  ggsave("Fastmnn_corrected_with_clusters_cell.cycle.phase.pdf")
plotReducedDim(filtered_exp, dimred="UMAP_fastMNN", colour_by = "CC.Phase", text_by = "cluster") +
  scale_fill_viridis_d(option = "D") +
  ggsave("UMAP_corrected_with_clusters_cell.cycle.phase.pdf")
plotReducedDim(filtered_exp, dimred="PHATE_fastMNN", colour_by = "CC.Phase", text_by = "cluster") +
  scale_fill_viridis_d(option = "D") +
  ggsave("PHATE_corrected_with_clusters_cell.cycle.phase.pdf")

# Save SCE object with CC annotations
if(place == "local") {
  saveRDS(filtered_exp, "Prefiltered_experiment_practice_merge_cluster_wCC.rds") # saves practice data if local
} else {
  saveRDS(filtered_exp, "Prefiltered_experiment_all_merge_cluster_wCC.rds") # saves whole dataset if wolfpack
}
