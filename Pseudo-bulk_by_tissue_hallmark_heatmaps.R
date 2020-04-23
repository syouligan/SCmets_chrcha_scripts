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

rowMetaData <- data.frame(filtered_exp@assays$RNA@meta.features)

# Load Camera hallmark enrichment results for each tissue
liver_primary_camera <- read.csv("pseudo-bulk_DGE/msigdb/Liver_Primary_Hallmark_camera.csv", header = TRUE)
liver_primary_camera_sig <- liver_primary_camera[liver_primary_camera$FDR < 0.01, ]

ln_primary_camera <- read.csv("pseudo-bulk_DGE/msigdb/LN_Primary_Hallmark_camera.csv", header = TRUE)
ln_primary_camera_sig <- ln_primary_camera[ln_primary_camera$FDR < 0.01, ]

lung_primary_camera <- read.csv("pseudo-bulk_DGE/msigdb/Lung_Primary_Hallmark_camera.csv", header = TRUE)
lung_primary_camera_sig <- lung_primary_camera[lung_primary_camera$FDR < 0.01, ]

hallmark_sig <- unique(c(as.character(liver_primary_camera_sig$X), as.character(ln_primary_camera_sig$X), as.character(lung_primary_camera_sig$X)))

# mSigDB lists
h_df <- msigdbr(species = "Homo sapiens", category = "H")

# Subset to genes active in any tissue
h_df <- h_df[is.element(h_df$gene_symbol, rowMetaData$GeneSymbol), ]

# Add rownames to hallmark lists
idx <- match(h_df$gene_symbol, rowMetaData$GeneSymbol)
h_df$Names <- rownames(rowMetaData) [idx]

# Split hallmark dataframe into lists of Names by geneset
h_list <- h_df %>% split(x = .$Names, f = .$gs_name)

# Make a list of only differentially regulated hallmark datasets 
h_sig_list <- h_list[names(h_list) %in% hallmark_sig]

# Plot heatmap of significant hallmark pathways
DefaultAssay(filtered_exp) <- "SCT"
for(pw in names(h_sig_list)) {
  hm <- DoHeatmap(object = filtered_exp, features = h_sig_list[[pw]], group.by = "Tissue", slot = "scale.data", assay = "SCT")
  ggsave(paste0("pseudo-bulk_DGE/msigdb/", pw, "_heatmap.png"), plot = hm, device = "png")
  }
