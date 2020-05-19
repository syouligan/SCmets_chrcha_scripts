#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Make imputed heatmaps of HALLMARK pathways.
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_by_organ/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/by_organ")
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
library('viridis')
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
library('Rmagic')
library('msigdbr')
library('colorblindr')
library('ggridges')

# Load prefiltered and clustered Seurat Object
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_experiment_practice_seurat_integrated_w_SiteSpecSig.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_experiment_all_seurat_integrated_w_SiteSpecSig.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

if(place == "local") {
  filtered_exp_SCT <- readRDS("Prefiltered_experiment_practice_seurat_integrated_whole_SCT.rds") # uses practice data if local
} else {
  filtered_exp_SCT <- readRDS("Prefiltered_experiment_all_seurat_integrated_whole_SCT.rds") # uses whole dataset if wolfpack
}

dir.create("plots")
rowMetaData <- data.frame(filtered_exp@assays$RNA@meta.features)
Idents(filtered_exp) <- "Tissue"
levels(filtered_exp) <- c("Primary", "LN", "Liver", "Lung")


# Heatmaps of significant HALLMARK pathways
# --------------------------------------------------------------------------

# Load Camera hallmark enrichment results for each tissue
if(place == "local") {
  liver_primary_camera <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Liver_Primary_Hallmark_camera.csv", header = TRUE)
} else {
  liver_primary_camera <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Liver_Primary_Hallmark_camera.csv", header = TRUE)
}
liver_primary_camera_sig <- liver_primary_camera[liver_primary_camera$FDR < 0.01, ]

if(place == "local") {
  ln_primary_camera <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/LN_Primary_Hallmark_camera.csv", header = TRUE)
} else {
  ln_primary_camera <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/LN_Primary_Hallmark_camera.csv", header = TRUE)
}
ln_primary_camera_sig <- ln_primary_camera[ln_primary_camera$FDR < 0.01, ]

if(place == "local") {
  lung_primary_camera <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Lung_Primary_Hallmark_camera.csv", header = TRUE)
} else {
  lung_primary_camera <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/msigdb/Lung_Primary_Hallmark_camera.csv", header = TRUE)
}
lung_primary_camera_sig <- lung_primary_camera[lung_primary_camera$FDR < 0.01, ]

hallmark_sig <- unique(c(as.character(liver_primary_camera_sig$X), as.character(ln_primary_camera_sig$X), as.character(lung_primary_camera_sig$X)))

# mSigDB lists
h_df <- msigdbr(species = "Homo sapiens", category = "H")

# Subset to genes active in any tissue
h_df <- h_df[is.element(h_df$gene_symbol, rowMetaData$GeneSymbol), ]

# Add rownames to hallmark lists
idx <- match(h_df$gene_symbol, rowMetaData$GeneSymbol)
h_df$Names <- rowMetaData$Row.names [idx]

# Split hallmark dataframe into lists of Names by geneset
h_list <- h_df %>% split(x = .$Names, f = .$gs_name)

# Make a list of only differentially regulated hallmark datasets 
h_sig_list <- h_list[names(h_list) %in% hallmark_sig]

# Impute counts for genes to be plotted
HM_GOIs <- as.character(unique(h_df$Names))
filtered_exp_magic <- magic(filtered_exp_SCT, assay = "SCT_whole", genes = c(HM_GOIs), t = "auto")
filtered_exp_magic@active.assay <- 'MAGIC_SCT_whole'
filtered_exp_magic <- ScaleData(filtered_exp_magic)
Idents(filtered_exp_magic) <- "Tissue"
levels(filtered_exp_magic) <- c("Primary", "LN", "Liver", "Lung")

for(pw in names(h_sig_list)) {
  hm <- DoHeatmap(object = filtered_exp_magic, features = h_sig_list[[pw]], slot = "scale.data", assay = "MAGIC_SCT_whole", group.bar = TRUE, group.colors = palette_OkabeIto[1:4], label = FALSE) +
    scale_fill_gradientn(colors = hcl.colors(n = 7, palette = "Berlin"), limits = c(-2.5, 2.5))
  ggsave(paste0("plots/", pw, "_heatmap.png"), plot = hm, device = "png")
}