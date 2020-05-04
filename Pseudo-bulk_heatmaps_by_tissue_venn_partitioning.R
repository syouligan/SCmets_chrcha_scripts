#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Look at GO terms enriched among DEGs between tissues (pseudobulk) with venn partitioning
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
library('ggplot2')
library('readr')
library('Matrix')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')
library('gplots')
library('msigdbr')

# Load prefiltered and clustered Seurat Object
if(place == "local") {
  filtered_exp <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/Prefiltered_experiment_all_seurat_integrated.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Normalise counts for plotting and save object
filtered_exp <- SCTransform(filtered_exp, verbose = TRUE, vars.to.regress = c("Replicate", "Lib_size", "S.Score", "G2M.Score", "digest_stress1"), variable.features.n = 5000, return.only.var.genes = FALSE, new.assay.name = "SCT_plotting")

if(place == "local") {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_practice_seurat_integrated_plotting.rds")
} else if(place == "wolfpack") {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_all_seurat_integrated_plotting.rds")
} else {
  print("Not overwritten")
}

# Make list of differentially expressed genes
liver_primary <- read.csv("pseudo-bulk_DGE/Liver_Primary_DEG_0LFC.csv", header = TRUE)
liver_primary_DGE_down <- liver_primary[liver_primary$adj.P.Val < 0.05 & liver_primary$logFC > 0, ]
liver_primary_DGE_up <- liver_primary[liver_primary$adj.P.Val < 0.05 & liver_primary$logFC < 0, ]

ln_primary <- read.csv("pseudo-bulk_DGE/LN_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary_DGE_down <- ln_primary[ln_primary$adj.P.Val < 0.05 & ln_primary$logFC > 0, ]
ln_primary_DGE_up <- ln_primary[ln_primary$adj.P.Val < 0.05 & ln_primary$logFC < 0, ]

lung_primary <- read.csv("pseudo-bulk_DGE/Lung_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary_DGE_down <- lung_primary[lung_primary$adj.P.Val < 0.05 & lung_primary$logFC > 0, ]
lung_primary_DGE_up <- lung_primary[lung_primary$adj.P.Val < 0.05 & lung_primary$logFC < 0, ]

DGE_list_Names_down <- list("Liver_down" = liver_primary_DGE_down$X, "LN_down" = ln_primary_DGE_down$X, "Lung_down" = lung_primary_DGE_down$X)
DGE_list_Names_up <- list("Liver_up" = liver_primary_DGE_down$X, "LN_up" = ln_primary_DGE_down$X, "Lung_up" = lung_primary_DGE_up$X)

# Heatmap of down regulated DEGs grouped by tissue and partitioned
# --------------------------------------------------------------------------

# Venn diagram overlap DEGs between metastatic sites and tissues
pdf("pseudo-bulk_DGE/DGE_psuedo-bulk_venn_down.pdf")
ItemsList <- venn(DGE_list_Names_down, show.plot = TRUE)
dev.off()

lists <- attributes(ItemsList)$intersections
hm <- DoHeatmap(object = filtered_exp, features = as.character(unlist(lists)), group.by = "Tissue", slot = "scale.data", assay = "SCT_plotting")
ggsave("pseudo-bulk_DGE/DGE_psuedo-bulk_venn_down_heatmap.png", plot = hm, device = "png")

# List of DEGs partitioned by their venn segment
DGE_list <- list("LiverOnly" = liver_primary[is.element(liver_primary$X, lists[["Liver"]]), ], "LungOnly" = lung_primary[is.element(lung_primary$X, lists[["Lung"]]), ], "All" = ln_primary[is.element(ln_primary$X, lists[["Liver:LN:Lung"]]), ], "LNOnly" = ln_primary[is.element(ln_primary$X, lists[["LN"]]), ], "Liver&Lung" = liver_primary[is.element(liver_primary$X, lists[["Liver:Lung"]]), ], "Liver&LN" = liver_primary[is.element(liver_primary$X, lists[["Liver:LN"]]), ], "LN&Lung" = ln_primary[is.element(ln_primary$X, lists[["LN:Lung"]]), ])

# Heatmap of up regulated DEGs grouped by tissue and partitioned
# --------------------------------------------------------------------------

# Venn diagram overlap DEGs between metastatic sites and tissues
pdf("pseudo-bulk_DGE/DGE_psuedo-bulk_venn_up.pdf")
ItemsList <- venn(DGE_list_Names_up, show.plot = TRUE)
dev.off()

lists <- attributes(ItemsList)$intersections
hm <- DoHeatmap(object = filtered_exp, features = as.character(unlist(lists)), group.by = "Tissue", slot = "scale.data", assay = "SCT_plotting")
ggsave("pseudo-bulk_DGE/DGE_psuedo-bulk_venn_up_heatmap.png", plot = hm, device = "png")

# List of DEGs partitioned by their venn segment
DGE_list <- list("LiverOnly" = liver_primary[is.element(liver_primary$X, lists[["Liver"]]), ], "LungOnly" = lung_primary[is.element(lung_primary$X, lists[["Lung"]]), ], "All" = ln_primary[is.element(ln_primary$X, lists[["Liver:LN:Lung"]]), ], "LNOnly" = ln_primary[is.element(ln_primary$X, lists[["LN"]]), ], "Liver&Lung" = liver_primary[is.element(liver_primary$X, lists[["Liver:Lung"]]), ], "Liver&LN" = liver_primary[is.element(liver_primary$X, lists[["Liver:LN"]]), ], "LN&Lung" = ln_primary[is.element(ln_primary$X, lists[["LN:Lung"]]), ])

