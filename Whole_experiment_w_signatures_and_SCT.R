#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Cluster samples merge samples by replicate within each tissue. Calculate site specificity scores for each cell.
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


# Load prefiltered and clustered Seurat Object
if(place == "local") {
  filtered_exp <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/Prefiltered_experiment_all_seurat_integrated.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Calculate tissue specificity score (correlation) for all cells. Correlation between signature and scaled normalised counts calculated across the whole experiment.
# --------------------------------------------------------------------------

# Load site specificity signatures
if(place == "local") {
  tissue_signatures <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Pseudo-bulk_tissue_specific_signature.rds") # uses practice data if local
} else {
  tissue_signatures <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Pseudo-bulk_tissue_specific_signature.rds") # uses whole dataset if wolfpack
}

# Calculate normalised, scaled counts across the whole experiment
DEG_GOIs <- as.character(unique(unlist(lapply(tissue_signatures, `[[`, "Gene_name"))))
DEG_GOIs <- DEG_GOIs[is.element(DEG_GOIs, rownames(filtered_exp@assays$RNA@counts))]
filtered_exp_SCT <- SCTransform(filtered_exp, verbose = TRUE, vars.to.regress = c("Replicate", "Lib_size", "S.Score", "G2M.Score", "digest_stress1"), variable.features.n = 5000, return.only.var.genes = FALSE, new.assay.name = "SCT_whole")
filtered_exp_magic <- magic(filtered_exp_SCT, assay = "SCT_whole", genes = c(DEG_GOIs), t = "auto")
filtered_exp_magic@active.assay <- "MAGIC_SCT_whole"
filtered_exp_magic <- ScaleData(filtered_exp_magic)

for(sig in names(tissue_signatures)) {
  signature_GOIs <- as.character(tissue_signatures[[sig]]$Gene_name)[is.element(tissue_signatures[[sig]]$Gene_name, rownames(filtered_exp_magic@assays$MAGIC_SCT_whole@scale.data))]
  SCT_signature_GOIs <- Matrix::as.matrix(filtered_exp_magic@assays$MAGIC_SCT_whole@scale.data)[signature_GOIs, ]
  unique(signature_GOIs == rownames(SCT_signature_GOIs))
  sig_tmp <- tissue_signatures[[sig]]
  rownames(sig_tmp) <- sig_tmp$Gene_name
  sig_tmp <- sig_tmp[rownames(SCT_signature_GOIs),]
  correlation <- apply(SCT_signature_GOIs, 2, cor.test, sig_tmp$Signature)
  filtered_exp[[paste0(sig, "_corr")]] <- unlist(lapply(correlation, `[[`, "estimate"))
  filtered_exp[[paste0(sig, "_p.adj")]] <- p.adjust(unlist(lapply(correlation, `[[`, "p.value")), method = "bonf")
}

# Save and remove large objects
if(place == "local") {
  saveRDS(filtered_exp, "Prefiltered_experiment_practice_seurat_integrated_w_SiteSpecSig.rds")
} else if(place == "wolfpack") {
  saveRDS(filtered_exp, "Prefiltered_experiment_all_seurat_integrated_w_SiteSpecSig.rds")
} else {
  print("Not overwritten")
}