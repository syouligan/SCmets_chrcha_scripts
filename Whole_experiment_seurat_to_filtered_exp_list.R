#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells using Seurat.
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data")
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

# Set up prefiltered object, add cell cycle difference info and split by sample
# --------------------------------------------------------------------------
filtered_exp_seurat <- as.Seurat(filtered_exp_sce, counts = "counts", data = NULL) # convert to Seurat

# Add Entrez IDs
rowData(filtered_exp_sce)$EntrezID <- mapIds(org.Hs.eg.db, keys=rowData(filtered_exp_sce)$Ensembl, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
rowMetaData <- data.frame(rowData(filtered_exp_sce))
filtered_exp_seurat <- AddMetaData(filtered_exp_seurat, rowMetaData)

# Add cell cycle
s.genes <- cc.genes.updated.2019$s.genes
s.genes <- c(rownames(rowMetaData[is.element(rowMetaData$GeneSymbol, s.genes), ]))
g2m.genes <- cc.genes.updated.2019$g2m.genes
g2m.genes <- c(rownames(rowMetaData[is.element(rowMetaData$GeneSymbol, g2m.genes), ]))
filtered_exp_seurat <- CellCycleScoring(filtered_exp_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

# Add score to idenitfy dying cells and stress scores with low mito content (as per DOI 10.1186/s13059-019-1830-0)
digest_stress <- c("FOS", "CXCL2", "ZFP36", "FOSB", "DUSP1", "ATF3", "CXCL8", "NR4A1", "CXCL3", "PPP1R15A", "JUNB", "EGR1", "HSPA1A", "HSPA1B", "SOCS3", "KLF6", "JUN", "IER2", "CXCL1", "NKFBIA", "HSPA6", "DNAJB1", "IER3", "CCNL1", "MTRNR2L2", "IER5", "ID1", "CEBPD", "KRT6A", "CYR61", "DEPP1", "CLDN4", "IRF1", "DUSP2", "BTG2", "PLAUR", "MAFF", "KLF4", "PHLDA2", "TNFAIP3")
digest_stress <- list(c(rownames(rowMetaData[is.element(rowMetaData$GeneSymbol, digest_stress), ])))
filtered_exp_seurat <- AddModuleScore(object = filtered_exp_seurat, features = digest_stress, name = 'digest_stress')

dying <- c("HLA-A", "HLA-B", "HLA-C", "B2M")
dying <- list(c(rownames(rowMetaData[is.element(rowMetaData$GeneSymbol, dying), ])))
filtered_exp_seurat <- AddModuleScore(object = filtered_exp_seurat, features = dying, name = 'dying')

# Normalise transform counts within each experiment. Note: will not overwrite if already exists.
# --------------------------------------------------------------------------

filtered_exp.list <- SplitObject(filtered_exp_seurat, split.by = "Replicate") # split into individual samples
all_features <- data.frame(rownames(x = filtered_exp_seurat))

# Perform SCT normalisation on each dataset individually
for (i in 1:length(filtered_exp.list)) {
  filtered_exp.list[[i]] <- SCTransform(filtered_exp.list[[i]], verbose = TRUE, vars.to.regress = c("Lib_size", "S.Score", "G2M.Score", "digest_stress1"), variable.features.n = 5000, return.only.var.genes = FALSE)
}

# Save seurat objects
# --------------------------------------------------------------------------

if(place == "local") {
  saveRDS(filtered_exp.list, "Prefiltered_QC_experiment_practice_filtered_exp_list.rds")
  saveRDS(all_features, "All_features.rds")
} else if(place == "wolfpack") {
  saveRDS(filtered_exp.list, "Prefiltered_QC_experiment_all_filtered_exp_list.rds")
  saveRDS(all_features, "All_features.rds")
} else {
  print("Not saved")
}
