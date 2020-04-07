#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells using Seurat within each tissue
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_by_organ") # Uses practice data (5% of cells from each sample) if running locally
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

# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/practice_all_data/Prefiltered_QC_experiment_practice.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/Prefiltered_QC_experiment_all.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Subset and integrate by tissue
for(tissue in unique(filtered_exp$Tissue)) {
  print(tissue)
  tissue_exp <- filtered_exp[,filtered_exp$Tissue == tissue]
  dir.create(tissue)
  
# Set up prefiltered object, add cell cycle difference info and split by sample
# --------------------------------------------------------------------------

tissue_exp_seurat <- as.Seurat(tissue_exp, counts = "counts", data = NULL) # convert to Seurat

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
tissue_exp_seurat <- CellCycleScoring(tissue_exp_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
tissue_exp_seurat$CC.Difference <- tissue_exp_seurat$S.Score - tissue_exp_seurat$G2M.Score

tissue_exp.list <- SplitObject(tissue_exp_seurat, split.by = "Sample") # split into individual samples

# Normalise transform counts within each experiment. Note: will not overwrite if already exists.
# --------------------------------------------------------------------------

# Integrate datasets based on highly correlated features
if(file.exists(paste0(tissue, "/Prefiltered_experiment_all_seurat_integrated_", tissue,".rds")) & place == "wolfpack") {
  tissue_exp.integrated <- readRDS(paste0(tissue, "/Prefiltered_experiment_all_seurat_integrated_", tissue,".rds"))
} else if(file.exists(paste0(tissue, "/Prefiltered_experiment_practice_seurat_integrated_", tissue,".rds")) & place == "local") {
  tissue_exp.integrated <- readRDS(paste0(tissue, "/Prefiltered_experiment_practice_seurat_integrated_", tissue,".rds"))
} else {
  # Perform SCT normalisation on each dataset individually
  for (i in 1:length(tissue_exp.list)) {
    tissue_exp.list[[i]] <- SCTransform(tissue_exp.list[[i]], verbose = TRUE, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score", "Mito_percent"))
  }
  # Integrate  
  tissue_exp.features <- SelectIntegrationFeatures(object.list = tissue_exp.list, nfeatures = 5000)
  tissue_exp.list <- PrepSCTIntegration(object.list = tissue_exp.list, verbose = TRUE, anchor.features = tissue_exp.features)
  # tissue_exp.list <- lapply(X = tissue_exp.list, FUN = RunPCA, verbose = TRUE, features = tissue_exp.features) # Perform PCA on each object individually (needed for rpca)
  tissue_exp.anchors <- FindIntegrationAnchors(object.list = tissue_exp.list, normalization.method = "SCT", anchor.features = tissue_exp.features, verbose = TRUE, reduction = "cca")
  tissue_exp.integrated <- IntegrateData(anchorset = tissue_exp.anchors, normalization.method = "SCT", verbose = TRUE)
}

# Run PCA on intergated dataset and determine clusters using Seurat
# --------------------------------------------------------------------------

# Seurat integration pipeline
tissue_exp.integrated <- RunPCA(tissue_exp.integrated, dims = 1:50, assay = "integrated", ndims.print = 1:5, nfeatures.print = 5)
tissue_exp.integrated <- FindNeighbors(tissue_exp.integrated, reduction = "pca", dims = 1:50)
tissue_exp.integrated <- FindClusters(tissue_exp.integrated, resolution = 0.05)
tissue_exp.integrated <- RunUMAP(tissue_exp.integrated, reduction = "pca", dims = 1:50)
tissue_exp.integrated[["seurat_PCA_clusters"]] <- Idents(object = tissue_exp.integrated)

for(i in c("pca", "umap")) {
  p2 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Phase")
  p5 <- DimPlot(tissue_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p2, p3, p5)
  ggsave(paste0(tissue, "/Seurat_clusters_", tissue, "_", i, ".png", plot = gridit), device = "png")
}

# Run PHATE on intergated dataset and determine clusters using kmeans
# --------------------------------------------------------------------------

# Run phate, find number of clusters in embeddings using silhouette
phate.out <- phate(Matrix::t(GetAssayData(tissue_exp.integrated, assay = "integrated", slot = "data")), ndim = 10) # Runs PHATE diffusion map
tissue_exp.integrated[["phate"]] <- CreateDimReducObject(embeddings = phate.out$embedding, key = "PHATE_", assay = "integrated")
ProjectDim(tissue_exp.integrated, reduction = "phate")
tissue_exp.integrated <- FindNeighbors(tissue_exp.integrated, reduction = "phate", dims = 1:10)
tissue_exp.integrated <- FindClusters(tissue_exp.integrated, resolution = 0.05)

for(i in c("pca", "umap", "phate")) {
  p2 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Phase")
  p5 <- DimPlot(tissue_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p2, p3, p5)
  ggsave(paste0(tissue, "/PHATE_clusters_", tissue, "_", i, ".png", plot = gridit), device = "png")
}

# Save seurat objects
# --------------------------------------------------------------------------

if(place == "local") {
  # saveRDS(phate.out, paste0(tissue, "/Prefiltered_experiment_practice_phate.out_", tissue,".rds"))
} else {
  saveRDS(phate.out, paste0(tissue, "/Prefiltered_experiment_all_phate.out_", tissue,".rds"))
}

if(place == "local") {
  # saveRDS(tissue_exp.list, paste0(tissue, "/Prefiltered_experiment_practice_seurat_list_", tissue,".rds"))
} else {
  saveRDS(tissue_exp.list, paste0(tissue, "/Prefiltered_experiment_all_seurat_list_", tissue,".rds"))
}

if(place == "local") {
  # saveRDS(tissue_exp.anchors, paste0(tissue, "/Prefiltered_experiment_practice_seurat_anchors_", tissue,".rds"))
} else {
  saveRDS(tissue_exp.anchors, paste0(tissue, "/Prefiltered_experiment_all_seurat_anchors_", tissue,".rds"))
}

if(place == "local") {
  saveRDS(tissue_exp.integrated, paste0(tissue, "/Prefiltered_experiment_practice_seurat_integrated_", tissue,".rds"))
} else {
  saveRDS(tissue_exp.integrated, paste0(tissue, "/Prefiltered_experiment_all_seurat_integrated_", tissue,".rds"))
}

}
