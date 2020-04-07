#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells using Seurat
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


# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/practice_all_data/Prefiltered_QC_experiment_practice.rds") # uses practice data if local
  } else {
  filtered_exp <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/Prefiltered_QC_experiment_all.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Set up prefiltered object, add cell cycle difference info and split by sample
# --------------------------------------------------------------------------
filtered_exp_seurat <- as.Seurat(filtered_exp, counts = "counts", data = NULL) # convert to Seurat

s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
filtered_exp_seurat <- CellCycleScoring(filtered_exp_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
filtered_exp_seurat$CC.Difference <- filtered_exp_seurat$S.Score - filtered_exp_seurat$G2M.Score

filtered_exp.list <- SplitObject(filtered_exp_seurat, split.by = "Sample") # split into individual samples

# Normalise transform counts within each experiment. Note: will not overwrite if already exists.
# --------------------------------------------------------------------------

# Integrate datasets based on highly correlated features
if(file.exists("Prefiltered_experiment_all_seurat_integrated.rds") & place == "wolfpack") {
  filtered_exp.integrated <- readRDS("Prefiltered_experiment_all_seurat_integrated.rds")
  } else if(file.exists("Prefiltered_experiment_practice_seurat_integrated.rds") & place == "local") {
  filtered_exp.integrated <- readRDS("Prefiltered_experiment_practice_seurat_integrated.rds")
  } else {
    # Perform SCT normalisation on each dataset individually
    for (i in 1:length(filtered_exp.list)) {
      filtered_exp.list[[i]] <- SCTransform(filtered_exp.list[[i]], verbose = TRUE, vars.to.regress = c("nCount_RNA", "S.Score", "G2M.Score", "Mito_percent"))
    }
  # Integrate  
  filtered_exp.features <- SelectIntegrationFeatures(object.list = filtered_exp.list, nfeatures = 5000)
  filtered_exp.list <- PrepSCTIntegration(object.list = filtered_exp.list, verbose = TRUE, anchor.features = filtered_exp.features)
  # filtered_exp.list <- lapply(X = filtered_exp.list, FUN = RunPCA, verbose = TRUE, features = filtered_exp.features) # Perform PCA on each object individually (needed for rpca)
  reference_datasets <- which(names(filtered_exp.list) == "LN_3" | names(filtered_exp.list) == "Liver_3" | names(filtered_exp.list) == "Lung_3"  | names(filtered_exp.list) == "Primary_3")
  filtered_exp.anchors <- FindIntegrationAnchors(object.list = filtered_exp.list, normalization.method = "SCT", anchor.features = filtered_exp.features, verbose = TRUE, reduction = "cca", reference = reference_datasets)
  filtered_exp.integrated <- IntegrateData(anchorset = filtered_exp.anchors, normalization.method = "SCT", verbose = TRUE)
  }

# Run PCA on intergated dataset and determine clusters using Seurat
# --------------------------------------------------------------------------

# Seurat integration pipeline
filtered_exp.integrated <- RunPCA(filtered_exp.integrated, dims = 1:50, assay = "integrated", ndims.print = 1:5, nfeatures.print = 5)
filtered_exp.integrated <- FindNeighbors(filtered_exp.integrated, reduction = "pca", dims = 1:50)
filtered_exp.integrated <- FindClusters(filtered_exp.integrated, resolution = 0.05)
filtered_exp.integrated <- RunUMAP(filtered_exp.integrated, reduction = "pca", dims = 1:50)
filtered_exp.integrated[["seurat_PCA_clusters"]] <- Idents(object = filtered_exp.integrated)

for(i in c("pca", "umap")) {
  p1 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Tissue")
  p2 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Phase")
  p5 <- DimPlot(filtered_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p1, p2, p3, p5)
  ggsave(paste0("Seurat_clusters_", i, ".png", plot = gridit), device = "png")
}

# Run PHATE on intergated dataset and determine clusters using kmeans
# --------------------------------------------------------------------------

# Run phate, find number of clusters in embeddings using silhouette
phate.out <- phate(Matrix::t(GetAssayData(filtered_exp.integrated, assay = "integrated", slot = "data")), ndim = 10) # Runs PHATE diffusion map
filtered_exp.integrated[["phate"]] <- CreateDimReducObject(embeddings = phate.out$embedding, key = "PHATE_", assay = "integrated")
ProjectDim(filtered_exp.integrated, reduction = "phate")
filtered_exp.integrated <- FindNeighbors(filtered_exp.integrated, reduction = "phate", dims = 1:10)
filtered_exp.integrated <- FindClusters(filtered_exp.integrated, resolution = 0.05)

for(i in c("pca", "umap", "phate")) {
  p1 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Tissue")
  p2 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Phase")
  p5 <- DimPlot(filtered_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p1, p2, p3, p5)
  ggsave(paste0("PHATE_clusters_", i, ".png", plot = gridit), device = "png")
}

# Save seurat objects
# --------------------------------------------------------------------------

if(place == "local") {
  # saveRDS(phate.out, "Prefiltered_experiment_practice_phate.out.rds")
} else {
  saveRDS(phate.out, "Prefiltered_experiment_all_phate.out.rds")
}

if(place == "local") {
  # saveRDS(filtered_exp.list, "Prefiltered_experiment_practice_seurat_list.rds")
} else {
  saveRDS(filtered_exp.list, "Prefiltered_experiment_all_seurat_list.rds")
}

if(place == "local") {
  # saveRDS(filtered_exp.anchors, "Prefiltered_experiment_practice_seurat_anchors.rds")
} else {
  saveRDS(filtered_exp.anchors, "Prefiltered_experiment_all_seurat_anchors.rds")
}

if(place == "local") {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_practice_seurat_integrated.rds")
} else {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_all_seurat_integrated.rds")
}
