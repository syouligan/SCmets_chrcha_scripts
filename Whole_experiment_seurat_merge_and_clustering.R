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


# Load prefiltered and clustered Seurat Object
if(place == "local") {
  filtered_exp_seurat <- readRDS("Prefiltered_experiment_practice_seurat_tissueDGE.rds") # uses practice data if local
} else {
  filtered_exp_seurat <- readRDS("Prefiltered_experiment_all_seurat_tissueDGE.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Add score to idenitfy dying cells with low mito content (as per DOI 10.1186/s13059-019-1830-0)
dying <- list(c("HLA-A", "HLA-B", "HLA-C", "B2M"))
filtered_exp_seurat <- AddModuleScore(object = filtered_exp_seurat, features = dying, name = 'dying')
filtered_exp.list <- SplitObject(filtered_exp_seurat, split.by = "Replicate") # split into individual samples

# Normalise transform counts within each experiment. Note: will not overwrite if already exists.
# --------------------------------------------------------------------------

# Perform SCT normalisation on each dataset individually
for (i in 1:length(filtered_exp.list)) {
  filtered_exp.list[[i]] <- SCTransform(filtered_exp.list[[i]], verbose = TRUE, vars.to.regress = c("Lib_size", "S.Score", "G2M.Score", "digest_stress1"))
}

# Integrate datasets based on highly correlated features
filtered_exp.features <- SelectIntegrationFeatures(object.list = filtered_exp.list, nfeatures = 5000)
filtered_exp.list <- PrepSCTIntegration(object.list = filtered_exp.list, verbose = TRUE, anchor.features = filtered_exp.features)
# filtered_exp.list <- lapply(X = filtered_exp.list, FUN = RunPCA, verbose = TRUE, features = filtered_exp.features) # Perform PCA on each object individually (needed for rpca)
reference_datasets <- which(names(filtered_exp.list) == "3")
filtered_exp.anchors <- FindIntegrationAnchors(object.list = filtered_exp.list, normalization.method = "SCT", anchor.features = filtered_exp.features, verbose = TRUE, reduction = "cca", reference = reference_datasets)
filtered_exp.integrated <- IntegrateData(anchorset = filtered_exp.anchors, normalization.method = "SCT", verbose = TRUE)

# Run PCA on intergated dataset and determine clusters using Seurat
# --------------------------------------------------------------------------

# Seurat integration pipeline
filtered_exp.integrated <- RunPCA(filtered_exp.integrated, dims = 1:50, assay = "integrated", ndims.print = 1:5, nfeatures.print = 5)
filtered_exp.integrated <- FindNeighbors(filtered_exp.integrated, reduction = "pca", dims = 1:50)
filtered_exp.integrated <- FindClusters(filtered_exp.integrated, resolution = 0.01)
filtered_exp.integrated <- RunUMAP(filtered_exp.integrated, reduction = "pca", dims = 1:50)
filtered_exp.integrated[["seurat_PCA_clusters"]] <- Idents(object = filtered_exp.integrated)

for(i in c("pca", "umap")) {
  p1 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Tissue")
  p2 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- FeaturePlot(filtered_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
  p4 <- FeaturePlot(filtered_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
  p5 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Phase")
  p6 <- DimPlot(filtered_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
  ggsave(paste0("Seurat_clusters_", i, ".png", plot = gridit), device = "png")
}

# Run PHATE on intergated dataset and determine clusters using kmeans
# --------------------------------------------------------------------------

# Run phate, find number of clusters in embeddings using silhouette
phate.out <- phate(Matrix::t(GetAssayData(filtered_exp.integrated, assay = "integrated", slot = "data")), ndim = 10) # Runs PHATE diffusion map
filtered_exp.integrated[["phate"]] <- CreateDimReducObject(embeddings = phate.out$embedding, key = "PHATE_", assay = "integrated")
ProjectDim(filtered_exp.integrated, reduction = "phate")
filtered_exp.integrated <- FindNeighbors(filtered_exp.integrated, reduction = "phate", dims = 1:10)
filtered_exp.integrated <- FindClusters(filtered_exp.integrated, resolution = 0.01)
phate_embed <- data.frame(Embeddings(filtered_exp.integrated, reduction = "phate"))

for(i in c("pca", "umap")) {
  p1 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Tissue")
  p2 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- FeaturePlot(filtered_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
  p4 <- FeaturePlot(filtered_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
  p5 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Phase")
  p6 <- DimPlot(filtered_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
  ggsave(paste0("PHATE_clusters_", i, ".png", plot = gridit), device = "png")
}

# No idea why phate needs coordinates set
for(i in c("phate")) {
  p1 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Tissue")
  p2 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- FeaturePlot(filtered_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
  p3 <- p3 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
  p4 <- FeaturePlot(filtered_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
  p4 <- p4 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
  p5 <- DimPlot(filtered_exp.integrated, reduction = i, group.by = "Phase")
  p6 <- DimPlot(filtered_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
  ggsave(paste0("PHATE_clusters_", i, ".png", plot = gridit), device = "png")
}

# Save seurat objects
# --------------------------------------------------------------------------

if(place == "local" & exists("phate.out")) {
  # saveRDS(phate.out, "Prefiltered_experiment_practice_phate.out.rds")
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(phate.out, "Prefiltered_experiment_all_phate.out.rds")
} else {
  print("Not overwritten")
}

if(place == "local" & exists("phate.out")) {
  # saveRDS(filtered_exp.list, "Prefiltered_experiment_practice_seurat_list.rds")
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(filtered_exp.list, "Prefiltered_experiment_all_seurat_list.rds")
} else {
  print("Not overwritten")
}

if(place == "local" & exists("phate.out")) {
  # saveRDS(filtered_exp.anchors, "Prefiltered_experiment_practice_seurat_anchors.rds")
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(filtered_exp.anchors, "Prefiltered_experiment_all_seurat_anchors.rds")
} else {
  print("Not overwritten")
}

if(place == "local" & exists("phate.out")) {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_practice_seurat_tissueDGE_integrated.rds")
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(filtered_exp.integrated, "Prefiltered_experiment_all_seurat_tissueDGE_integrated.rds")
} else {
  print("Not overwritten")
}
