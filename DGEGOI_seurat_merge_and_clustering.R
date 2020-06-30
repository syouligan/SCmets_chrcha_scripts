#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells using Seurat.
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/DGEGOI/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/DGEGOI/")
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

# Load prefiltered and clustered Seurat Object
if(place == "local") {
  filtered_exp <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/Prefiltered_experiment_all_seurat_integrated.rds") # uses whole dataset if wolfpack
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Load site specificity signatures
if(place == "local") {
  tissue_signatures <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Pseudo-bulk_tissue_specific_signature.rds") # uses practice data if local
} else {
  tissue_signatures <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Pseudo-bulk_tissue_specific_signature.rds") # uses whole dataset if wolfpack
}

# Calculate normalised, scaled counts across the whole experiment
DEG_GOIs <- as.character(unique(unlist(lapply(tissue_signatures, `[[`, "Gene_name"))))
DEG_GOIs <- DEG_GOIs[is.element(DEG_GOIs, rownames(filtered_exp@assays$RNA@counts))]

# Set up prefiltered object, add cell cycle difference info and split by sample
# --------------------------------------------------------------------------

DefaultAssay(filtered_exp) <- "RNA"
DGEGOI_exp <- subset(filtered_exp, features = DEG_GOIs)
n.features <- nrow(DGEGOI_exp)
DGEGOI_exp[['integrated']] <- NULL
DGEGOI_exp[['SCT']] <- NULL

DGEGOI_exp$Liver <- DGEGOI_exp$Tissue == "Liver"
DGEGOI_exp$LN <- DGEGOI_exp$Tissue == "LN"
DGEGOI_exp$Lung <- DGEGOI_exp$Tissue == "Lung"
DGEGOI_exp$Primary <- DGEGOI_exp$Tissue == "Primary"

# Normalise transform counts within each experiment. Note: will not overwrite if already exists.
# --------------------------------------------------------------------------

DGEGOI_exp.list <- SplitObject(DGEGOI_exp, split.by = "Replicate") # split into individual samples

# Perform SCT normalisation on each dataset individually
for (i in 1:length(DGEGOI_exp.list)) {
  DGEGOI_exp.list[[i]] <- SCTransform(DGEGOI_exp.list[[i]], verbose = TRUE, vars.to.regress = c("Lib_size", "S.Score", "G2M.Score", "digest_stress1"), variable.features.n = n.features, return.only.var.genes = FALSE)
}

# Integrate datasets based on highly correlated features
DGEGOI_exp.features <- SelectIntegrationFeatures(object.list = DGEGOI_exp.list, nfeatures = n.features)
DGEGOI_exp.list <- PrepSCTIntegration(object.list = DGEGOI_exp.list, verbose = TRUE, anchor.features = DGEGOI_exp.features)
# DGEGOI_exp.list <- lapply(X = DGEGOI_exp.list, FUN = RunPCA, verbose = TRUE, features = DGEGOI_exp.features) # Perform PCA on each object individually (needed for rpca)
reference_datasets <- which(names(DGEGOI_exp.list) == "3")
DGEGOI_exp.anchors <- FindIntegrationAnchors(object.list = DGEGOI_exp.list, normalization.method = "SCT", anchor.features = DGEGOI_exp.features, verbose = TRUE, reduction = "cca", reference = reference_datasets)
DGEGOI_exp.integrated <- IntegrateData(anchorset = DGEGOI_exp.anchors, normalization.method = "SCT", verbose = TRUE)

# Run PCA on intergated dataset and determine clusters using Seurat
# --------------------------------------------------------------------------

# Seurat integration pipeline
DGEGOI_exp.integrated <- RunPCA(DGEGOI_exp.integrated, dims = 1:50, assay = "integrated", ndims.print = 1:5, nfeatures.print = 5)
DGEGOI_exp.integrated <- FindNeighbors(DGEGOI_exp.integrated, reduction = "pca", dims = 1:50)
DGEGOI_exp.integrated <- FindClusters(DGEGOI_exp.integrated, resolution = 0.03)
DGEGOI_exp.integrated <- RunUMAP(DGEGOI_exp.integrated, reduction = "pca", dims = 1:50)
DGEGOI_exp.integrated[["whole_experiment_seurat_PCA_clusters"]] <- Idents(object = DGEGOI_exp.integrated)

# for(i in c("pca", "umap")) {
#   p1 <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Tissue")
#   p2 <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Replicate")
#   p3 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
#   p4 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
#   p5 <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Phase")
#   p6 <- DimPlot(DGEGOI_exp.integrated, reduction = i, label = TRUE)
#   gridit <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
#   ggsave(paste0("Seurat_clusters_", i, ".png"), plot = gridit, device = "png")
# }

# Run PHATE on intergated dataset and determine clusters using kmeans
# --------------------------------------------------------------------------

# Run phate, find number of clusters in embeddings using silhouette
phate.out <- phate(Matrix::t(GetAssayData(DGEGOI_exp.integrated, assay = "integrated", slot = "data")), ndim = 10) # Runs PHATE diffusion map
DGEGOI_exp.integrated[["phate"]] <- CreateDimReducObject(embeddings = phate.out$embedding, key = "PHATE_", assay = "integrated")
ProjectDim(DGEGOI_exp.integrated, reduction = "phate")
DGEGOI_exp.integrated <- FindNeighbors(DGEGOI_exp.integrated, reduction = "phate", dims = 1:10)
DGEGOI_exp.integrated <- FindClusters(DGEGOI_exp.integrated, resolution = 0.03)
DGEGOI_exp.integrated[["whole_experiment_PHATE_clusters"]] <- Idents(object = DGEGOI_exp.integrated)

phate_embed <- data.frame(Embeddings(DGEGOI_exp.integrated, reduction = "phate"))

for(i in c("pca", "umap")) {
  pLiver <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Liver", label = "Liver", order = "TRUE")
  pLN <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "LN", label = "LN", order = "TRUE")
  pLung <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Lung", label = "Lung", order = "TRUE")
  pPrimary <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Primary", label = "Primary", order = "TRUE")
  p1 <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Tissue")
  p2 <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
  p4 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
  p5 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("Mito_percent"), sort.cell = TRUE)
  p6 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("Lib_size"), sort.cell = TRUE)
  p7 <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Phase")
  p8 <- DimPlot(DGEGOI_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
  gridit_tissue <- plot_grid(pLiver, pLN, pLung, pPrimary, nrow = 2)
  ggsave(paste0("PHATE_clusters_", i, "_DGEGOI_QC.png"), plot = gridit, device = "png")
  ggsave(paste0("PHATE_clusters_", i, "_DGEGOI_tissues.png"), plot = gridit_tissue, device = "png")
}

# No idea why phate needs coordinates set
for(i in c("phate")) {
  pLiver <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Liver", label = "Liver", order = "TRUE")
  pLN <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "LN", label = "LN", order = "TRUE")
  pLung <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Lung", label = "Lung", order = "TRUE")
  pPrimary <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Primary", label = "Primary", order = "TRUE")
  p1 <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Tissue")
  p2 <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
  p3 <- p3 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
  p4 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
  p4 <- p4 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
  p5 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("Mito_percent"), sort.cell = TRUE)
  p5 <- p5 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
  p6 <- FeaturePlot(DGEGOI_exp.integrated, reduction = i, features = c("Lib_size"), sort.cell = TRUE)
  p6 <- p6 + xlim(c(min(phate_embed$PHATE_1), max(phate_embed$PHATE_1))) + ylim(c(min(phate_embed$PHATE_2), max(phate_embed$PHATE_2)))
  p7 <- DimPlot(DGEGOI_exp.integrated, reduction = i, group.by = "Phase")
  p8 <- DimPlot(DGEGOI_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, nrow = 4)
  gridit_tissue <- plot_grid(pLiver, pLN, pLung, pPrimary, nrow = 2)
  ggsave(paste0("PHATE_clusters_", i, "_DGEGOI_QC.png"), plot = gridit, device = "png")
  ggsave(paste0("PHATE_clusters_", i, "_DGEGOI_tissues.png"), plot = gridit_tissue, device = "png")
}

DGEGOI_exp.integrated@assays$RNA@meta.features <- merge(DGEGOI_exp.integrated@assays$RNA@meta.features, rowMetaData, by.x = 0, by.y = 0)

# Save seurat objects
# --------------------------------------------------------------------------

if(place == "local" & exists("phate.out")) {
  # saveRDS(phate.out, "Prefiltered_experiment_practice_phate.out.rds")
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(phate.out, "DGEGOI_experiment_all_phate.out.rds")
} else {
  print("Not overwritten")
}

if(place == "local" & exists("phate.out")) {
  # saveRDS(DGEGOI_exp.list, "Prefiltered_experiment_practice_seurat_list.rds")
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(DGEGOI_exp.list, "DGEGOI_experiment_all_seurat_list.rds")
} else {
  print("Not overwritten")
}

if(place == "local" & exists("phate.out")) {
  # saveRDS(DGEGOI_exp.anchors, "Prefiltered_experiment_practice_seurat_anchors.rds")
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(DGEGOI_exp.anchors, "DGEGOI_experiment_all_seurat_anchors.rds")
} else {
  print("Not overwritten")
}

if(place == "local" & exists("phate.out")) {
  saveRDS(DGEGOI_exp.integrated, "DGEGOI_experiment_practice_seurat_integrated.rds")
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(DGEGOI_exp.integrated, "DGEGOI_experiment_all_seurat_integrated.rds")
} else {
  print("Not overwritten")
}
