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
digest_stress <- list(c("FOS", "CXCL2", "ZFP36", "FOSB", "DUSP1", "ATF3", "CXCL8", "NR4A1", "CXCL3", "PPP1R15A", "JUNB", "EGR1", "HSPA1A", "HSPA1B", "SOCS3", "KLF6", "JUN", "IER2", "CXCL1", "NKFBIA", "HSPA6", "DNAJB1", "IER3", "CCNL1", "MTRNR2L2", "IER5", "ID1", "CEBPD", "KRT6A", "CYR61", "DEPP1", "CLDN4", "IRF1", "DUSP2", "BTG2", "PLAUR", "MAFF", "KLF4", "PHLDA2", "TNFAIP3"))
tissue_exp_seurat <- AddModuleScore(object = tissue_exp_seurat, features = digest_stress, name = 'digest_stress')
Idents(object = tissue_exp_seurat) <- "Tissue"

# Add score to idenitfy dying cells with low mito content (as per DOI 10.1186/s13059-019-1830-0)
dying <- list(c("HLA-A", "HLA-B", "HLA-C", "B2M"))
tissue_exp_seurat <- AddModuleScore(object = tissue_exp_seurat, features = dying, name = 'dying')

tissue_exp.list <- SplitObject(tissue_exp_seurat, split.by = "Sample") # split into individual samples

# Normalise transform counts within each experiment. Note: will not overwrite if already exists.
# --------------------------------------------------------------------------

# Perform SCT normalisation on each dataset individually
  for (i in 1:length(tissue_exp.list)) {
    tissue_exp.list[[i]] <- SCTransform(tissue_exp.list[[i]], verbose = TRUE, vars.to.regress = c("Replicate", "Lib_size", "S.Score", "G2M.Score", "digest_stress1"))
  }

# Integrate datasets based on highly correlated features
tissue_exp.features <- SelectIntegrationFeatures(object.list = tissue_exp.list, nfeatures = 5000)
tissue_exp.list <- PrepSCTIntegration(object.list = tissue_exp.list, verbose = TRUE, anchor.features = tissue_exp.features)
# tissue_exp.list <- lapply(X = tissue_exp.list, FUN = RunPCA, verbose = TRUE, features = tissue_exp.features) # Perform PCA on each object individually (needed for rpca)
reference_datasets <- which(names(filtered_exp.list) == "3")
tissue_exp.anchors <- FindIntegrationAnchors(object.list = tissue_exp.list, normalization.method = "SCT", anchor.features = tissue_exp.features, verbose = TRUE, reduction = "cca", reference = reference_datasets)
tissue_exp.integrated <- IntegrateData(anchorset = tissue_exp.anchors, normalization.method = "SCT", verbose = TRUE)

# Run PCA on intergated dataset and determine clusters using Seurat
# --------------------------------------------------------------------------

# Seurat integration pipeline
tissue_exp.integrated <- RunPCA(tissue_exp.integrated, dims = 1:50, assay = "integrated", ndims.print = 1:5, nfeatures.print = 5)
tissue_exp.integrated <- FindNeighbors(tissue_exp.integrated, reduction = "pca", dims = 1:50)
tissue_exp.integrated <- FindClusters(tissue_exp.integrated, resolution = 0.01)
tissue_exp.integrated <- RunUMAP(tissue_exp.integrated, reduction = "pca", dims = 1:50)
tissue_exp.integrated[["seurat_PCA_clusters"]] <- Idents(object = tissue_exp.integrated)

for(i in c("pca", "umap")) {
  p1 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Tissue")
  p2 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
  p4 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
  p5 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Phase")
  p6 <- DimPlot(tissue_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
  ggsave(paste0(tissue, "/Seurat_clusters_", tissue, "_", i, ".png", plot = gridit), device = "png")
}

# Run PHATE on intergated dataset and determine clusters using kmeans
# --------------------------------------------------------------------------

# Run phate, find number of clusters in embeddings using silhouette
phate.out <- phate(Matrix::t(GetAssayData(tissue_exp.integrated, assay = "integrated", slot = "data")), ndim = 10) # Runs PHATE diffusion map
tissue_exp.integrated[["phate"]] <- CreateDimReducObject(embeddings = phate.out$embedding, key = "PHATE_", assay = "integrated")
ProjectDim(tissue_exp.integrated, reduction = "phate")
tissue_exp.integrated <- FindNeighbors(tissue_exp.integrated, reduction = "phate", dims = 1:10)
tissue_exp.integrated <- FindClusters(tissue_exp.integrated, resolution = 0.01)
phate_embed <- data.frame(Embeddings(tissue_exp.integrated, reduction = "phate"))

for(i in c("pca", "umap", "phate")) {
  p1 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Tissue")
  p2 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Replicate")
  p3 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("digest_stress1"), sort.cell = TRUE)
  p4 <- FeaturePlot(tissue_exp.integrated, reduction = i, features = c("dying1"), sort.cell = TRUE)
  p5 <- DimPlot(tissue_exp.integrated, reduction = i, group.by = "Phase")
  p6 <- DimPlot(tissue_exp.integrated, reduction = i, label = TRUE)
  gridit <- plot_grid(p1, p2, p3, p4, p5, p6, nrow = 3)
  ggsave(paste0(tissue, "/PHATE_clusters_", tissue, "_", i, ".png", plot = gridit), device = "png")
}

# Save seurat objects
# --------------------------------------------------------------------------

if(place == "local" & exists("phate.out")) {
  # saveRDS(phate.out, paste0(tissue, "/Prefiltered_experiment_practice_phate.out_", tissue,".rds"))
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(phate.out, paste0(tissue, "/Prefiltered_experiment_all_phate.out_", tissue,".rds"))
} else {
  print("Not overwritten")
}

if(place == "local" & exists("tissue_exp.list")) {
  # saveRDS(tissue_exp.list, paste0(tissue, "/Prefiltered_experiment_practice_seurat_list_", tissue,".rds"))
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(tissue_exp.list, paste0(tissue, "/Prefiltered_experiment_all_seurat_list_", tissue,".rds"))
} else {
  print("Not overwritten")
}

if(place == "local" & exists("tissue_exp.anchors")) {
  # saveRDS(tissue_exp.anchors, paste0(tissue, "/Prefiltered_experiment_practice_seurat_anchors_", tissue,".rds"))
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(tissue_exp.anchors, paste0(tissue, "/Prefiltered_experiment_all_seurat_anchors_", tissue,".rds"))
} else {
  print("Not overwritten")
}

if(place == "local" & exists("tissue_exp.integrated")) {
  saveRDS(tissue_exp.integrated, paste0(tissue, "/Prefiltered_experiment_practice_seurat_integrated_", tissue,".rds"))
} else if(place == "wolfpack" & exists("phate.out")) {
  saveRDS(tissue_exp.integrated, paste0(tissue, "/Prefiltered_experiment_all_seurat_integrated_", tissue,".rds"))
} else {
  print("Not overwritten")
}
}
