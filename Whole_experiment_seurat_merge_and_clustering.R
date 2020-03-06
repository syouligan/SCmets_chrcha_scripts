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

set.seed(100)
options(future.globals.maxSize = 200000*1024^2)

# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/practice_all_data/Prefiltered_experiment_practice_merge_cluster_wCC.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/Prefiltered_experiment_all_merge_cluster_wCC.rds") # uses whole dataset if wolfpack
}

# Set up prefiltered object, normalise transform and scale counts
filtered_exp_seurat <- as.Seurat(filtered_exp, counts = "counts", data = NULL) # convert to Seurat
filtered_exp.list <- SplitObject(filtered_exp_seurat, split.by = "Sample") # split into individual samples

# Normalise transform counts within each experiment
for (i in 1:length(filtered_exp.list)) {
  filtered_exp.list[[i]] <- SCTransform(filtered_exp.list[[i]], verbose = TRUE, vars.to.regress = "nCount_RNA")
}

# Integrate datasets based on highly correlated features
filtered_exp.features <- SelectIntegrationFeatures(object.list = filtered_exp.list, nfeatures = 3000)
filtered_exp.list <- PrepSCTIntegration(object.list = filtered_exp.list, verbose = TRUE, anchor.features = filtered_exp.features)
# filtered_exp.list <- lapply(X = filtered_exp.list, FUN = RunPCA, verbose = TRUE, features = filtered_exp.features) # Perform PCA on each object individually (needed for rpca)
filtered_exp.anchors <- FindIntegrationAnchors(object.list = filtered_exp.list, normalization.method = "SCT", anchor.features = filtered_exp.features, verbose = TRUE, reduction = "cca", reference = which(names(filtered_exp.list) == c("LN_1", "Liver_1", "Lung_1", "Primary_1")))
filtered_exp.integrated <- IntegrateData(anchorset = filtered_exp.anchors, normalization.method = "SCT", verbose = TRUE)

# Run PCA on intergated dataset and determine clusters
filtered_exp.integrated <- RunPCA(filtered_exp.integrated, npcs = 50, ndims.print = 1:5, nfeatures.print = 5)
DimHeatmap(filtered_exp.integrated, dims = c(1:10), cells = 500, balanced = TRUE) +
  ggsave("PCs_seurat_with_clusters.pdf")
filtered_exp.integrated <- FindNeighbors(filtered_exp.integrated, reduction = "pca", dims = 1:50)
filtered_exp.integrated <- FindClusters(filtered_exp.integrated)

filtered_exp.integrated_sce <- as.SingleCellExperiment(filtered_exp.integrated)
reducedDim(filtered_exp, "PCA_seurat") <- reducedDim(filtered_exp.integrated_sce, "PCA")
filtered_exp$seurat_clusters <- filtered_exp.integrated_sce$ident

phate.out <- phate(Matrix::t(assay(filtered_exp.integrated_sce, "logcounts"))) # Runs PHATE diffusion map
reducedDim(filtered_exp, "PHATE_seurat") <- phate.out$embedding
filtered_exp <- runUMAP(filtered_exp, dimred = "PCA_seurat", name = "UMAP_seurat")

for(i in c("PCA_seurat", "UMAP_seurat", "PHATE_seurat")){
  plotReducedDim(filtered_exp, dimred=i, colour_by = "Mito_percent", text_by = "seurat_clusters") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "_seurat_with_clusters_Mito.pdf"))
  plotReducedDim(filtered_exp, dimred=i, colour_by = "Genes_detected", text_by = "seurat_clusters") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "_seurat_with_clusters_NGenes.pdf"))
  plotReducedDim(filtered_exp, dimred=i, colour_by = "Lib_size", text_by = "seurat_clusters") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "_seurat_with_clusters_Lib_size.pdf"))
  plotReducedDim(filtered_exp, dimred=i, colour_by = "Tissue", text_by = "seurat_clusters") +
    scale_fill_viridis_d(option = "B") +
    ggsave(paste0(i, "_seurat_with_clusters_tissue.pdf"))
  plotReducedDim(filtered_exp, dimred=i, colour_by = "Replicate", text_by = "seurat_clusters") +
    scale_fill_viridis_d(option = "B") +
    ggsave(paste0(i, "_seurat_with_clusters_replicate.pdf"))
  plotReducedDim(filtered_exp, dimred=i, colour_by = "CC.Phase", text_by = "seurat_clusters") +
    scale_fill_viridis_d(option = "B") +
    ggsave(paste0(i, "_seurat_with_clusters_cell.cycle.pdf"))
  plotReducedDim(filtered_exp, dimred=i, colour_by = "cluster", text_by = "cluster") +
    scale_fill_viridis_d(option = "B") +
    ggsave(paste0(i, "_seurat_with_scran_clusters.pdf"))
  plotReducedDim(filtered_exp, dimred=i, colour_by = "seurat_clusters", text_by = "seurat_clusters") +
    scale_fill_viridis_d(option = "B") +
    ggsave(paste0(i, "_seurat_with_seurat_clusters.pdf"))
}

# Make dataframe used for plotting
# --------------------------------------------------------------------------

plotting <- data.frame(reducedDim(filtered_exp, "PHATE_seurat"))
plotting$Tissue <- filtered_exp$Tissue
plotting$cluster <- filtered_exp$seurat_clusters

# Quantitative expression
plotting$OCT4_expression <- assay(filtered_exp, "logcounts")["POU5F2", ]
plotting$SOX2_expression <- assay(filtered_exp, "logcounts")["SOX2", ]

# Scaled expression
plotting$OCT4_zscore <- scale(assay(filtered_exp, "logcounts")["POU5F2", ])
plotting$SOX2_zscore <- scale(assay(filtered_exp, "logcounts")["SOX2", ])
plotting$OCT4ySOX2_zscore <- plotting$OCT4_expression + plotting$SOX2_expression

# Qualitative expression
plotting$OCT4ySOX2_expressed <- plotting$OCT4_expression > 0 & plotting$SOX2_expression > 0
plotting$OCT4_expressed <- plotting$OCT4_expression > 0
plotting$SOX2_expressed <- plotting$SOX2_expression > 0

# Cluster labelling
text_out <- retrieveCellInfo(filtered_exp, "seurat_clusters", search="colData")
text_out$val <- as.factor(text_out$val)
by_text_x <- vapply(split(plotting$PHATE1, text_out$val), median, FUN.VALUE=0)
by_text_y <- vapply(split(plotting$PHATE2, text_out$val), median, FUN.VALUE=0)

# Plot cells within each cluster
# --------------------------------------------------------------------------

plotting_t <- plotting[plotting$OCT4ySOX2_expressed == TRUE, ]
ggplot() +
  geom_point(data = plotting, aes(x = PHATE1, y = PHATE2), color = "grey80") +
  geom_point(data = plotting_t, aes(x = PHATE1, y = PHATE2, color = OCT4ySOX2_zscore)) +
  scale_color_viridis_c(option = "C") +
  theme(aspect.ratio = 1) +
  theme_classic() +
  annotate("text", x = by_text_x, y = by_text_y, label = names(by_text_x)) +
  ggsave("select_genes/All_cells_OCT4ySOX2_scaled_all.pdf", useDingbats = FALSE)

ggplot() +
  geom_point(data = plotting, aes(x = PHATE1, y = PHATE2, color = OCT4ySOX2_expressed)) +
  scale_color_viridis_d(option = "D") +
  theme(aspect.ratio = 1) +
  theme_classic() +
  annotate("text", x = by_text_x, y = by_text_y, label = names(by_text_x)) +
  ggsave("select_genes/All_cells_OCT4_expressed_all.pdf", useDingbats = FALSE)

OCT4ySOX2_by_cluster <- table(cluster=plotting$cluster, OCT4ySOX2=plotting$OCT4ySOX2_expressed)
round(OCT4ySOX2_by_cluster/rowSums(OCT4ySOX2_by_cluster)*100)
OCT4ySOX2_by_cluster

# Save total filtered dataset
if(place == "local") {
  saveRDS(filtered_exp, "Prefiltered_experiment_practice_seurat_merge_cluster_sce.rds")
} else {
  saveRDS(filtered_exp, "Prefiltered_experiment_all_seurat_merge_cluster_sce.rds")
}

# Save intermediate seurat objects
if(place == "local") {
  saveRDS(filtered_exp.list, "Prefiltered_experiment_practice_seurat_list.rds")
} else {
  saveRDS(filtered_exp.list, "Prefiltered_experiment_all_seurat_list.rds")
}

if(place == "local") {
  saveRDS(filtered_exp.anchors, "Prefiltered_experiment_practice_seurat_anchors.rds")
} else {
  saveRDS(filtered_exp.anchors, "Prefiltered_experiment_all_seurat_anchors.rds")
}

if(place == "local") {
  saveRDS(filtered_exp.integrated_sce, "Prefiltered_experiment_practice_seurat_integrated.rds")
} else {
  saveRDS(filtered_exp.integrated_sce, "Prefiltered_experiment_all_seurat_integrated.rds")
}
