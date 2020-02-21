#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells within each tissue
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCmets_chrcha/project_results/prefiltered/by_organ") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCmets_chrcha/project_results/prefiltered/by_organ")
  place <- "wolfpack"
}

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('pheatmap')
library('phateR')
library('ggplot2')
library('readr')
library('Rmagic')
library('batchelor')
library('Matrix')

set.seed(100)
# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("/Users/mac/cloudstor/sarah_projects/SCmets_chrcha/project_results/prefiltered/practice_all_data/Prefiltered_experiment_Practice_merge_cluster_wCC.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCmets_chrcha/project_results/prefiltered/all_data/Prefiltered_experiment_All_merge_cluster_wCC.rds") # uses whole dataset if wolfpack
}

# Subset filtered experiment to each tissue
for(i in unique(filtered_exp$Tissue)) {
  tissue_exp <- filtered_exp[,filtered_exp$Tissue == i]
  
  # Normalised to adjust for differences between samples
  tissue_exp <- multiBatchNorm(tissue_exp, batch = tissue_exp$Sample)
  
  # Select HVG based on combined variance across all samples
  tissue_exp.dec <- modelGeneVarByPoisson(tissue_exp, block = tissue_exp$Sample)
  HVG <- tissue_exp.dec$bio > 0
  rowData(tissue_exp)$bio_var <- tissue_exp.dec$bio
  rowData(tissue_exp)$bio_var_FDR <- tissue_exp.dec$FDR
  rowData(tissue_exp)$HVG <- tissue_exp.dec$bio > 0
  rowData(tissue_exp)$HVG_sig <- tissue_exp.dec$bio > 0 & tissue_exp.dec$FDR < 0.05
  
  pdf(paste0(i, "/Mean_variance_model_fit_poisson_", i, ".pdf"))
  plot(tissue_exp.dec$mean, tissue_exp.dec$total, pch=16, xlab="Mean of log-expression", ylab="Variance of log-expression")
  # curve(metadata(tissue_exp.dec)$trend(x), col="dodgerblue", add=TRUE)
  dev.off()
  
  # Run multibatch PCA and clustering without correction
  sample_details <- data.frame(unique(tissue_exp$Sample)) %>%
    separate(1, c("Tissue", "Replicate"), "_", remove = FALSE)
  
  multiPCA <- multiBatchPCA(tissue_exp,
                            batch = tissue_exp$Sample,
                            subset.row=HVG,
                            get.all.genes = TRUE,
                            preserve.single = TRUE,
                            get.variance = TRUE,
                            BSPARAM=BiocSingular::IrlbaParam(deferred=TRUE)) # PCs used later for reducedMNN
  
  PCA50 <- as.matrix(data.frame((multiPCA@listData)))
  colnames(PCA50) <- paste0("PC", 1:50)
  attr(PCA50, 'percentVar') <- multiPCA@metadata$var.explained
  reducedDim(tissue_exp, "PCA_tissue") <- PCA50
  
  # tissue_exp <- denoisePCA(tissue_exp, technical=tissue_exp.dec, subset.row=HVG) # Keep PCs which associated with significant biological variation
  
  snn.gr <- buildSNNGraph(tissue_exp, use.dimred="PCA_tissue", k=20)
  clusters <- igraph::cluster_walktrap(snn.gr)$membership
  uncorrected_tab <- table(Cluster=clusters, Batch=tissue_exp$Sample)
  write.csv(uncorrected_tab, paste0(i, "/Uncorrected_batch_cell_cluster_membership_", i, ".csv"))
  uncorrected_tab
  
  # Visualise uncorrected clusters using PCA, UMAP and PHATE
  tissue_exp$uncorrected_cluster_tissue <- factor(clusters)
  plotReducedDim(tissue_exp, dimred="PCA_tissue", colour_by = "Mito_percent", text_by = "uncorrected_cluster_tissue") +
    scale_fill_viridis_c(option = "C") +
    ggsave(paste0(i, "/PCA_uncorrected_with_clusters_Mito_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="PCA_tissue", colour_by = "Replicate", text_by = "uncorrected_cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/PCA_uncorrected_with_clusters_replicate_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="PCA_tissue", colour_by = "uncorrected_cluster_tissue", text_by = "uncorrected_cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/PCA_uncorrected_with_clusters_", i, ".pdf"))
  
  tissue_exp <- runUMAP(tissue_exp, dimred="PCA_tissue", name = "UMAP_tissue")
  plotReducedDim(tissue_exp, dimred="UMAP_tissue", colour_by = "Mito_percent", text_by = "uncorrected_cluster_tissue") +
    scale_fill_viridis_c(option = "C") +
    ggsave(paste0(i, "/UMAP_uncorrected_with_clusters_Mito_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="UMAP_tissue", colour_by = "Replicate", text_by = "uncorrected_cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/UMAP_uncorrected_with_clusters_replicate_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="UMAP_tissue", colour_by = "uncorrected_cluster_tissue", text_by = "uncorrected_cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/UMAP_uncorrected_with_clusters_", i, ".pdf"))
  
  phate.tree <- phate(Matrix::t(assay(tissue_exp, "logcounts"))) # Runs PHATE diffusion map
  reducedDim(tissue_exp, "PHATE_tissue") <- phate.tree$embedding
  plotReducedDim(tissue_exp, dimred="PHATE_tissue", colour_by = "Mito_percent", text_by = "uncorrected_cluster_tissue") +
    scale_fill_viridis_c(option = "C") +
    ggsave(paste0(i, "/PHATE_uncorrected_with_clusters_Mito_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="PHATE_tissue", colour_by = "Replicate", text_by = "uncorrected_cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/PHATE_uncorrected_with_clusters_replicate_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="PHATE_tissue", colour_by = "uncorrected_cluster_tissue", text_by = "uncorrected_cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/PHATE_uncorrected_with_clusters_", i, ".pdf"))
  
  # Run clustering with correction for batch
  merge_order <- list(LN = list(c("LN_3", "LN_4", "LN_2", "LN_1")),
                      Primary = list(c("Primary_4", "Primary_3", "Primary_1", "Primary_2")),
                      Liver = list(c("Liver_3", "Liver_2", "Liver_4", "Liver_1")),
                      Lung = list(c("Lung_3", "Lung_4", "Lung_1", "Lung_2")))
  
  fastMNN.sce <- fastMNN(tissue_exp,
                         subset.row=HVG,
                         cos.norm = FALSE,
                         k=20,
                         correct.all = TRUE,
                         batch = tissue_exp$Sample,
                         merge.order = merge_order[i])
  
  snn.gr <- buildSNNGraph(fastMNN.sce, use.dimred = "corrected", k=20)
  clusters <- igraph::cluster_walktrap(snn.gr)$membership
  corrected_tab <- table(Cluster=clusters, Batch=fastMNN.sce$batch)
  write.csv(corrected_tab, paste0(i, "/Corrected_batch_cell_cluster_membership_", i, ".csv"))
  corrected_tab
  
  colSums(metadata(fastMNN.sce)$merge.info$lost.var)
  
  # Visualise corrected clusters using PCA, UMAP and PHATE
  tissue_exp$cluster_tissue <- factor(clusters)
  reducedDim(tissue_exp, "corrected_fastMNN_tissue") <- reducedDim(fastMNN.sce, "corrected")
  assay(tissue_exp, "reconstructed_fastMNN_tissue") <- assay(fastMNN.sce, "reconstructed")
  
  plotReducedDim(tissue_exp, dimred="corrected_fastMNN_tissue", colour_by = "Mito_percent", text_by = "cluster_tissue") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "/Fastmnn_corrected_with_clusters_Mito_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="corrected_fastMNN_tissue", colour_by = "Genes_detected", text_by = "cluster_tissue") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "/Fastmnn_corrected_with_clusters_NGenes_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="corrected_fastMNN_tissue", colour_by = "Lib_size", text_by = "cluster_tissue") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "/Fastmnn_corrected_with_clusters_Lib_size_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="corrected_fastMNN_tissue", colour_by = "CC.Phase", text_by = "cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/Fastmnn_corrected_with_clusters_cell.cycle.phase_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="corrected_fastMNN_tissue", colour_by = "Replicate", text_by = "cluster_tissue") +
    scale_fill_viridis_d(option = "C") +
    ggsave(paste0(i, "/Fastmnn_corrected_with_clusters_replicate_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="corrected_fastMNN_tissue", colour_by = "cluster_tissue", text_by = "cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/Fastmnn_corrected_with_clusters_", i, ".pdf"))
  
  tissue_exp <- runUMAP(tissue_exp, dimred="corrected_fastMNN_tissue", name = "UMAP_fastMNN_tissue")
  plotReducedDim(tissue_exp, dimred="UMAP_fastMNN_tissue", colour_by = "Mito_percent", text_by = "cluster_tissue") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "/UMAP_corrected_with_clusters_Mito_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="UMAP_fastMNN_tissue", colour_by = "Genes_detected", text_by = "cluster_tissue") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "/UMAP_corrected_with_clusters_NGenes_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="UMAP_fastMNN_tissue", colour_by = "Lib_size", text_by = "cluster_tissue") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "/UMAP_corrected_with_clusters_Lib_size_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="UMAP_fastMNN_tissue", colour_by = "CC.Phase", text_by = "cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/UMAP_corrected_with_clusters_cell.cycle.phase_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="UMAP_fastMNN_tissue", colour_by = "Replicate", text_by = "cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/UMAP_corrected_with_clusters_replicate_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="UMAP_fastMNN_tissue", colour_by = "cluster_tissue", text_by = "cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/UMAP_corrected_with_clusters_", i, ".pdf"))
  
  phate.out <- phate(Matrix::t(assay(tissue_exp, "reconstructed_fastMNN_tissue"))) # Runs PHATE diffusion map
  reducedDim(tissue_exp, "PHATE_fastMNN_tissue") <- phate.out$embedding
  plotReducedDim(tissue_exp, dimred="PHATE_fastMNN_tissue", colour_by = "Mito_percent", text_by = "cluster_tissue") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "/PHATE_corrected_with_clusters_Mito_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="PHATE_fastMNN_tissue", colour_by = "Genes_detected", text_by = "cluster_tissue") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "/PHATE_corrected_with_clusters_NGenes_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="PHATE_fastMNN_tissue", colour_by = "Lib_size", text_by = "cluster_tissue") +
    scale_fill_viridis_c(option = "B") +
    ggsave(paste0(i, "/PHATE_corrected_with_clusters_Lib_size_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="PHATE_fastMNN_tissue", colour_by = "CC.Phase", text_by = "cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/PHATE_corrected_with_clusters_cell.cycle.phase_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="PHATE_fastMNN_tissue", colour_by = "Replicate", text_by = "cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/PHATE_corrected_with_clusters_replicate_", i, ".pdf"))
  plotReducedDim(tissue_exp, dimred="PHATE_fastMNN_tissue", colour_by = "cluster_tissue", text_by = "cluster_tissue") +
    scale_fill_viridis_d(option = "D") +
    ggsave(paste0(i, "/PHATE_corrected_with_clusters_", i, ".pdf"))
  
  # Save total filtered dataset
  if(place == "local") {
    saveRDS(tissue_exp, paste0(i, "/Prefiltered_experiment_Practice_merge_cluster_wCC_", i, ".rds"))
  } else {
    saveRDS(tissue_exp, paste0(i, "/Prefiltered_experiment_All_merge_cluster_wCC_", i, ".rds"))
  }
  
}
