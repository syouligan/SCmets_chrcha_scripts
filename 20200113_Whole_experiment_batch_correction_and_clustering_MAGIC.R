#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCmets_chrcha/project_results/prefiltered/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCmets_chrcha/project_results/prefiltered/all_data")
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
  filtered_exp <- readRDS("Prefiltered_experiment_Practice.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_experiment_All.rds") # uses whole dataset if wolfpack
}

# Normalised to adjust for differences between samples
filtered_exp <- multiBatchNorm(filtered_exp, batch = filtered_exp$Sample)
assay(filtered_exp, "logcounts") <- Matrix(t(as.matrix(magic(t(as.matrix(assay(filtered_exp, "logcounts"))), genes = "all_genes", t = 'auto', n.jobs = 4))), sparse = TRUE)

# Select HVG based on combined variance across all samples
filtered_exp.dec <- modelGeneVarByPoisson(filtered_exp, block = filtered_exp$Sample)
HVG <- filtered_exp.dec$bio > 0
rowData(filtered_exp)$bio_var <- filtered_exp.dec$bio
rowData(filtered_exp)$bio_var_FDR <- filtered_exp.dec$FDR
rowData(filtered_exp)$HVG <- filtered_exp.dec$bio > 0
rowData(filtered_exp)$HVG_sig <- filtered_exp.dec$bio > 0 & filtered_exp.dec$FDR < 0.05

pdf("Mean_variance_model_fit_poisson_MAGIC.pdf")
plot(filtered_exp.dec$mean, filtered_exp.dec$total, pch=16, xlab="Mean of log-expression", ylab="Variance of log-expression")
# curve(metadata(filtered_exp.dec)$trend(x), col="dodgerblue", add=TRUE)
dev.off()

# Run multibatch PCA and clustering without correction
sample_details <- data.frame(unique(filtered_exp$Sample)) %>%
  separate(1, c("Tissue", "Met", "Replicate"), "_", remove = FALSE)

met_weights <- ifelse(sample_details$Met == "NA", 1, 2) # Weight PCA for the number of mets from each sample
names(met_weights) <- unique(filtered_exp$Sample)
met_weights

multiPCA <- multiBatchPCA(filtered_exp,
                          batch = filtered_exp$Sample,
                          subset.row=HVG,
                          get.all.genes = TRUE,
                          preserve.single = TRUE,
                          weights = 1/met_weights,
                          get.variance = TRUE,
                          BSPARAM=BiocSingular::IrlbaParam(deferred=TRUE)) # PCs used later for reducedMNN

PCA50 <- as.matrix(data.frame((multiPCA@listData)))
colnames(PCA50) <- paste0("PC", 1:50)
attr(PCA50, 'percentVar') <- multiPCA@metadata$var.explained
reducedDim(filtered_exp, "PCA") <- PCA50

# filtered_exp <- denoisePCA(filtered_exp, technical=filtered_exp.dec, subset.row=HVG) # Keep PCs which associated with significant biological variation

snn.gr <- buildSNNGraph(filtered_exp, use.dimred="PCA")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
uncorrected_tab <- table(Cluster=clusters, Batch=filtered_exp$Sample)
write.csv(uncorrected_tab, "Uncorrected_batch_cell_cluster_membership_MAGIC.csv")
uncorrected_tab

# Visualise uncorrected clusters using PCA, UMAP and PHATE
filtered_exp$uncorrected_cluster <- clusters
plotReducedDim(filtered_exp, dimred="PCA", colour_by = "Tissue", text_by = "uncorrected_cluster") +
  ggsave("PCA_uncorrected_with_clusters_tissue_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="PCA", colour_by = "Replicate", text_by = "uncorrected_cluster") +
  ggsave("PCA_uncorrected_with_clusters_replicate_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="PCA", colour_by = "uncorrected_cluster", text_by = "uncorrected_cluster") +
  ggsave("PCA_uncorrected_with_clusters_MAGIC.pdf")

filtered_exp <- runUMAP(filtered_exp, dimred="PCA")
plotReducedDim(filtered_exp, dimred="UMAP", colour_by = "Tissue", text_by = "cluster") +
  ggsave("UMAP_uncorrected_with_clusters_tissue_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="UMAP", colour_by = "Replicate", text_by = "cluster") +
  ggsave("UMAP_uncorrected_with_clusters_replicate_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="UMAP", colour_by = "cluster", text_by = "cluster") +
  ggsave("UMAP_uncorrected_with_clusters_MAGIC.pdf")

phate.tree <- phate(t(as.matrix(assay(filtered_exp, "logcounts")))) # Runs PHATE diffusion map
reducedDim(filtered_exp, "PHATE") <- phate.tree$embedding
plotReducedDim(filtered_exp, dimred="PHATE", colour_by = "Tissue", text_by = "uncorrected_cluster") +
  ggsave("PHATE_uncorrected_with_clusters_tissue_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="PHATE", colour_by = "Replicate", text_by = "uncorrected_cluster") +
  ggsave("PHATE_uncorrected_with_clusters_replicate_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="PHATE", colour_by = "uncorrected_cluster", text_by = "uncorrected_cluster") +
  ggsave("PHATE_uncorrected_with_clusters_MAGIC.pdf")

# Run clustering with correction for batch
merge_order <- list(list(unique(filtered_exp$Sample)[sample_details$Tissue == "Primary"]),
                    list(unique(filtered_exp$Sample)[sample_details$Tissue == "Liver"]),
                    list(unique(filtered_exp$Sample)[sample_details$Tissue == "Lung"]),
                    list(unique(filtered_exp$Sample)[sample_details$Tissue == "LN"]))

fastMNN.sce <- fastMNN(filtered_exp,
                       subset.row=HVG,
                       cos.norm = FALSE,
                       correct.all = TRUE,
                       batch = filtered_exp$Sample,
                       merge.order = merge_order,
                       weights = 1/met_weights)

snn.gr <- buildSNNGraph(fastMNN.sce, use.dimred="corrected")
clusters <- igraph::cluster_walktrap(snn.gr)$membership
corrected_tab <- table(Cluster=clusters, Batch=fastMNN.sce$batch)
write.csv(corrected_tab, "Corrected_batch_cell_cluster_membership_MAGIC.csv")
corrected_tab

colSums(metadata(fastMNN.sce)$merge.info$lost.var)

# Group samples based on cell composition with and without batch correction
sample_details <- data.frame(colnames(corrected_tab)) %>%
  separate(1, c("Tissue", "Met", "Replicate"), "_", remove = FALSE)

pc <- data.frame(prcomp(as.matrix(uncorrected_tab/colSums(uncorrected_tab)), scale. = TRUE)$rotation)
rownames(pc) <- colnames(uncorrected_tab)
pc$Tissue <- sample_details$Tissue
pc$Replicate <- sample_details$Replicate

ggplot(pc, aes(x=PC1, y=PC2, label = rownames(pc))) +
  stat_ellipse(aes(fill = Tissue), alpha = 0.3, geom = "polygon", type = "norm", level = 0.5, color = "black") +
  geom_point(aes(color = Tissue)) +
  geom_text(aes(color = Tissue)) +
  theme_classic() +
  ggsave("Sample_clustering_before_fastmnn_cell_composition_Tissue_MAGIC.pdf", useDingbats = FALSE)

ggplot(pc, aes(x=PC1, y=PC2, label = rownames(pc))) +
  stat_ellipse(aes(fill = Replicate), alpha = 0.3, geom = "polygon", type = "norm", level = 0.5, color = "black") +
  geom_point(aes(color = Replicate)) +
  geom_text(aes(color = Replicate)) +
  theme_classic() +
  ggsave("Sample_clustering_before_fastmnn_cell_composition_Replicate_MAGIC.pdf", useDingbats = FALSE)

pc <- data.frame(prcomp(as.matrix(corrected_tab/colSums(corrected_tab)), scale. = TRUE)$rotation)
rownames(pc) <- colnames(corrected_tab)
pc$Tissue <- sample_details$Tissue
pc$Replicate <- sample_details$Replicate

ggplot(pc, aes(x=PC1, y=PC2, label = rownames(pc))) +
  stat_ellipse(aes(fill = Tissue), alpha = 0.3, geom = "polygon", type = "norm", level = 0.5, color = "black") +
  geom_point(aes(color = Tissue)) +
  geom_text(aes(color = Tissue)) +
  theme_classic() +
  ggsave("Sample_clustering_after_fastmnn_cell_composition_Tissue_MAGIC.pdf", useDingbats = FALSE)

ggplot(pc, aes(x=PC1, y=PC2, label = rownames(pc))) +
  stat_ellipse(aes(fill = Replicate), alpha = 0.3, geom = "polygon", type = "norm", level = 0.5, color = "black") +
  geom_point(aes(color = Replicate)) +
  geom_text(aes(color = Replicate)) +
  theme_classic() +
  ggsave("Sample_clustering_after_fastmnn_cell_composition_Replicate_MAGIC.pdf", useDingbats = FALSE)

# Visualise corrected clusters using PCA, UMAP and PHATE
filtered_exp$cluster <- clusters
reducedDim(filtered_exp, "corrected_fastMNN") <- reducedDim(fastMNN.sce, "corrected")
assay(filtered_exp, "reconstructed_fastMNN") <- assay(fastMNN.sce, "reconstructed")

plotReducedDim(filtered_exp, dimred="corrected_fastMNN", colour_by = "Tissue", text_by = "cluster") +
  ggsave("Fastmnn_corrected_with_clusters_tissue_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="corrected_fastMNN", colour_by = "Replicate", text_by = "cluster") +
  ggsave("Fastmnn_corrected_with_clusters_replicate_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="corrected_fastMNN", colour_by = "cluster", text_by = "cluster") +
  ggsave("Fastmnn_corrected_with_clusters_MAGIC.pdf")

filtered_exp <- runUMAP(filtered_exp, dimred="corrected_fastMNN")
plotReducedDim(filtered_exp, dimred="corrected_fastMNN", colour_by = "Tissue", text_by = "cluster") +
  ggsave("UMAP_corrected_with_clusters_tissue_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="corrected_fastMNN", colour_by = "Replicate", text_by = "cluster") +
  ggsave("UMAP_corrected_with_clusters_replicate_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="corrected_fastMNN", colour_by = "cluster", text_by = "cluster") +
  ggsave("UMAP_corrected_with_clusters_MAGIC.pdf")

phate.tree <- phate(t(as.matrix(assay(filtered_exp, "reconstructed_fastMNN")))) # Runs PHATE diffusion map
reducedDim(filtered_exp, "PHATE_fastMNN") <- phate.tree$embedding
plotReducedDim(filtered_exp, dimred="PHATE_fastMNN", colour_by = "Tissue", text_by = "cluster") +
  ggsave("PHATE_corrected_with_clusters_tissue_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="PHATE_fastMNN", colour_by = "Replicate", text_by = "cluster") +
  ggsave("PHATE_corrected_with_clusters_replicate_MAGIC.pdf")
plotReducedDim(filtered_exp, dimred="PHATE_fastMNN", colour_by = "cluster", text_by = "cluster") +
  ggsave("PHATE_corrected_with_clusters_MAGIC.pdf")

# Save total filtered dataset
saveRDS(filtered_exp, "Prefiltered_experiment_All_merge_cluster_MAGIC.rds")

