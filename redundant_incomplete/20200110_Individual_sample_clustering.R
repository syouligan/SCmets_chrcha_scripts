#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise and cluster individual samples
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  location <- "/Users/mac/cloudstor/"
} else {
  location <- "/share/ScratchGeneral/scoyou/"
}
setwd(paste0(location, "sarah_projects/SCmets_chrcha/project_results/prefiltered/individual"))

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('pheatmap')
library('destiny')
library('phateR')
library('ggplot2')
library('readr')
library('Rmagic')

# Perform on each sample
set.seed(100)
samples <- c("Liver_NA_1", "LN_NA_1", "Lung_NA_1", "Liver_NA_2", "Lung_NA_2", "LN_NA_2", "Liver_A_3", "Liver_B_3", "LN_A_3", "LN_B_3", "Lung_A_3", "Lung_B_3", "Liver_NA_4", "Lung_NA_4", "LN_NA_4", "Primary_NA_1", "Primary_NA_2", "Primary_NA_3", "Primary_NA_4")

for(i in samples) {
# i <- "Liver_A_3" # for development

# Load prefiltered SingleCellExperiment
sample_exp <- readRDS(paste0(i,"/Prefiltered_experiment_", i, ".rds"))

# Remove genes without counts in at least 3 cells for this sample
GT_three <- rowSums(counts(sample_exp) > 0) > 3
sample_exp <- sample_exp[which(GT_three), ]

# Normalise counts (by deconvolution) and identify variable features
clust.sample_exp <- quickCluster(sample_exp) 
sample_exp <- computeSumFactors(sample_exp, cluster=clust.sample_exp, min.mean=0.1)
sample_exp <- logNormCounts(sample_exp)
dec.sample_exp <- modelGeneVar(sample_exp)
rowData(sample_exp)$bio_var <- dec.sample_exp$bio
rowData(sample_exp)$bio_var_FDR <- dec.sample_exp$FDR
rowData(sample_exp)$HVG <- dec.sample_exp$bio > 0
rowData(sample_exp)$HVG_sig <- dec.sample_exp$bio > 0 & dec.sample_exp$FDR < 0.05

fit.sample_exp <- metadata(dec.sample_exp)
pdf(paste0(i, "/Mean_variance_model_fit_", i, ".pdf"))
plot(fit.sample_exp$mean, fit.sample_exp$var, xlab="Mean of log-expression", ylab="Variance of log-expression")
curve(fit.sample_exp$trend(x), col="dodgerblue", add=TRUE, lwd=2)
dev.off()

# Dimensionality reduction and visualisation
HVGs <- getTopHVGs(dec.sample_exp)
sample_exp <- runPCA(sample_exp, subset_row=HVGs)
sample_exp <- denoisePCA(sample_exp, technical=dec.sample_exp, subset.row=HVGs) # Keep PCs which associated with significant biological variation
sample_exp <- runUMAP(sample_exp, dimred="PCA")
sample_exp <- runDiffusionMap(sample_exp, dimred="PCA")
phate.tree <- phate(t(as.matrix(assay(sample_exp, "logcounts")))) # Runs PHATE diffusion map
reducedDim(sample_exp, "PHATE") <- phate.tree$embedding

# Cluster definition and check stability
myClusterFUN <- function(x) {
  g <- buildSNNGraph(x, use.dimred="PCA", k=15)
  igraph::cluster_walktrap(g)$membership
} # Function used for clustering to be bootstrapped. scran default.

cluster <- myClusterFUN(sample_exp)
sample_exp$cluster <- factor(cluster)
table(sample_exp$cluster)
plotReducedDim(sample_exp, dimred="PCA", colour_by = "cluster") +
  ggsave(paste0(i, "/PCA_with_clusters_", i, ".pdf"))
plotReducedDim(sample_exp, dimred="UMAP", colour_by="cluster") +
  ggsave(paste0(i, "/UMAP_with_clusters_", i, ".pdf"))
plotReducedDim(sample_exp, dimred="DiffusionMap", colour_by = "cluster") +
  ggsave(paste0(i, "/DiffusionMap_with_clusters_", i, ".pdf"))
plotReducedDim(sample_exp, dimred="PHATE", colour_by = "cluster") +
  ggsave(paste0(i, "/PHATE_with_clusters_", i, ".pdf"))

coassign <- bootstrapCluster(sample_exp, FUN=myClusterFUN, clusters=cluster) # Bootstraps clustering and evaluates how frequently genes are assigned to each cluster. Similar clusters have high coassignment.
pdf(paste0(i, "/Cluster_coassignment_stability_", i, ".pdf"))
pheatmap(coassign, cluster_row=FALSE, cluster_col=FALSE, color=rev(viridis::magma(100)))
dev.off()

g <- buildSNNGraph(sample_exp, use.dimred="PCA")
ratio <- clusterModularity(g, cluster, as.ratio=TRUE) # Checks the ratio of observed relative to expected edge weights between cells in each cluster.
pdf(paste0(i, "/Cluster_modularity_heatmap_", i, ".pdf"))
pheatmap(log2(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE, color=rev(viridis::magma(100)))
dev.off()

cluster.gr <- igraph::graph_from_adjacency_matrix(log2(ratio+1), mode="upper", weighted=TRUE, diag=FALSE)
pdf(paste0(i, "/Cluster_modularity_node-edges_", i, ".pdf"))
plot(cluster.gr, edge.width=igraph::E(cluster.gr)$weight*5, layout=igraph::layout_with_lgl)
dev.off()

# Assign cell cycle phases
# hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
# cyc.values <- cyclone(sample_exp, hs.pairs, gene.names = rowData(sample_exp)$Ensembl)
# colData(sample_exp)$CC.Phase <- cyc.values$phases
# colData(sample_exp)$CC.G1 <- cyc.values$normalized.scores$G1
# colData(sample_exp)$CC.S <- cyc.values$normalized.scores$S
# colData(sample_exp)$CC.G2M <- cyc.values$normalized.scores$G2M  
# plotReducedDim(sample_exp, dimred="UMAP", colour_by = "CC.Phase") +
#   ggsave(paste0(i, "/UMAP_with_cell.cycle.phase_", i, ".pdf"))
# plotReducedDim(sample_exp, dimred="PCA", colour_by = "CC.Phase") +
#   ggsave(paste0(i, "/PCA_with_cell.cycle.phase_", i, ".pdf"))
# plotReducedDim(sample_exp, dimred="PHATE", colour_by = "CC.Phase") +
#   ggsave(paste0(i, "/PHATE_with_cell.cycle.phase_", i, ".pdf"))
# plotReducedDim(sample_exp, dimred="DiffusionMap", colour_by = "CC.Phase") +
#   ggsave(paste0(i, "/DiffusionMap_with_cell.cycle.phase_", i, ".pdf"))

plotReducedDim(sample_exp, dimred="PCA", colour_by = "CDK1") +
  ggsave(paste0(i, "/PCA_with_CDK1_", i, ".pdf"))
plotReducedDim(sample_exp, dimred="UMAP", colour_by="CDK1") +
  ggsave(paste0(i, "/UMAP_with_CDK1_", i, ".pdf"))
plotReducedDim(sample_exp, dimred="DiffusionMap", colour_by = "CDK1") +
  ggsave(paste0(i, "/DiffusionMap_with_CDK1_", i, ".pdf"))
plotReducedDim(sample_exp, dimred="PHATE", colour_by = "CDK1") +
  ggsave(paste0(i, "/PHATE_with_CDK1_", i, ".pdf"))

saveRDS(sample_exp, paste0(i,"/Prefiltered_experiment_with_clusters_", i, ".rds"))
rm(sample_exp)
}

