#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Characterise cluster by SOX2 and OCT4, markers of cancer stem cells http://dx.doi.org/10.1016/j.stemcr.2014.11.002
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results_W_strict_filtering_DEP/seurat/all_data")
  place <- "wolfpack"
}

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('ggridges')
library('readr')
library('Matrix')
library('viridis')
library('reshape2')
library('edgeR')
library('Rmagic')

set.seed(100)
# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_experiment_practice_seurat_merge_cluster_sce.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_experiment_all_seurat_merge_cluster_sce.rds") # uses whole dataset if wolfpack
}

# Make dataframe used for plotting CSC markers
# --------------------------------------------------------------------------

plotting <- data.frame(reducedDim(filtered_exp, "PHATE_seurat"))
plotting$Tissue <- filtered_exp$Tissue
plotting$cluster <- filtered_exp$seurat_clusters

# Quantitative expression
plotting$OCT4_expression <- logcounts(filtered_exp)["POU5F2",]
plotting$SOX2_expression <- logcounts(filtered_exp)["SOX2",]

# Scaled expression
plotting$OCT4_zscore <- scale(plotting$OCT4_expression)
plotting$SOX2_zscore <- scale(plotting$SOX2_expression)
plotting$OCT4ySOX2_zscore <- plotting$OCT4_expression + plotting$SOX2_expression

# Qualitative expression
plotting$OCT4_expressed <- plotting$OCT4_expression > 0
plotting$SOX2_expressed <- plotting$SOX2_expression > 0
plotting$OCT4ySOX2_expressed <- plotting$OCT4_expressed & plotting$SOX2_expressed

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
  ggsave("select_genes/All_cells_OCT4ySOX2_expressed_all.pdf", useDingbats = FALSE)

plotting %>%
  select(cluster, OCT4ySOX2_expressed) %>%
  group_by(cluster) %>%
  summarize(OCT4ySOX2=sum(OCT4ySOX2_expressed)) %>%
  ggplot(aes(x=cluster, y=OCT4ySOX2)) +
  geom_bar(stat="identity", color="black") +
  ggtitle("OCT4ySOX2_expressed") +
  scale_color_viridis(discrete=TRUE) +
  coord_flip() +
  theme_classic() +
  ggsave("select_genes/All_cells_OCT4ySOX2_expressed_all_barplot.pdf", useDingbats = FALSE)

# Make dataframe used for plotting EMT sorting markers
# --------------------------------------------------------------------------

# Quantitative expression
plotting$CD44_expression <- logcounts(filtered_exp)["CD44",]
plotting$ITGB4_expression <- logcounts(filtered_exp)["ITGB4",]

# Scaled expression
plotting$CD44_zscore <- scale(plotting$CD44_expression)
plotting$ITGB4_zscore <- scale(plotting$ITGB4_expression)
plotting$CD44yITGB4_zscore <- plotting$CD44_expression + plotting$ITGB4_expression

# Qualitative expression
plotting$CD44_expressed <- plotting$CD44_expression > 0
plotting$ITGB4_expressed <- plotting$ITGB4_expression > 0
plotting$CD44yITGB4_expressed <- plotting$CD44_expressed & plotting$ITGB4_expressed

# Cluster labelling
text_out <- retrieveCellInfo(filtered_exp, "seurat_clusters", search="colData")
text_out$val <- as.factor(text_out$val)
by_text_x <- vapply(split(plotting$PHATE1, text_out$val), median, FUN.VALUE=0)
by_text_y <- vapply(split(plotting$PHATE2, text_out$val), median, FUN.VALUE=0)

# Plot cells within each cluster
# --------------------------------------------------------------------------

plotting_t <- plotting[plotting$CD44yITGB4_expressed == TRUE, ]
ggplot() +
  geom_point(data = plotting, aes(x = PHATE1, y = PHATE2), color = "grey80") +
  geom_point(data = plotting_t, aes(x = PHATE1, y = PHATE2, color = CD44yITGB4_zscore)) +
  scale_color_viridis_c(option = "C") +
  theme(aspect.ratio = 1) +
  theme_classic() +
  annotate("text", x = by_text_x, y = by_text_y, label = names(by_text_x)) +
  ggsave("select_genes/All_cells_CD44yITGB4_scaled_all.pdf", useDingbats = FALSE)

ggplot() +
  geom_point(data = plotting, aes(x = PHATE1, y = PHATE2, color = CD44yITGB4_expressed)) +
  scale_color_viridis_d(option = "D") +
  theme(aspect.ratio = 1) +
  theme_classic() +
  annotate("text", x = by_text_x, y = by_text_y, label = names(by_text_x)) +
  ggsave("select_genes/All_cells_CD44yITGB4_expressed_all.pdf", useDingbats = FALSE)

CD44yITGB4_by_cluster <- table(cluster=plotting$cluster, CD44yITGB4=plotting$CD44yITGB4_expressed)
round(CD44yITGB4_by_cluster/rowSums(CD44yITGB4_by_cluster)*100)

plotting %>%
  select(cluster, CD44yITGB4_expressed) %>%
  group_by(cluster) %>%
  summarize(CD44yITGB4=sum(CD44yITGB4_expressed)/n(cluster)*100) %>%
  ggplot(aes(x=cluster, y=CD44yITGB4)) +
  geom_bar(stat="identity", color="black") +
  ggtitle("CD44yITGB4_expressed") +
  scale_color_viridis(discrete=TRUE) +
  coord_flip() +
  theme_classic() +
  ggsave("select_genes/All_cells_CD44yITGB4_expressed_all_barplot.pdf", useDingbats = FALSE)








