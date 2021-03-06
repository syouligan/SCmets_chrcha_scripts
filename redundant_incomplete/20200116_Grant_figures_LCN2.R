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
library('viridis')
library('reshape2')

set.seed(100)
# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_experiment_Practice_merge_cluster.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_experiment_All_merge_cluster.rds") # uses whole dataset if wolfpack
}

# Summarize the number of cells per tissue for each cluster
cell_clusters <- read.csv("Corrected_batch_cell_cluster_membership.csv", row.names = 1)
liver_sum <- rowSums(cell_clusters[,1:5])
ln_sum <- rowSums(cell_clusters[,6:10])
lung_sum <- rowSums(cell_clusters[,11:15])
primary_sum <- rowSums(cell_clusters[,16:19])
cell_sum <- data.frame(Primary = primary_sum, Liver = liver_sum, Lymph_node = ln_sum, Lung = lung_sum)
cell_percent <- cell_sum/rowSums(cell_sum)*100
cell_percent$Cluster <- rownames(cell_clusters)
cell_percent <- melt(cell_percent)
colnames(cell_percent) <- c("Cluster", "Tissue", "Percent")

# Make cumulative bargraphs for the number of cells per tissue each cluster
cell_percent$Cluster <- factor(cell_percent$Cluster, levels = rev(c(1:length(unique(cell_percent$Cluster)))))
cell_percent$Tissue <- factor(cell_percent$Tissue, levels = rev(c("Primary", "Lung", "Lymph_node", "Liver")))
cell_percent <- arrange(cell_percent, Tissue, Cluster) 

ggplot(data=cell_percent, aes(x=Cluster, y=Percent, fill=Tissue)) +
  geom_bar(stat="identity") +
  scale_fill_viridis_d(option = "C") +
  coord_flip() +
  theme_classic() +
  ggsave("Total_cumsum_cluster_counts.pdf", useDingbats = FALSE)

# Make individual bargraphs for the number of cells per tissue each cluster
for(i in as.character(unique(cell_percent$Tissue))) {
  cell_tissue <- cell_percent[cell_percent$Tissue == i, ]
  
  ggplot(data=cell_tissue, aes(x=Cluster, y=Percent)) +
    geom_bar(stat="identity") +
    scale_fill_viridis_c(option = "C") +
    ylim(0, 100) +
    coord_flip() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    ggsave(paste0(i, "_cluster_counts.pdf"), useDingbats = FALSE)
}

# Plot cells from each tissue on total PHATE background
color <- viridis(4, option = "C")
names(color) <- as.character(unique(filtered_exp$Tissue))

# color_vec <- filtered_exp$Tissue
# for(i in as.character(unique(filtered_exp$Tissue))){
# color_vec[color_vec == i] <- color[[i]]
# } # Makes vector of colors according to tissue
# filtered_exp$Tissue_colors <- color_vec

plotReducedDim(filtered_exp, dimred="PHATE_fastMNN", colour_by = "Tissue", text_by = "cluster" ) +
  scale_fill_viridis_d(option = "C") +
  ggsave("PHATE_corrected_with_clusters_tissue.pdf")

plotting <- data.frame(reducedDim(filtered_exp, "PHATE_fastMNN"))
plotting$Tissue <- filtered_exp$Tissue
plotting$Mito_percent <- filtered_exp$Mito_percent
text_out <- retrieveCellInfo(filtered_exp, "cluster", search="colData")
text_out$val <- as.factor(text_out$val)
by_text_x <- vapply(split(plotting$PHATE1, text_out$val), median, FUN.VALUE=0)
by_text_y <- vapply(split(plotting$PHATE2, text_out$val), median, FUN.VALUE=0)

for(i in as.character(unique(filtered_exp$Tissue))) {
  plotting_t <- plotting[plotting$Tissue == i, ]
  ggplot() +
    geom_point(data = plotting, aes(x = PHATE1, y = PHATE2), color = "grey80") +
    geom_point(data = plotting_t, aes(x = PHATE1, y = PHATE2), color = color[[i]], alpha = 0.5) +
    theme(aspect.ratio = 1) +
    theme_classic() +
    annotate("text", x = by_text_x, y = by_text_y, label = names(by_text_x))
    ggsave(paste0("All_cells_",i ,".pdf"), useDingbats = FALSE)
}

# Plot cells from each tissue on total PHATE background and based on mito-content
ggplot() +
  geom_point(data = plotting, aes(x = PHATE1, y = PHATE2, color = Mito_percent)) +
  scale_color_viridis_c(option = "B") +
  theme(aspect.ratio = 1) +
  theme_classic() +
  ggsave("All_cells_Mito_content.pdf", useDingbats = FALSE)

for(i in as.character(unique(filtered_exp$Tissue))) {
  plotting_t <- plotting[plotting$Tissue == i, ]
  ggplot() +
    geom_point(data = plotting, aes(x = PHATE1, y = PHATE2), color = "grey80") +
    geom_point(data = plotting_t, aes(x = PHATE1, y = PHATE2, color = Mito_percent)) +
    scale_color_viridis_c(option = "B") +
    theme(aspect.ratio = 1) +
    theme_classic() +
    ggsave(paste0("All_cells_Mito_",i ,".pdf"), useDingbats = FALSE)
}

# Find genes positive for tumor intiating cell (TI cell) markers in different combinations
TIC_groups <- list(TIC1 = c("ZEB1", "CD44", "AR"),
     TIC2 = c("ZEB1", "CD44", "ITGB4", "AR"),
     TIC3 = c("ZEB1", "CD44", "LCN2"),
     TIC4 = c("ZEB1", "CD44", "LCN2", "ITGB4")
     )

for(i in 1:length(TIC_groups)) {
  name <- names(TIC_groups[i])
  GOI.sce <- filtered_exp[is.element(rowData(filtered_exp)$GeneSymbol, TIC_groups[[i]]),]
  rowSums(counts(GOI.sce) > 0)
  GOI.sce$TI_cells <- colSums(counts(GOI.sce) > 0) == length(TIC_groups[[i]])
  TI_cells_tab <- table(Cluster=GOI.sce$TI_cells, Batch=GOI.sce$Sample)
  write.csv(TI_cells_tab, paste0(name, "_cells_counts_per_sample.csv"))
  TI_cells_tissues_tab <- table(Cluster=GOI.sce$TI_cells, Batch=GOI.sce$Tissue)
  write.csv(TI_cells_tissues_tab, paste0(name, "_cells_counts_per_tissue.csv"))
  
  # Plot number of TI cells per sample
  TI_cells_tab <- data.frame(Sample = colnames(TI_cells_tab), TI_cells = TI_cells_tab["TRUE",])
  ggplot(data=TI_cells_tab, aes(x=Sample, y=TI_cells)) +
    geom_bar(stat="identity", color="black") +
    ggtitle(paste0(name," Tumor intiating cells")) +
    scale_color_viridis(discrete=TRUE) +
    coord_flip() +
    theme_classic() +
    ggsave(paste0(name, "_cells_counts.pdf"), useDingbats = FALSE)
  
  # Plot number of TI cells per tissue
  TI_cells_tissues_tab <- data.frame(Tissue = colnames(TI_cells_tissues_tab), TI_cells = TI_cells_tissues_tab["TRUE",])
  ggplot(data=TI_cells_tissues_tab, aes(x=Tissue, y=TI_cells)) +
    geom_bar(stat="identity", color="black") +
    ggtitle(paste0(name," Tumor intiating cells")) +
    scale_color_viridis(discrete=TRUE) +
    coord_flip() +
    theme_classic() +
    ggsave(paste0(name, "_cells_counts_per_tissue.pdf"), useDingbats = FALSE)
  
  # Plot the TI cell score
  GOI.sce$Total_TI <- colSums(counts(GOI.sce))
  plotting <- data.frame(reducedDim(GOI.sce, "PHATE_fastMNN"))
  plotting$TI_cells <- GOI.sce$TI_cells
  plotting$Total_TI <- GOI.sce$Total_TI
  plotting$Tissue <- GOI.sce$Tissue
  plotting_w <- plotting[plotting$TI_cells == TRUE, ]
  plotting_wo <- plotting[plotting$TI_cells == FALSE, ]
  
  ggplot() +
    geom_point(data = plotting_wo, aes(x = PHATE1, y = PHATE2), color = "grey80") +
    geom_point(data = plotting_w, aes(x = PHATE1, y = PHATE2, color = Total_TI)) +
    scale_color_viridis_c(option = "B") +
    theme(aspect.ratio = 1) +
    theme_classic() +
    ggsave(paste0(name, "_cells_total.pdf"), useDingbats = FALSE)
}

# Impute counts for ZEB1, CD44, LCN2 genes with MAGIC (NOT USED AS DIDNT HAVE TIME)
# data_MAGIC <- magic(Matrix::t(assay(filtered_exp, "logcounts")), genes=c("ZEB1", "CD44", "LCN2", "CD104", "CDH1", "SULF1", "AR"))
# filtered_exp$ZEB1_magic <- scale(data_MAGIC$result$ZEB1)
# filtered_exp$CD44_magic <- scale(data_MAGIC$result$CD44)
# filtered_exp$LCN2_magic <- scale(data_MAGIC$result$LCN2)
# filtered_exp$Total_TI <- filtered_exp$ZEB1_magic + filtered_exp$CD44_magic + filtered_exp$LCN2_magic

# plotting <- data.frame(reducedDim(filtered_exp, "PHATE_fastMNN"))
# plotting$TI_cells <- filtered_exp$TI_cells
# plotting$Total_TI <- filtered_exp$Total_TI
# plotting$Tissue <- filtered_exp$Tissue
# plotting_w <- plotting[plotting$TI_cells == TRUE, ]
# plotting_wo <- plotting[plotting$TI_cells == FALSE, ]
# 
# ggplot() +
#   geom_point(data = plotting_wo, aes(x = PHATE1, y = PHATE2), color = "grey80") +
#   geom_point(data = plotting_w, aes(x = PHATE1, y = PHATE2, color = Total_TI)) +
#   scale_color_viridis_c(option = "B") +
#   theme(aspect.ratio = 1) +
#   theme_classic() +
#   ggsave("TI_cells_total.pdf", useDingbats = FALSE)
