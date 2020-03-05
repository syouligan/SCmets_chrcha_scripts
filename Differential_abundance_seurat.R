#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Normalise, batch-correct and cluster cells
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
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('readr')
library('Matrix')
library('viridis')
library('reshape2')
library('edgeR')

set.seed(100)
# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_experiment_practice_seurat_merge_cluster_sce.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_experiment_all_seurat_merge_cluster_sce.rds") # uses whole dataset if wolfpack
}

# Summarize the number of cells per tissue for each cluster
cell_clusters <- table(Cluster=filtered_exp$seurat_clusters, Sample=filtered_exp$Sample)
liver_sum <- rowSums(cell_clusters[,1:4])
ln_sum <- rowSums(cell_clusters[,5:8])
lung_sum <- rowSums(cell_clusters[,9:12])
primary_sum <- rowSums(cell_clusters[,13:16])
cell_sum <- data.frame(Primary = primary_sum, Liver = liver_sum, Lymph_node = ln_sum, Lung = lung_sum)
cell_percent <- cell_sum/rowSums(cell_sum)*100
cell_percent$Cluster <- rownames(cell_clusters)
cell_percent <- melt(cell_percent)
colnames(cell_percent) <- c("Cluster", "Tissue", "Percent")

# Make cumulative bargraphs for the number of cells per tissue each cluster
cell_percent$Cluster <- factor(cell_percent$Cluster, levels = rev(c(0:(length(unique(tissue_percent$Cluster))-1))))
cell_percent$Tissue <- factor(cell_percent$Tissue, levels = rev(c("Primary", "Lung", "Lymph_node", "Liver")))
cell_percent <- arrange(cell_percent, Tissue, Cluster) 

ggplot(data=cell_percent, aes(x=Cluster, y=Percent, fill=Tissue)) +
  geom_bar(stat="identity") +
  scale_fill_viridis_d(option = "C") +
  coord_flip() +
  theme_classic() +
  ggsave("cell_composition/Total_cumsum_cluster_counts.pdf", useDingbats = FALSE)

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
    ggsave(paste0("cell_composition/", i, "_cluster_counts.pdf"), useDingbats = FALSE)
}

tissue_percent <- data.frame(t(apply(cell_sum, 1, function(x) x/colSums(cell_sum)*100)))
tissue_percent$Cluster <- rownames(tissue_percent)
tissue_percent <- melt(tissue_percent)
colnames(tissue_percent) <- c("Cluster", "Tissue", "Percent")

tissue_percent$Cluster <- factor(tissue_percent$Cluster, levels = rev(c(0:(length(unique(tissue_percent$Cluster))-1))))
tissue_percent$Tissue <- factor(tissue_percent$Tissue, levels = rev(c("Primary", "Lung", "Lymph_node", "Liver")))
tissue_percent <- arrange(tissue_percent, Tissue, Cluster) 

ggplot(data=tissue_percent, aes(x=Tissue, y=Percent, fill=Cluster)) +
  geom_bar(stat="identity") +
  scale_fill_viridis_d(option = "D") +
  coord_flip() +
  theme_classic() +
  ggsave("cell_composition/Total_cumsum_tissue_counts.pdf", useDingbats = FALSE)

# Make individual bargraphs for the number of cells per tissue each cluster
for(i in as.character(unique(tissue_percent$Tissue))) {
  cell_tissue <- tissue_percent[tissue_percent$Tissue == i, ]
  
  ggplot(data=cell_tissue, aes(x=Cluster, y=Percent)) +
    geom_bar(stat="identity") +
    scale_fill_viridis_c(option = "C") +
    ylim(0, 100) +
    coord_flip() +
    theme_classic() +
    theme(aspect.ratio = 1) +
    ggsave(paste0("cell_composition/", i, "_tissue_counts.pdf"), useDingbats = FALSE)
}

# Plot cells from each tissue on total PHATE background
color <- viridis(4, option = "C")
names(color) <- as.character(unique(filtered_exp$Tissue))

# color_vec <- filtered_exp$Tissue
# for(i in as.character(unique(filtered_exp$Tissue))){
# color_vec[color_vec == i] <- color[[i]]
# } # Makes vector of colors according to tissue
# filtered_exp$Tissue_colors <- color_vec

plotting <- data.frame(reducedDim(filtered_exp, "PHATE_seurat"))
plotting$Tissue <- filtered_exp$Tissue
plotting$Mito_percent <- filtered_exp$Mito_percent
text_out <- retrieveCellInfo(filtered_exp, "seurat_clusters", search="colData")
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
  ggsave(paste0("cell_composition/All_cells_",i ,".pdf"), useDingbats = FALSE)
}

# Plot cells from each tissue on total PHATE background and based on mito-content

for(i in as.character(unique(filtered_exp$Tissue))) {
  plotting_t <- plotting[plotting$Tissue == i, ]
  ggplot() +
    geom_point(data = plotting, aes(x = PHATE1, y = PHATE2), color = "grey80") +
    geom_point(data = plotting_t, aes(x = PHATE1, y = PHATE2, color = Mito_percent)) +
    scale_color_viridis_c(option = "B") +
    theme(aspect.ratio = 1) +
    theme_classic() +
    ggsave(paste0("cell_composition/All_cells_Mito_",i ,".pdf"), useDingbats = FALSE)
}

# Calculate differential abundance of cells in each cluster between tissues
cell_clusters <- unclass(cell_clusters) 
head(cell_clusters)
extra.info <- colData(filtered_exp)[match(colnames(cell_clusters), filtered_exp$Sample),]
extra.info <- extra.info[,c("Tissue", "Replicate")]
y.ab <- DGEList(cell_clusters, samples=extra.info)
tissue <- factor(extra.info$Tissue, levels = c("Primary", "Liver", "LN", "Lung"))
repli <- factor(extra.info$Replicate)
design <- model.matrix(~0 + tissue + repli, y.ab$Samples)
y.ab <- estimateDisp(y.ab, design, trend="none")
fit.ab <- glmQLFit(y.ab, design, robust=TRUE, abundance.trend=FALSE)
my.contrasts <- makeContrasts(LivervsPrimary=tissueLiver-tissuePrimary, LNvsPrimary=tissueLN-tissuePrimary, LungvsPrimary=tissueLung-tissuePrimary, LungvsLN=tissueLung-tissueLN, LivervsLN=tissueLiver-tissueLN, LungvsLiver=tissueLung-tissueLiver, levels=design)

for(i in colnames(my.contrasts)) {
  res <- glmQLFTest(fit.ab, contrast=my.contrasts[,i])
  summary(decideTests(res))
  write.csv(data.frame(topTags(res, n = Inf)), paste0("cell_composition/Differential_abundance_", i,".csv"))
  }



