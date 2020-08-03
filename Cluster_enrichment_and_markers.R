#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Enrichment analysis of clusters identified by MELD
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/meld/all_data")
  place <- "wolfpack"
  set.seed(100)
  options(future.globals.maxSize = 200000*1024^2)
}

# Libraries
library('Seurat')
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('batchelor')
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
library('msigdbr')
library('fgsea')
library('clusterProfiler')

# Load Seurat object and metadata generated in python.
# --------------------------------------------------------------------------

if (place == 'local'){
  metadata_MELD <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/practice_all_data/Prefiltered_experiment_practice_seurat_integrated_MELD_metadata_final_clusters.csv", header = TRUE, row.names = 1)
} else {
  metadata_MELD <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/all_data/Prefiltered_experiment_up_seurat_integrated_MELD_metadata_final_clusters.csv", header = TRUE, row.names = 1)
}

if(place == "local") {
  filtered_exp.integrated <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Prefiltered_experiment_practice_seurat_integrated_for_MELD.rds") # uses practice data if local
} else {
  filtered_exp.integrated <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/Prefiltered_experiment_up_seurat_integrated_for_MELD.rds") # uses whole dataset if wolfpack
}

if(place == "local") {
  rowMetaData <- readRDS("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/practice_all_data/Feature_metadata.rds") # uses practice data if local
} else {
  rowMetaData <- readRDS("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/seurat/all_data/Feature_metadata.rds") # uses whole dataset if wolfpack
}

if(place == "local") {
imputed_MAGIC <- read.csv("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/meld/practice_all_data/Prefiltered_experiment_practice_seurat_integrated_MAGIC_counts.csv", header = TRUE, row.names = 1) # uses practice data if local
} else {
imputed_MAGIC <- read.csv("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/meld/all_data/Prefiltered_experiment_all_seurat_integrated_MAGIC_counts.csv", header = TRUE, row.names = 1) # uses whole dataset if wolfpack
}

# Plot relative abundance of cells from each tissue type in each cluster
# --------------------------------------------------------------------------

# Colour scheme for tissues
tissue_cmap = c('Primary' = '#E69F00', 'LN' = '#56B4E9', 'Liver' = '#009E73', 'Lung' = '#F0E442')

# Number of cells in each cluster by tissue
cell_states <- data.frame(table(metadata_MELD$Cluster_final, metadata_MELD$Tissue))
cell_states <- spread(data = cell_states, key = "Var2", value = "Freq")
rownames(cell_states) <- cell_states[,1]
cell_states <- cell_states[,-1]

# Normalise counts by number of cells in tissue
total_sample <- colSums(cell_states)
weighted_cell_states <- apply(cell_states, 1, FUN = function(x){x/total_sample})

# Plot of cluster composition by tissue
total_cell_stated <- colSums(weighted_cell_states)
percent_cell_states <- data.frame(apply(weighted_cell_states, 1, FUN = function(x){x/total_cell_stated}))
percent_cell_states_long <- gather(percent_cell_states)
percent_cell_states_long$cluster <- factor(rep(rownames(percent_cell_states), 4), levels = rownames(percent_cell_states))
percent_cell_states_long$tissue <- factor(percent_cell_states_long$key, levels = rev(c("Primary", "LN", "Liver", "Lung")))

cluster_comp <- ggplot(percent_cell_states_long, aes(x = value, y = cluster, fill = tissue)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = tissue_cmap) +
  theme_classic() +
  coord_polar("y", start=0) +
  ggsave("Cluster_composition_by_weighted_tissue.pdf")

# Plot relative abundance of cells from each tissue type in each cluster
# --------------------------------------------------------------------------

# Colour scheme for clusters
cluster_cmap <- hcl.colors(n = length(unique(metadata_MELD$Cluster_final)), palette = "Spectral")
names(cluster_cmap) <- unique(metadata_MELD$Cluster_final)

cell_states <- data.frame(table(metadata_MELD$Tissue, metadata_MELD$Cluster_final))
cell_states <- spread(data = cell_states, key = "Var2", value = "Freq")
rownames(cell_states) <- cell_states[,1]
cell_states <- cell_states[,-1]

total_sample <- rowSums(cell_states)

percent_cell_states <- data.frame(apply(cell_states, 2, FUN = function(x){x/total_sample}))
colnames(percent_cell_states) <- colnames(cell_states)
percent_cell_states_long <- gather(percent_cell_states)
percent_cell_states_long$cluster <- percent_cell_states_long$key
percent_cell_states_long$tissue <- factor(rep(rownames(cell_states), length(colnames(percent_cell_states))), levels = rev(c("Primary", "LN", "Liver", "Lung")))

tissue_comp <- ggplot(percent_cell_states_long, aes(x = value, y = tissue, fill = cluster)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = cluster_cmap) +
  theme_classic() +
  coord_polar("y", start=0) +
  ggsave("Tissue_composition_by_cluster.pdf")


tissue_dot <- ggplot(percent_cell_states_long, aes(x = tissue, y = cluster)) +
  geom_point(aes(size = value, fill = cluster), shape = 21, color = "black", show.legend = c('size' = TRUE, 'fill' = FALSE)) +
  scale_radius(range = c(0, 15), limits = c(0, max(percent_cell_states_long$value))) +
  scale_fill_manual(values = cluster_cmap) +
  coord_flip() +
  theme_minimal() +
  ggsave("Tissue_composition_by_cluster_dotplot.pdf", useDingbats = FALSE)

# Convert to SCE object to perform differential marker detection.
# --------------------------------------------------------------------------
DefaultAssay(filtered_exp.integrated) <- "RNA"
filtered_exp.sce <- as.SingleCellExperiment(filtered_exp.integrated)
filtered_exp.sce <- filtered_exp.sce[rownames(filtered_exp.integrated@assays$integrated@scale.data),]

idx <- match(rownames(colData(filtered_exp.sce)), rownames(metadata_MELD))
filtered_exp.sce$MELD_final_clusters <- as.factor(metadata_MELD$Cluster_final) [idx]

reducedDim(filtered_exp.sce, "PHATE") <- metadata_MELD[,grepl("PHATE", colnames(metadata_MELD))]
assay(filtered_exp.sce, "SCT_integrated") <- filtered_exp.integrated@assays$integrated@scale.data
assay(filtered_exp.sce, "MAGIC_integrated") <- t(imputed_MAGIC)

# Log normalise counts, scaled by library size, across batches.
filtered_exp.sce <- multiBatchNorm(filtered_exp.sce, batch = filtered_exp.sce$Replicate, normalize.all = TRUE)

# Idenfity markers that are upregulated and down regulated in each cluster
# --------------------------------------------------------------------------

# Plot colored by final clusters
final_clusters <- plotReducedDim(filtered_exp.sce, dimred = "PHATE", colour_by = "MELD_final_clusters") +
  scale_fill_manual(values = cluster_cmap) +
  ggsave("Final_MELD_clusters.pdf")


# Plot colored by tissue
tissues <- plotReducedDim(filtered_exp.sce, dimred = "PHATE", colour_by = "Tissue") +
  scale_fill_manual(values = tissue_cmap)

# Find the top genes based on the median of their scaled expression
top_median <- function(x){
  tmp <- assay(filtered_exp.sce[,filtered_exp.sce$MELD_final_clusters == x], "SCT_integrated")
  rm_tmp <- rowMedians(tmp)
  names(rm_tmp) <- rownames(tmp)
  rm_tmp <- rm_tmp[order(rm_tmp, decreasing = TRUE)]
  rm_tmp <- rm_tmp[1:50]
  return(rm_tmp)
}

top_genes <- lapply(as.list(as.character(unique(filtered_exp.sce$MELD_final_clusters))), top_median)
names(top_genes) <- as.character(unique(filtered_exp.sce$MELD_final_clusters))

# Cluster samples based on median scaled expression profiles
cluster_median <- function(x){
  tmp <- assay(filtered_exp.sce[,filtered_exp.sce$MELD_final_clusters == x], "SCT_integrated")
  rm_tmp <- rowMedians(tmp)
  return(rm_tmp)
}
cluster_medians <- lapply(as.list(as.character(unique(filtered_exp.sce$MELD_final_clusters))), cluster_median)
names(cluster_medians) <- as.character(unique(filtered_exp.sce$MELD_final_clusters))
plot(hclust(dist(t(as.data.frame(cluster_medians))), method = "median"))

# Idenitify genes that are significantly up or down regulated in each cluster
up_markers <- findMarkers(filtered_exp.sce, groups = filtered_exp.sce$MELD_final_clusters, test.type = "wilcox", full.stats = TRUE, pval.type = "some", direction = "up", block = filtered_exp.sce$Replicate, lfc = 0.25, min.prop = 0.7)
down_markers <- findMarkers(filtered_exp.sce, groups = filtered_exp.sce$MELD_final_clusters, test.type = "wilcox", full.stats = TRUE, pval.type = "some", direction = "down", block = filtered_exp.sce$Replicate, lfc = 0.25, min.prop = 0.7)

# Find genes with the highest cumulative AUC relative to other clusters 
cumulative_AUC <-  function(x, y){
  tmp <- data.frame(x[[y]])
  tmp_auc <- tmp[,grepl("AUC", colnames(tmp))]
  cum_auc <- rowSums(tmp_auc)
  return(cum_auc)
}

for(i in names(up_markers)){
  print(i)
  interesting_up <- data.frame(up_markers[[i]])
  interesting_up$cum_auc <- cumulative_AUC(up_markers, i)
  interesting_up <- interesting_up[order(interesting_up$cum_auc, decreasing = T),]
  interesting_up <- interesting_up[1:10, ]
  print("Up")
  print(data.frame("Gene" = rownames(interesting_up), "Cum_AUC" = interesting_up$cum_auc))

  for(gene in rownames(interesting_up)) {
    prd <- plotReducedDim(filtered_exp.sce[,order(logcounts(filtered_exp.sce[gene,]))], dimred = "PHATE", colour_by = gene, by_exprs_values = "MAGIC_integrated") +
      # scale_fill_gradientn(colors = hcl.colors(n = 7, palette = "TealRose"), limits = c(-3, 3), oob = scales::squish) # Uncomment if using SCT values
      scale_fill_gradientn(colors = hcl.colors(n = 7, palette = "TealRose"), oob = scales::squish)
    pe <- plotExpression(filtered_exp.sce, gene, x = "MELD_final_clusters")
    tpe <- plotExpression(filtered_exp.sce, gene, x = "Tissue")
    mp <- plot_grid(final_clusters, cluster_comp, tissue_dot, prd, pe, tpe, ncol = 2)
    ggsave(paste0("MELD_cluster_", i, "_upregulated_marker_", gene, ".png"), plot = mp, device = "png")
  }
  
  interesting_down <- data.frame(down_markers[[i]])
  interesting_down$cum_auc <- cumulative_AUC(down_markers, i)
  interesting_down <- interesting_down[order(interesting_down$cum_auc, decreasing = T),]
  interesting_down <- interesting_down[1:10, ]
  print("Down")
  print(data.frame("Gene" = rownames(interesting_down), "Cum_AUC" = interesting_down$cum_auc))

  for(gene in rownames(interesting_down)) {
    prd <- plotReducedDim(filtered_exp.sce[,order(logcounts(filtered_exp.sce[gene,]))], dimred = "PHATE", colour_by = gene, by_exprs_values = "SCT_integrated") +
      scale_fill_gradientn(colors = hcl.colors(n = 7, palette = "TealRose"), limits = c(-3, 3), oob = scales::squish)
    pe <- plotExpression(filtered_exp.sce, gene, x = "MELD_final_clusters")
    tpe <- plotExpression(filtered_exp.sce, gene, x = "Tissue")
    mp <- plot_grid(final_clusters, cluster_comp, tissue_dot, prd, pe, tpe, ncol = 2)
    ggsave(paste0("MELD_cluster_", i, "_downregulated_marker_", gene, ".png"), plot = mp, device = "png")
  }
}

# Perform pathway analysis on marker genes in each cluster
# --------------------------------------------------------------------------

# Make gene universe(s)
universe <- as.character(unique(rowMetaData[rowMetaData$Any_Active, "Ensembl"]))
universe_entrez <- as.character(unique(rowMetaData[rowMetaData$Any_Active, "EntrezID"]))
universe_genesymbol <- as.character(unique(rowMetaData[rowMetaData$Any_Active, "GeneSymbol"]))

# Makes hallmark geneset for enrichment testing
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_t2g <- h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()

# Makes hallmark geneset for enrichment testing
metabolic_pathways <- read.csv("~/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Metabolic_pathways_genes.csv", header = TRUE)
metabolic_pathways$Gene <- as.character(metabolic_pathways$Gene)
metabolic_pathways$Metabolic_pathway <- as.character(metabolic_pathways$Metabolic_pathway)

# Perform enrichment analysis on markers up-regulated markers in each cluster
# --------------------------------------------------------------------------

dir.create("markers")

for(i in names(up_markers)){
  dir.create(paste0("markers/MELD_", i))
  interesting_up <- data.frame(up_markers[[i]])
  interesting_up$cum_auc <- cumulative_AUC(up_markers, i)
  interesting_up <- interesting_up[order(interesting_up$cum_auc, decreasing = T),]
  interesting_up <- interesting_up[1:50, ]
  idx <- match(rownames(interesting_up), rowMetaData$Seurat_IDs)
  interesting_up$Ensembl <- rowMetaData$Ensembl [idx]
  interesting_up$EntrezID <- rowMetaData$EntrezID [idx]
  interesting_up$GeneSymbol <- rowMetaData$GeneSymbol [idx]

  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting_up[, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/MELD_", i,"/GOBP_markers_up_MELD_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting_up[, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/MELD_", i,"/GOMF_markers_up_MELD_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting_up[, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/MELD_", i,"/GOCC_markers_up_MELD_", i, ".csv"))
  
  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting_up[, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = universe_genesymbol)
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("markers/MELD_", i,"/HALLMARK_markers_up_MELD_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting_up[, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/MELD_", i,"/KEGG_markers_up_MELD_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting_up[, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/MELD_", i,"/Reactome_markers_up_MELD_", i, ".csv"))

  # Perform Metabolic_pathways enrichment analysis
  Metabolic_pathwaysEnrich <- enricher(gene = interesting_up[, "GeneSymbol"], 
                                       TERM2GENE = metabolic_pathways,
                                       pvalueCutoff = 1,
                                       qvalueCutoff = 1,
                                       pAdjustMethod = "bonferroni",
                                       minGSSize = 10,
                                       maxGSSize = 500,
                                       universe = universe_genesymbol)
  Metabolic_pathways_GOI <- as.data.frame(Metabolic_pathwaysEnrich@result)
  write.csv(Metabolic_pathways_GOI, paste0("markers/MELD_", i,"/Metabolic_pathways_up_MELD_", i, ".csv"))
}

# Perform enrichment analysis on markers down-regulated markers in each cluster
# --------------------------------------------------------------------------

dir.create("markers")

for(i in names(down_markers)){
  dir.create(paste0("markers/MELD_", i))
  interesting_down <- data.frame(down_markers[[i]])
  interesting_down$cum_auc <- cumulative_AUC(down_markers, i)
  interesting_down <- interesting_down[order(interesting_down$cum_auc, decreasing = T),]
  interesting_down <- interesting_down[1:50, ]
  idx <- match(rownames(interesting_down), rowMetaData$Seurat_IDs)
  interesting_down$Ensembl <- rowMetaData$Ensembl [idx]
  interesting_down$EntrezID <- rowMetaData$EntrezID [idx]
  interesting_down$GeneSymbol <- rowMetaData$GeneSymbol [idx]
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting_down[, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "BP",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOBP_GOI <- as.data.frame(BPenrich@result)
  write.csv(GOBP_GOI, paste0("markers/MELD_", i,"/GOBP_markers_down_MELD_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting_down[, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "MF",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOMF_GOI <- as.data.frame(MFenrich@result)
  write.csv(GOMF_GOI, paste0("markers/MELD_", i,"/GOMF_markers_down_MELD_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting_down[, "Ensembl"],
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENSEMBL",
                       ont = "CC",
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "bonferroni",
                       minGSSize = 10,
                       maxGSSize = 500,
                       readable = TRUE,
                       universe = universe)
  GOCC_GOI <- as.data.frame(CCenrich@result)
  write.csv(GOCC_GOI, paste0("markers/MELD_", i,"/GOCC_markers_down_MELD_", i, ".csv"))
  
  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting_down[, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = universe_genesymbol)
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("markers/MELD_", i,"/HALLMARK_markers_down_MELD_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting_down[, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/MELD_", i,"/KEGG_markers_down_MELD_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting_down[, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/MELD_", i,"/Reactome_markers_down_MELD_", i, ".csv"))
  
  # Perform Metabolic_pathways enrichment analysis
  Metabolic_pathwaysEnrich <- enricher(gene = interesting_down[, "GeneSymbol"], 
                                       TERM2GENE = metabolic_pathways,
                                       pvalueCutoff = 1,
                                       qvalueCutoff = 1,
                                       pAdjustMethod = "bonferroni",
                                       minGSSize = 10,
                                       maxGSSize = 500,
                                       universe = universe_genesymbol)
  Metabolic_pathways_GOI <- as.data.frame(Metabolic_pathwaysEnrich@result)
  write.csv(Metabolic_pathways_GOI, paste0("markers/MELD_", i,"/Metabolic_pathways_down_MELD_", i, ".csv"))
}

