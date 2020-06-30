#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Plot heatmaps of DEGs between tissues
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/")
  place <- "wolfpack"
}

# Libraries
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('Matrix')
library('limma')
library('edgeR')
library('gtools')
library('gplots')
library('msigdbr')
library('org.Hs.eg.db')
library('fgsea')
library('data.table')
library('sva')
library('viridis')
library('dendextend')
library('clusterProfiler')
library('ReactomePA')
library('Matrix')
library('phateR')

set.seed(100)

# Load Pseudo-bulk counts and row data SingleCellExperiment
exprs <- readRDS("Pseudo-bulk_whole_experiment_all.rds") # uses whole dataset if wolfpack
rowData_summed <- read.csv("Pseudo-bulk_rowData_all.csv", header = TRUE, row.names = 1)

# Load DGE data
# --------------------------------------------------------------------------
# Make list of top 100 differentially expressed genes
liver_primary <- read.csv("Liver_Primary_DEG_0LFC.csv", header = TRUE)
liver_primary <- liver_primary[liver_primary$adj.P.Val < 0.05, ]
liver_primary <- liver_primary[order(liver_primary$t),]
liver_primary <- as.character(liver_primary[, "X"])

ln_primary <- read.csv("LN_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary <- ln_primary[ln_primary$adj.P.Val < 0.05, ]
ln_primary <- ln_primary[order(ln_primary$t),]
ln_primary <- as.character(ln_primary[, "X"])

lung_primary <- read.csv("Lung_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary <- lung_primary[lung_primary$adj.P.Val < 0.05, ]
lung_primary <- lung_primary[order(lung_primary$t),]
lung_primary <- as.character(lung_primary[, "X"])

liver_lung <- read.csv("Liver_Lung_DEG_0LFC.csv", header = TRUE)
liver_lung <- liver_lung[liver_lung$adj.P.Val < 0.05, ]
liver_lung <- liver_lung[order(liver_lung$t),]
liver_lung <- as.character(liver_lung[, "X"])

liver_ln <- read.csv("Liver_LN_DEG_0LFC.csv", header = TRUE)
liver_ln <- liver_ln[liver_ln$adj.P.Val < 0.05, ]
liver_ln <- liver_ln[order(liver_ln$t),]
liver_ln <- as.character(liver_ln[, "X"])

ln_lung <- read.csv("LN_Lung_DEG_0LFC.csv", header = TRUE)
ln_lung <- ln_lung[ln_lung$adj.P.Val < 0.05, ]
ln_lung <- ln_lung[order(ln_lung$t),]
ln_lung <- as.character(ln_lung[, "X"])

total_DGE <- as.character(unique(c(liver_primary, ln_primary, lung_primary, liver_lung, ln_lung, liver_ln)))

DGE_list <- list("Total" = total_DGE, "Liver_Primary" = liver_primary, "LN_Primary" = ln_primary, "Lung_Primary" = lung_primary, "Liver_Lung" = liver_lung, "LN_Lung" = ln_lung, "Liver_LN" = liver_ln)

# Make heatmap of DEGs between tissues and primary tumors
# --------------------------------------------------------------------------

# Find genes detected in all replicates of a given tissue type
active_liver <- rowSums(exprs[,1:4] >= 1) == 4
active_LN <- rowSums(exprs[,5:8] >= 1) == 4
active_lung <- rowSums(exprs[,9:12] >= 1) == 4
active_primary <- rowSums(exprs[,13:16] >= 1) == 4
geneActivity <- data.frame("Liver" = active_liver, "LN" = active_LN, "Lung" = active_lung, "Primary" = active_primary)

# Make sample info table
sampleName <- colnames(exprs)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTable <- data.frame(sampleName, sampleType, sampleRep)

# Make DGEList object and normalise using TMM
allDGEList <- DGEList(counts = exprs, group = sampleType)
allDGEList <- calcNormFactors(allDGEList, method = "TMM")

# Make design matrix, fit linear models, list comparisons
type <- factor(sampleTable$sampleType)
rep <- factor(sampleTable$sampleRep, levels = c(1:4))
design <- model.matrix(~0 + type)
rownames(design) <- sampleName
allDGEList <- voom(allDGEList, design, plot = FALSE)

exprs_no_batch <- removeBatchEffect(allDGEList$E, sampleRep)

# Perform PHATE on normalised bulk counts
# --------------------------------------------------------------------------

phate.out <- phate(Matrix::t(exprs_no_batch), ndim = 10) # Runs PHATE diffusion map
phate.out.embed <- data.frame(phate.out$embedding)

Group <- as.character(gsub( "\\_[0-9]*$", "", rownames(phate.out.embed)))
Names <- as.character(rownames(phate.out.embed))
phate.out.embed$Group <- factor(Group, levels = unique(Group))
phate.out.embed$Names <- Names

# Make plots of the first two principle components
ggplot(phate.out.embed, aes(x = phate.out.embed[,1], y = phate.out.embed[,2])) +
  stat_ellipse(data = phate.out.embed, aes(x = phate.out.embed[,1], y = phate.out.embed[,2], fill = Group), alpha = 0.3, geom = "polygon", type = "norm", level = 0.5, color = "black") +
  geom_point(aes(fill = Group), alpha = 0.8, size = 6, shape = 21, colour = "black") +
  geom_text(aes(label = Names), alpha = 0.8, size = 4, colour = "black") +
  scale_fill_manual("Sample type", values = c("#fdae61", "#f46d43", "#d53e4f", "#3288bd")) +
  scale_x_continuous(position = "top") +
  xlab("PHATE1") +
  scale_y_continuous(position = "left") +
  ylab("PHATE2") +
  theme_classic() +
  theme(plot.title = element_blank(), axis.text = element_text(color = "black"), panel.background = element_rect(fill = "#f0f0f0"), axis.line = element_blank(), aspect.ratio = 1) +
  ggsave("Pseudo-bulk_PHATE1_PHATE2_all.pdf")

# Identify clusters of DEGs based on expression
# --------------------------------------------------------------------------

# Scale gene expression values
zscores <- t(apply(exprs_no_batch, 1, scale))
colnames(zscores) <- colnames(allDGEList)
zscores <- zscores[! is.na(zscores[ ,1]), ]

# Subset to GOIs
zscores_GOI <- zscores[total_DGE, ]

# Perform heirarchical clustering and return gene cluster membership
hr <- hclust(dist(zscores_GOI))
mycl <- cutree(hr, k = 6)
clusterCols <- viridis(6)
myClusterSideBar <- clusterCols[mycl]

pdf(paste0("DGE_total_genes_HC_expression_heatmap.pdf"))
heatmap.2(zscores_GOI, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="row", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-2.5, 2.5, 1), col = magma(5, direction = -1), RowSideColors= myClusterSideBar)
dev.off()

# Identify clusters of DEGs based on expression
# --------------------------------------------------------------------------

zscores_GOI_annot <- data.frame(zscores_GOI)
idx <- match(rownames(zscores_GOI_annot), rownames(rowData_summed))
zscores_GOI_annot$Ensembl <- rowData_summed$Ensembl [idx]
zscores_GOI_annot$EntrezID <- rowData_summed$EntrezID [idx]
zscores_GOI_annot$GeneSymbol <- rowData_summed$GeneSymbol [idx]

idx <- match(rownames(zscores_GOI_annot), names(mycl))
zscores_GOI_annot$Cluster <- mycl [idx]

# Make boxplots of z-scores in each cluster for each sample
cluster1 <- zscores_GOI[zscores_GOI_annot$Cluster == 1, ]
cluster2 <- zscores_GOI[zscores_GOI_annot$Cluster == 2, ]
cluster3 <- zscores_GOI[zscores_GOI_annot$Cluster == 3, ]
cluster4 <- zscores_GOI[zscores_GOI_annot$Cluster == 4, ]
cluster5 <- zscores_GOI[zscores_GOI_annot$Cluster == 5, ]
cluster6 <- zscores_GOI[zscores_GOI_annot$Cluster == 6, ]

clusterList <- list("cluster1" = cluster1, "cluster2" = cluster2, "cluster3" = cluster3, "cluster4" = cluster4, "cluster5" = cluster5, "cluster6" = cluster6)
cluster_eigengenes <- lapply(X = names(clusterList), function(x){ return(svd(clusterList[[x]])$v) })
names(cluster_eigengenes) <- names(clusterList)

eigens <- c()
for( i in names(cluster_eigengenes)){
  test <- cluster_eigengenes[[i]][,1]
  eigens <- cbind(eigens, test)
}
colnames(eigens) <- names(cluster_eigengenes)
rownames(eigens) <- colnames(zscores_GOI)
write.csv(eigens, "HC_cluster_eigengenes.csv")

for(i in names(clusterList)) {
  melted <- gather(data.frame(clusterList[[i]]))
  # melted$Days <- c(rep(3, nrow(clusterList[[i]])*3), rep(14, nrow(clusterList[[i]])*3), rep(35, nrow(clusterList[[i]])*3))
  
  ggplot(data = melted, aes(x=key, y=value)) +
    geom_hline(yintercept = 0) +
    # geom_boxplot(aes(fill=factor(Days))) +
    geom_boxplot() +
    scale_fill_viridis_d(option = "inferno") +
    theme_classic() +
    ggsave(paste0("markers/HC_", i,"/Pseudo-bulk_DGE_HC_", i, "_boxplot.pdf"), useDingbats = FALSE)
}

# Identify genes most highly correlated with cluster eigengene

cor_list <- lapply(names(cluster_eigengenes), function(x){apply(zscores_GOI, 1, cor, cluster_eigengenes[[x]][,1])})
print(cor_list[[1]][order(-cor_list[[1]])][1:20])

# Make gene universe(s)
universe <- as.character(unique(rowData_summed[rowData_summed$Any_Active, "Ensembl"]))
universe_entrez <- as.character(unique(rowData_summed[rowData_summed$Any_Active, "EntrezID"]))
universe_genesymbol <- as.character(unique(rowData_summed[rowData_summed$Any_Active, "GeneSymbol"]))

# Makes hallmark geneset for enrichment testing
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_t2g <- h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()

# Makes hallmark geneset for enrichment testing
metabolic_pathways <- read.csv("~/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Metabolic_pathways_genes.csv", header = TRUE)
metabolic_pathways$Gene <- as.character(metabolic_pathways$Gene)
metabolic_pathways$Metabolic_pathway <- as.character(metabolic_pathways$Metabolic_pathway)

# Perform enrichment analysis on markers
# --------------------------------------------------------------------------
dir.create("markers")

for(i in names(clusterList)){
  dir.create(paste0("markers/HC_", i))
  interesting <- data.frame(clusterList[[i]])
  idx <- match(rownames(interesting), rownames(rowData_summed))
  interesting$Ensembl <- rowData_summed$Ensembl [idx]
  interesting$EntrezID <- rowData_summed$EntrezID [idx]
  interesting$GeneSymbol <- rowData_summed$GeneSymbol [idx]
  correlation <- apply(zscores_GOI[rownames(interesting),], 1, cor.test, cluster_eigengenes[[i]][,1])
  corr.df <- data.frame("correlation" = unlist(lapply(correlation, `[[`, "estimate")), row.names = names(correlation))
  corr.df$p.value <- unlist(lapply(correlation, `[[`, "p.value"))
  corr.df$p.adj <- p.adjust(corr.df$p.value, method = "bonf")
  
  idx <- match(rownames(interesting), rownames(corr.df))
  interesting$correlation <- corr.df$correlation [idx]
  interesting$p.value <- corr.df$p.value [idx]
  interesting$p.adj <- corr.df$p.adj [idx]
  interesting <- interesting[order(interesting$correlation), ]
  write.csv(interesting, paste0("markers/HC_", i,"/HC_", i, "_markers.csv"))
  
  # Make heatmap
  interesting_plot <- interesting[, ]
  interesting_plot <- interesting_plot[order(interesting_plot$correlation),]
  interesting_plot <- as.character(rownames(interesting_plot))
  zscores_GOI_HC <- zscores[interesting_plot, ]
  pdf(paste0("markers/HC_", i,"/HC_", i,"_expression_heatmap.pdf"))
  heatmap.2(zscores_GOI_HC, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="none", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-2.5, 2.5, 1), col = magma(5, direction = -1))
  dev.off()
  
  # All-regulated markers in each tissue
  # --------------------------------------------------------------------------
  
  # Perform GO enrichment analysis
  BPenrich <- enrichGO(interesting[interesting$p.adj < 0.05, "Ensembl"],
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
  write.csv(GOBP_GOI, paste0("markers/HC_", i,"/GOBP_markers_all_HC_", i, ".csv"))
  
  MFenrich <- enrichGO(interesting[interesting$p.adj < 0.05, "Ensembl"],
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
  write.csv(GOMF_GOI, paste0("markers/HC_", i,"/GOMF_markers_all_HC_", i, ".csv"))
  
  CCenrich <- enrichGO(interesting[interesting$p.adj < 0.05, "Ensembl"],
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
  write.csv(GOCC_GOI, paste0("markers/HC_", i,"/GOCC_markers_all_HC_", i, ".csv"))
  
  # Perform HALLMARK enrichment analysis
  HallmarkEnrich <- enricher(gene = interesting[interesting$p.adj < 0.05, "GeneSymbol"], 
                             TERM2GENE = h_t2g,
                             pvalueCutoff = 1,
                             qvalueCutoff = 1,
                             pAdjustMethod = "bonferroni",
                             minGSSize = 15,
                             maxGSSize = 500,
                             universe = universe_genesymbol)
  Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
  write.csv(Hallmark_GOI, paste0("markers/HC_", i,"/HALLMARK_markers_all_HC_", i, ".csv"))

  # Perform Metabolic_pathways enrichment analysis
  Metabolic_pathwaysEnrich <- enricher(gene = interesting[interesting$p.adj < 0.05, "GeneSymbol"], 
                                       TERM2GENE = metabolic_pathways,
                                       pvalueCutoff = 1,
                                       qvalueCutoff = 1,
                                       pAdjustMethod = "bonferroni",
                                       minGSSize = 10,
                                       maxGSSize = 500,
                                       universe = universe_genesymbol)
  Metabolic_pathways_GOI <- as.data.frame(Metabolic_pathwaysEnrich@result)
  write.csv(Metabolic_pathways_GOI, paste0("markers/HC_", i,"/Metabolic_pathways_all_HC_", i, ".csv"))
    
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[interesting$p.adj < 0.05, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/HC_", i,"/KEGG_markers_all_HC_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[interesting$p.adj < 0.05, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/HC_", i,"/Reactome_markers_all_HC_", i, ".csv"))
}
