#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Make tissue specific SVD signatures
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
library('clues')

set.seed(100)

# Load Pseudo-bulk counts and row data SingleCellExperiment
exprs <- readRDS("Pseudo-bulk_whole_experiment_all.rds") # uses whole dataset if wolfpack
rowData_summed <- read.csv("Pseudo-bulk_rowData_all.csv", header = TRUE, row.names = 1)

# Load DGE data
# --------------------------------------------------------------------------
# Make list of differentially expressed genes
liver_primary <- read.csv("Liver_Primary_DEG_0LFC.csv", header = TRUE)
liver_primary <- liver_primary[liver_primary$adj.P.Val < 0.05, ]
liver_primary <- liver_primary[order(liver_primary$t),]
liver_primary <- as.character(liver_primary[, "X"])

lung_primary <- read.csv("Lung_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary <- lung_primary[lung_primary$adj.P.Val < 0.05, ]
lung_primary <- lung_primary[order(lung_primary$t),]
lung_primary <- as.character(lung_primary[, "X"])

liver_lung <- read.csv("Liver_Lung_DEG_0LFC.csv", header = TRUE)
liver_lung <- liver_lung[liver_lung$adj.P.Val < 0.05, ]
liver_lung <- liver_lung[order(liver_lung$t),]
liver_lung <- as.character(liver_lung[, "X"])

total_DGE <- as.character(unique(c(liver_primary, lung_primary, liver_lung)))
total_Liver <- as.character(unique(c(liver_primary, liver_lung)))
total_Lung <- as.character(unique(c(lung_primary, liver_lung)))
total_Primary <- as.character(unique(c(liver_primary, lung_primary)))

DGE_list <- list("Total_Liver" = total_Liver, "Total_Lung" = total_Lung, "Total_Primary" = total_Primary)

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

# Scale gene expression values
zscores <- t(apply(exprs_no_batch, 1, scale))
colnames(zscores) <- colnames(allDGEList)
zscores <- zscores[! is.na(zscores[ ,1]), ]

# Identify tissue specific signature based on DEG expression
# --------------------------------------------------------------------------

for(y in names(DGE_list)){
  # Subset to GOIs
  zscores_GOI <- zscores[DGE_list[[y]], ]
  
  # Perform heirarchical clustering and return gene cluster membership
  hr <- hclust(dist(zscores_GOI))
  # itissil <- clues(zscores_GOI, n0 = 6, strengthMethod = "sil", itmax = 100)
  # mycl <- cutree(hr, k = itissil$g)
  mycl <- cutree(hr, k = 2)
  clusterCols <- viridis(6)
  myClusterSideBar <- clusterCols[mycl]
  
  pdf(paste0("markers/HC_", y,"/HC_", y,"_expression_heatmap.pdf"))
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

  clusterList <- lapply(X = unique(zscores_GOI_annot$Cluster), FUN = function(x){
    return(zscores_GOI[zscores_GOI_annot$Cluster == x, ])
  })
  
  names(clusterList) <- paste0("cluster", 1:length(clusterList))
  
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
      ggsave(paste0("markers/HC_", y,"/Pseudo-bulk_DGE_HC_", i, "_boxplot.pdf"), useDingbats = FALSE)
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
    write.csv(interesting, paste0("markers/HC_", y,"/HC_", i, "_markers.csv"))
    
    # Make heatmap
    interesting_plot <- interesting[, ]
    interesting_plot <- interesting_plot[order(interesting_plot$correlation),]
    interesting_plot <- as.character(rownames(interesting_plot))
    zscores_GOI_HC <- zscores_GOI[interesting_plot, ]
    pdf(paste0("markers/HC_", y,"/HC_", i,"_expression_heatmap.pdf"))
    heatmap.2(zscores_GOI_HC, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="none", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-2.5, 2.5, 1), col = magma(5, direction = -1))
    dev.off()
    
    # All-regulated markers in each tissue
    # --------------------------------------------------------------------------
    
    # Perform GO enrichment analysis
    BPenrich <- enrichGO(interesting[, "Ensembl"],
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
    write.csv(GOBP_GOI, paste0("markers/HC_", y,"/GOBP_markers_all_HC_", i, ".csv"))
    
    MFenrich <- enrichGO(interesting[, "Ensembl"],
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
    write.csv(GOMF_GOI, paste0("markers/HC_", y,"/GOMF_markers_all_HC_", i, ".csv"))
    
    CCenrich <- enrichGO(interesting[, "Ensembl"],
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
    write.csv(GOCC_GOI, paste0("markers/HC_", y,"/GOCC_markers_all_HC_", i, ".csv"))
    
    # Perform HALLMARK enrichment analysis
    HallmarkEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                               TERM2GENE = h_t2g,
                               pvalueCutoff = 1,
                               qvalueCutoff = 1,
                               pAdjustMethod = "bonferroni",
                               minGSSize = 15,
                               maxGSSize = 500,
                               universe = universe_genesymbol)
    Hallmark_GOI <- as.data.frame(HallmarkEnrich@result)
    write.csv(Hallmark_GOI, paste0("markers/HC_", y,"/HALLMARK_markers_all_HC_", i, ".csv"))
    
    # Perform Metabolic_pathways enrichment analysis
    Metabolic_pathwaysEnrich <- enricher(gene = interesting[, "GeneSymbol"], 
                                         TERM2GENE = metabolic_pathways,
                                         pvalueCutoff = 1,
                                         qvalueCutoff = 1,
                                         pAdjustMethod = "bonferroni",
                                         minGSSize = 10,
                                         maxGSSize = 500,
                                         universe = universe_genesymbol)
    Metabolic_pathways_GOI <- as.data.frame(Metabolic_pathwaysEnrich@result)
    write.csv(Metabolic_pathways_GOI, paste0("markers/HC_", y,"/Metabolic_pathways_all_HC_", i, ".csv"))
    
    # Perform KEGG enrichment analysis
    KEGGenrichsig <- enrichKEGG(interesting[, "EntrezID"],
                                organism = "hsa",
                                keyType = "kegg",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "bonferroni",
                                universe = universe_entrez)
    write.csv(KEGGenrichsig, paste0("markers/HC_", y,"/KEGG_markers_all_HC_", i, ".csv"))
    
    # Perform REACTOME enrichment analysis
    REACTOMEenrichsig <- enrichPathway(interesting[, "EntrezID"],
                                       organism = "human",
                                       pvalueCutoff = 0.05,
                                       pAdjustMethod = "bonferroni",
                                       readable = TRUE,
                                       universe = universe_entrez)
    write.csv(REACTOMEenrichsig, paste0("markers/HC_", y,"/Reactome_markers_all_HC_", i, ".csv"))
  }
}
