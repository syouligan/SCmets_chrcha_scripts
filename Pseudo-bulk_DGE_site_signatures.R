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

# Idenitfy clusters of DEGs based on expression
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

# Idenitfy tissue specific signature based on DEG expression
# --------------------------------------------------------------------------

# Subset to GOIs, transpose matrix and create lists for each tissue
t_zscores_GOI <- t(zscores_GOI)

Liver_Sig <- t_zscores_GOI[sampleType == "Liver", ]
Lung_Sig <- t_zscores_GOI[sampleType == "Lung", ]
LN_Sig <- t_zscores_GOI[sampleType == "LN", ]
Primary_Sig <- t_zscores_GOI[sampleType == "Primary", ]

# Make signature for each tissue based on SVD of DEG expression in each sample
tissueList <- list("Liver_Sig" = Liver_Sig, "Lung_Sig" = Lung_Sig, "LN_Sig" = LN_Sig, "Primary_Sig" = Primary_Sig)
tissue_signatures <- lapply(X = names(tissueList), function(x){
  signature_tmp <- data.frame("Signature" = svd(tissueList[[x]])$v[,1])
  signature_tmp$Gene_name <- colnames(tissueList[[x]])
  signature_tmp$HC_Cluster <- factor(mycl[signature_tmp$Gene_name == names(mycl)])
  signature_tmp <- signature_tmp[order(signature_tmp$Signature), , drop = FALSE]
  signature_tmp$Gene_name <- factor(signature_tmp$Gene_name, levels = signature_tmp$Gene_name)
  return(signature_tmp)
})
names(tissue_signatures) <- names(tissueList)

# Plot each signature
for(i in names(tissue_signatures)) {
  tissue_signatures[[i]] %>%
    ggplot(aes(x = Gene_name, y = Signature)) +
    geom_point(aes(color = HC_Cluster)) +
    scale_colour_viridis_d() +
    ggsave(paste0("Pseudo-bulk_tissue_specific_", i, "_dotplot.pdf"), useDingbats = FALSE)
}

# Save R object containing signatures
saveRDS(tissue_signatures, "Pseudo-bulk_tissue_specific_signature.csv")

# Find processes in top and bottom markers
# --------------------------------------------------------------------------

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

for(i in names(tissue_signatures)){
  dir.create(paste0("markers/Sig300up_", i))
  interesting <- data.frame(tissue_signatures[[i]])[1:300, ]
  idx <- match(interesting$Gene_name, rownames(rowData_summed))
  interesting$Ensembl <- rowData_summed$Ensembl [idx]
  interesting$EntrezID <- rowData_summed$EntrezID [idx]
  interesting$GeneSymbol <- rowData_summed$GeneSymbol [idx]
  write.csv(interesting, paste0("markers/Sig300up_", i,"/Sig300up_", i, "_markers.csv"))
  
  # Make heatmap
  interesting_plot <- interesting[, ]
  interesting_plot <- interesting_plot[order(interesting_plot$Signature),]
  interesting_plot <- as.character(interesting_plot$Gene_name)
  zscores_GOI_HC <- zscores[interesting_plot, ]
  pdf(paste0("markers/Sig300up_", i,"/Sig300up_", i,"_expression_heatmap.pdf"))
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
  write.csv(GOBP_GOI, paste0("markers/Sig300up_", i,"/GOBP_markers_all_Sig300up_", i, ".csv"))
  
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
  write.csv(GOMF_GOI, paste0("markers/Sig300up_", i,"/GOMF_markers_all_Sig300up_", i, ".csv"))
  
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
  write.csv(GOCC_GOI, paste0("markers/Sig300up_", i,"/GOCC_markers_all_Sig300up_", i, ".csv"))
  
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
  write.csv(Hallmark_GOI, paste0("markers/Sig300up_", i,"/HALLMARK_markers_all_Sig300up_", i, ".csv"))
  
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
  write.csv(Metabolic_pathways_GOI, paste0("markers/Sig300up_", i,"/Metabolic_pathways_all_Sig300up_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/Sig300up_", i,"/KEGG_markers_all_Sig300up_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/Sig300up_", i,"/Reactome_markers_all_Sig300up_", i, ".csv"))
}

for(i in names(tissue_signatures)){
  dir.create(paste0("markers/Sig300down_", i))
  interesting <- data.frame(tissue_signatures[[i]])[(nrow(tissue_signatures[[i]]) - 300):nrow(tissue_signatures[[i]]), ]
  idx <- match(interesting$Gene_name, rownames(rowData_summed))
  interesting$Ensembl <- rowData_summed$Ensembl [idx]
  interesting$EntrezID <- rowData_summed$EntrezID [idx]
  interesting$GeneSymbol <- rowData_summed$GeneSymbol [idx]
  write.csv(interesting, paste0("markers/Sig300down_", i,"/Sig300down_", i, "_markers.csv"))
  
  # Make heatmap
  interesting_plot <- interesting[, ]
  interesting_plot <- interesting_plot[order(interesting_plot$Signature),]
  interesting_plot <- as.character(interesting_plot$Gene_name)
  zscores_GOI_HC <- zscores[interesting_plot, ]
  pdf(paste0("markers/Sig300down_", i,"/Sig300down_", i,"_expression_heatmap.pdf"))
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
  write.csv(GOBP_GOI, paste0("markers/Sig300down_", i,"/GOBP_markers_all_Sig300down_", i, ".csv"))
  
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
  write.csv(GOMF_GOI, paste0("markers/Sig300down_", i,"/GOMF_markers_all_Sig300down_", i, ".csv"))
  
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
  write.csv(GOCC_GOI, paste0("markers/Sig300down_", i,"/GOCC_markers_all_Sig300down_", i, ".csv"))
  
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
  write.csv(Hallmark_GOI, paste0("markers/Sig300down_", i,"/HALLMARK_markers_all_Sig300down_", i, ".csv"))
  
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
  write.csv(Metabolic_pathways_GOI, paste0("markers/Sig300down_", i,"/Metabolic_pathways_all_Sig300down_", i, ".csv"))
  
  # Perform KEGG enrichment analysis
  KEGGenrichsig <- enrichKEGG(interesting[, "EntrezID"],
                              organism = "hsa",
                              keyType = "kegg",
                              pvalueCutoff = 0.05,
                              pAdjustMethod = "bonferroni",
                              universe = universe_entrez)
  write.csv(KEGGenrichsig, paste0("markers/Sig300down_", i,"/KEGG_markers_all_Sig300down_", i, ".csv"))
  
  # Perform REACTOME enrichment analysis
  REACTOMEenrichsig <- enrichPathway(interesting[, "EntrezID"],
                                     organism = "human",
                                     pvalueCutoff = 0.05,
                                     pAdjustMethod = "bonferroni",
                                     readable = TRUE,
                                     universe = universe_entrez)
  write.csv(REACTOMEenrichsig, paste0("markers/Sig300down_", i,"/Reactome_markers_all_Sig300down_", i, ".csv"))
}

