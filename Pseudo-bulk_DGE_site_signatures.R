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

for(i in names(DGE_list)){
  zscores_GOI <- zscores[DGE_list[[i]], ]
  pdf(paste0("markers/HC_", i,"/HC_", i,"_expression_heatmap.pdf"))
  heatmap.2(zscores_GOI, Rowv=TRUE, Colv="none", offsetRow = 0.01, labCol = colnames(zscores_GOI), labRow = FALSE, dendrogram="none", symm=FALSE, trace = "none", density.info = "none", breaks = seq(-2.5, 2.5, 1), col = magma(5, direction = -1))
  dev.off()
}

tissue_signatures <- lapply(X = names(DGE_list), function(x){
  # Subset to GOIs
  zscores_GOI <- zscores[DGE_list[[x]], ]
  
  # Perform heirarchical clustering and return gene cluster membership
  hr <- hclust(dist(zscores_GOI))
  
  # Subset to GOIs, transpose matrix and create lists for each tissue
  t_zscores_GOI <- t(zscores_GOI)
  
  Tissue_Sig <- t_zscores_GOI[sampleType == gsub( "^.*?_","", x), ]
  
  # Make signature for each tissue based on SVD of DEG expression in each sample
  signature_tmp <- data.frame("Signature" = colMedians(Tissue_Sig))
  signature_tmp$Gene_name <- colnames(Tissue_Sig)
  signature_tmp <- signature_tmp[order(signature_tmp$Signature), ]
  signature_tmp$Gene_name <- factor(signature_tmp$Gene_name, levels = signature_tmp$Gene_name)
  return(signature_tmp)
})
names(tissue_signatures) <- paste0(names(DGE_list), "_sig")
  
# Plot each signature
for(i in names(tissue_signatures)) {
  tissue_signatures[[i]] %>%
    ggplot(aes(x = Gene_name, y = Signature)) +
    geom_point(aes(color = Signature), alpha = 0.3) +
    scale_color_gradientn(colors = hcl.colors(n = 7, palette = "Blue-Red 2")) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    coord_flip() +
    ggsave(paste0("Pseudo-bulk_tissue_specific_", i, "_dotplot.pdf"), useDingbats = FALSE)
}

# Save R object containing signatures
saveRDS(tissue_signatures, "Pseudo-bulk_tissue_specific_signature.rds")
