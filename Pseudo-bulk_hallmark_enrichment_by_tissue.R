#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Identify differentially hallmark genes sets between each tissue
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data")
  place <- "wolfpack"
}

# Libraries
library('dplyr')
library('tidyr')
library('scater')
library('scran')
library('ggplot2')
library('Matrix')
library('SAVER')
library('limma')
library('edgeR')
library('gtools')
library('msigdbr')

set.seed(100)

# Load prefiltered SingleCellExperiment
if(place == "local") {
  filtered_exp <- readRDS("Prefiltered_QC_experiment_practice.rds") # uses practice data if local
} else {
  filtered_exp <- readRDS("Prefiltered_QC_experiment_all.rds") # uses whole dataset if wolfpack
}

# filtered_exp <- readRDS("~/cloudstor/sarah_projects/DEP_hg19_SCmets_chrcha/project_results/prefiltered/practice_all_data/redundant_incomplete/Prefiltered_experiment_Practice_merge_cluster_wCC.rds")

# Sum counts across cells from each tissue within each replicate
summed <- aggregateAcrossCells(filtered_exp, ids=DataFrame(label=filtered_exp$Tissue, sample=filtered_exp$Replicate))
exprs <- counts(summed)
colnames(exprs) <- unique(paste0(filtered_exp$Tissue, "_", filtered_exp$Replicate))
exprs <- exprs[ ,order(colnames(exprs))]

# Find genes detected in all replicates of a given tissue type
active_liver <- rowSums(exprs[,1:4] >= 1) == 4
active_LN <- rowSums(exprs[,5:8] >= 1) == 4
active_lung <- rowSums(exprs[,9:12] >= 1) == 4
active_primary <- rowSums(exprs[,13:16] >= 1) == 4
geneActivity <- data.frame("Liver" = active_liver, "LN" = active_LN, "Lung" = active_lung, "Primary" = active_primary)

# Identify differentially expressed genes in MDA cells from different tissues
# --------------------------------------------------------------------------

# Make sample info table
sampleName <- colnames(exprs)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleTable <- data.frame(sampleName, sampleType, sampleRep)

# Make table of comparisons across tissue
comparisons <- combinations(4, 2, unique(as.character(sampleType)))
colnames(comparisons) <- c("Tissue_1", "Tissue_2")
rownames(comparisons) <- 1:nrow(comparisons)

h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_list <- h_df %>% split(x = .$human_gene_symbol, f = .$gs_name)

cgp_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
cgp_list <- cgp_df %>% split(x = .$gene_symbol, f = .$gs_name)

pid_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID")
pid_list <- pid_df %>% split(x = .$gene_symbol, f = .$gs_name)

tft_df <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT")
tft_list <- tft_df %>% split(x = .$gene_symbol, f = .$gs_name)

CGN_df <- msigdbr(species = "Homo sapiens", category = "C4", subcategory = "CGN")
CGN_list <- CGN_df %>% split(x = .$gene_symbol, f = .$gs_name)

gobp_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gobp_list <- gobp_df %>% split(x = .$gene_symbol, f = .$gs_name)

c6_df <- msigdbr(species = "Homo sapiens", category = "C6")
c6_list <- c6_df %>% split(x = .$gene_symbol, f = .$gs_name)

mdbsig <- list("Hallmark" = h_list, "PID_pathways" = pid_list, "Chem_genetic_interventions" = cgp_list, "Transcription_factors" = tft_list, "Cancer_neighbourhoods" = CGN_list, "GO_BP" = gobp_list, "Oncogenic_signatures" = c6_list)

dir.create("pseudo-bulk_DGE/msigdb")
# Find DEGs for each comparison
for(i in 1:nrow(comparisons)) {
  
  # Make vectors for each tissue in the comparison
  tissue1 <- comparisons[i, 1]
  tissue2 <- comparisons[i, 2]
  
  # Subset expression data to those genes actively expressed in either sample type
  expression <- exprs[geneActivity[ ,tissue1] | geneActivity[ ,tissue2], c(sampleType == tissue1 | sampleType == tissue2)]
  
  # Make DGEList object and normalise using TMM
  allDGEList <- DGEList(counts = expression, group = sampleType[sampleType == tissue1 | sampleType == tissue2])
  allDGEList <- calcNormFactors(allDGEList, method = "TMM")
  
  # Make design matrix, fit linear models, list comparisons
  type <- factor(sampleTable$sampleType, levels = c(tissue1, tissue2))
  rep <- factor(sampleTable$sampleRep, levels = c(1:4))
  design <- model.matrix(~0 + rep + type)
  rownames(design) <- sampleName[sampleType == tissue1 | sampleType == tissue2]
  allDGEList <- voom(allDGEList, design, plot = FALSE)
  
  # Test gene sets
  for(sig in names(mdbsig)) {
    camera.out <- camera(allDGEList$E, mdbsig[[sig]], design, inter.gene.cor=0.01, trend.var = TRUE)
    write.csv(camera.out, paste0("pseudo-bulk_DGE/msigdb/", tissue1, "_", tissue2, "_", sig, "_camera.csv"))
}}  
