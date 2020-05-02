#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Identify differentially expressed genes and differential cell composition between each tissue
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
library('msigdbr')
library('org.Hs.eg.db')
library('fgsea')
library('data.table')
library('sva')

set.seed(100)

# Load Pseudo-bulk counts and row data SingleCellExperiment
exprs <- readRDS("Pseudo-bulk_whole_experiment_all.rds") # uses whole dataset if wolfpack
rowData_summed <- read.csv("Pseudo-bulk_rowData_all.csv", header = TRUE, row.names = 1)

# Make lists for enrichment testing
# --------------------------------------------------------------------------

# Make mdbsig lists for CAMERA analysis
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

metabolic_pathways <- read.csv("~/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/all_data/pseudo-bulk_DGE/Metabolic_pathways_genes.csv", header = TRUE)
metabolic_pathways$Gene <- as.character(metabolic_pathways$Gene)
metabolic_pathways$Metabolic_pathway <- as.character(metabolic_pathways$Metabolic_pathway)
metabolic_list <- metabolic_pathways %>% split(x = .$Gene, f = .$Metabolic_pathway)

mdbsig <- list("Hallmark" = h_list, "PID_pathways" = pid_list, "Chem_genetic_interventions" = cgp_list, "Transcription_factors" = tft_list, "Cancer_neighbourhoods" = CGN_list, "GO_BP" = gobp_list, "Oncogenic_signatures" = c6_list, "Metabolic_pathways" = metabolic_list)
dir.create("Primary-Mets")
dir.create("Primary-Mets/msigdb")

# Find DEGs between mets and primary
# --------------------------------------------------------------------------

# Make sample info table
sampleName <- colnames(exprs)
sampleType <- factor(c(gsub( "\\_[0-9]*$", "", sampleName)))
sampleRep <- factor(c(gsub( "^.*?_","", sampleName)))
sampleSite <- factor(c(rep("Metastasis", 12), rep("Primary", 4)))
sampleTable <- data.frame(sampleName, sampleType, sampleRep, sampleSite)

# Find genes detected in all replicates of a given tissue type
active_liver <- rowSums(exprs[,1:4] >= 1) == 4
active_LN <- rowSums(exprs[,5:8] >= 1) == 4
active_lung <- rowSums(exprs[,9:12] >= 1) == 4
active_primary <- rowSums(exprs[,13:16] >= 1) == 4
geneActivity <- data.frame("Liver" = active_liver, "LN" = active_LN, "Lung" = active_lung, "Primary" = active_primary)

# Subset expression data to those genes actively expressed in either sample type
expression <- exprs[geneActivity[ , "Primary"] | (geneActivity[ , "Liver"] & geneActivity[ , "Lung"] & geneActivity[ , "LN"]), ]
genes <- data.frame(rowData_summed[is.element(rownames(rowData_summed), rownames(expression)), ])

# Make DGEList object and normalise using TMM
allDGEList <- DGEList(counts = expression, group = sampleSite, genes = genes)
allDGEList <- calcNormFactors(allDGEList, method = "TMM")

# Make design matrix, fit linear models, list comparisons
type <- factor(sampleTable$sampleSite, levels = c("Primary", "Metastasis"))
rep <- factor(sampleTable$sampleRep, levels = c(1:4))
design <- model.matrix(~0 + rep + type)
rownames(design) <- sampleName
allDGEList <- voom(allDGEList, design, plot = TRUE)

# Identify significantly altered mSigDB genesets
for(sig in names(mdbsig)) {
  camera.out <- camera(allDGEList$E, mdbsig[[sig]], design, contrast = 5, inter.gene.cor=0.01, trend.var = TRUE)
  write.csv(camera.out, paste0("Primary-Mets/msigdb/Primary-Mets_", sig, "_camera.csv"))
}

# Identify differentially expressed genes at 0 LFC threshold.
lbFit <- lmFit(allDGEList, design)
lbFit3 <- treat(lbFit, lfc = 0)
DEG3 <- decideTests(lbFit3)
print(summary(DEG3))
results <- data.frame(topTreat(lbFit3, coef = 5, n = Inf, p.value = Inf))
write.csv(results, "Primary-Mets/Primary-Mets_DEG_0LFC.csv")

# Perform fgsea enrichment analysis across mSigDB genesets
interesting <- results
interesting$GeneSymbol <- as.character(interesting$GeneSymbol)
interesting_in <- interesting[order(interesting$t),]
interesting_in <- interesting_in[! is.na(interesting_in$GeneSymbol), ]
interesting_in_vec <- interesting_in$t
names(interesting_in_vec) <- interesting_in$GeneSymbol

for(sig in names(mdbsig)) {
  fgseaRes <- fgsea(pathways = mdbsig[[sig]], 
                    stats = interesting_in_vec,
                    minSize=10,
                    maxSize=500,
                    nperm=10000)
  fgseaRes <- fgseaRes[order(pval), ]
  fwrite(fgseaRes, file = paste0("Primary-Mets/msigdb/Primary-Mets_", sig, "_fgsea.csv"), sep=",", sep2=c("", " ", ""))
}
