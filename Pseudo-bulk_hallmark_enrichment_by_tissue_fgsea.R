#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Identify differentially regulated msigdb genes sets between each tissue
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  setwd("/Users/mac/cloudstor/sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/practice_all_data") # Uses practice data (5% of cells from each sample) if running locally
  place <- "local"
} else {
  setwd("/share/ScratchGeneral/scoyou/sarah_projects/SCMDA231mets_chrcha/prefiltering/all_data")
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
library('fgsea')
library('org.Hs.eg.db')
library('data.table')

# Make list of differentially expressed genes
liver_primary <- read.csv("pseudo-bulk_DGE/Liver_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary <- read.csv("pseudo-bulk_DGE/LN_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary <- read.csv("pseudo-bulk_DGE/Lung_Primary_DEG_0LFC.csv", header = TRUE)
liver_lung <- read.csv("pseudo-bulk_DGE/Liver_Lung_DEG_0LFC.csv", header = TRUE)
liver_ln <- read.csv("pseudo-bulk_DGE/Liver_LN_DEG_0LFC.csv", header = TRUE)
ln_lung <- read.csv("pseudo-bulk_DGE/LN_Lung_DEG_0LFC.csv", header = TRUE)

DGE_list <- list("Liver" = liver_primary, "LN" = ln_primary, "Lung" = lung_primary, "Liver_Lung" = liver_lung, "Liver_LN" = liver_ln, "LN_Lung" = ln_lung)

# Make lists of mSigDB datasets for fgsea enrichment
h_df <- msigdbr(species = "Homo sapiens", category = "H")
h_list <- h_df %>% split(x = .$gene_symbol, f = .$gs_name)
h_list.char <- lapply(names(h_list), 
                    function(x) {
                      b <- as.character(h_list[[x]])
                      b})
names(h_list.char) <- names(h_list)

cgp_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CGP")
cgp_list <- cgp_df %>% split(x = .$gene_symbol, f = .$gs_name)
cgp_list.char <- lapply(names(cgp_list), 
                        function(x) {
                          b <- as.character(cgp_list[[x]])
                          b})
names(cgp_list.char) <- names(cgp_list)

pid_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:PID")
pid_list <- pid_df %>% split(x = .$gene_symbol, f = .$gs_name)
pid_list.char <- lapply(names(pid_list), 
                        function(x) {
                          b <- as.character(pid_list[[x]])
                          b})
names(pid_list.char) <- names(pid_list)


tft_df <- msigdbr(species = "Homo sapiens", category = "C3", subcategory = "TFT")
tft_list <- tft_df %>% split(x = .$gene_symbol, f = .$gs_name)
tft_list.char <- lapply(names(tft_list), 
                        function(x) {
                          b <- as.character(tft_list[[x]])
                          b})
names(tft_list.char) <- names(tft_list)

CGN_df <- msigdbr(species = "Homo sapiens", category = "C4", subcategory = "CGN")
CGN_list <- CGN_df %>% split(x = .$gene_symbol, f = .$gs_name)
CGN_list.char <- lapply(names(CGN_list), 
                        function(x) {
                          b <- as.character(CGN_list[[x]])
                          b})
names(CGN_list.char) <- names(CGN_list)

gobp_df <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
gobp_list <- gobp_df %>% split(x = .$gene_symbol, f = .$gs_name)
gobp_list.char <- lapply(names(gobp_list), 
                        function(x) {
                          b <- as.character(gobp_list[[x]])
                          b})
names(gobp_list.char) <- names(gobp_list)

c6_df <- msigdbr(species = "Homo sapiens", category = "C6")
c6_list <- c6_df %>% split(x = .$gene_symbol, f = .$gs_name)
c6_list.char <- lapply(names(c6_list), 
                        function(x) {
                          b <- as.character(c6_list[[x]])
                          b})
names(c6_list.char) <- names(c6_list)

mdbsig <- list("Hallmark" = h_list.char, "PID_pathways" = pid_list.char, "Chem_genetic_interventions" = cgp_list.char, "Transcription_factors" = tft_list.char, "Cancer_neighbourhoods" = CGN_list.char, "GO_BP" = gobp_list.char, "Oncogenic_signatures" = c6_list.char)

dir.create("pseudo-bulk_DGE/msigdb")

# Perform fgsea enrichment analysis across hallmark gene sets
for(i in names(DGE_list)) {
  interesting <- DGE_list[[i]]
  interesting$GeneSymbol <- as.character(interesting$GeneSymbol)
  interesting_in <- interesting[order(interesting$t),]
  interesting_in <- interesting_in[! is.na(interesting_in$GeneSymbol), ]
  interesting_in_vec <- interesting_in$t
  names(interesting_in_vec) <- interesting_in$GeneSymbol

 for(sig in names(mdbsig)) {
  fgseaRes <- fgsea(pathways = mdbsig[[sig]], 
                  stats = interesting_in_vec,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)
  fgseaRes <- fgseaRes[order(pval), ]
  fwrite(fgseaRes, file=paste0("pseudo-bulk_DGE/msigdb/", i, "_", sig, "_fgsea.csv"), sep=",", sep2=c("", " ", ""))
}
}

