#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Make volcano plots of top DEGs between tissues (pseudobulk)
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
library('ggplot2')
library('ggrepel')
library('readr')
library('Matrix')
library('clusterProfiler')
library('ReactomePA')
library('org.Hs.eg.db')
library('gplots')
library('msigdbr')

# Make list of differentially expressed genes
liver_primary <- read.csv("Liver_Primary_DEG_0LFC.csv", header = TRUE)
liver_primary$MinusLog10 <- -log10(liver_primary$adj.P.Val)
liver_primary$Primary <- liver_primary$logFC > 0

ln_primary <- read.csv("LN_Primary_DEG_0LFC.csv", header = TRUE)
ln_primary$MinusLog10 <- -log10(ln_primary$adj.P.Val)
ln_primary$Primary <- ln_primary$logFC > 0

lung_primary <- read.csv("Lung_Primary_DEG_0LFC.csv", header = TRUE)
lung_primary$MinusLog10 <- -log10(lung_primary$adj.P.Val)
lung_primary$Primary <- lung_primary$logFC > 0

ggplot(liver_primary) +
  geom_point(aes(x = logFC, y = MinusLog10, size = abs(t)), color = "grey80", alpha = 0.5) +
  geom_point(data = liver_primary[1:20,], aes(x = logFC, y = MinusLog10, color = Primary, size = abs(t)), alpha = 0.5) +
  geom_text_repel(data = liver_primary[1:20,], aes(x = logFC, y = MinusLog10, label = GeneSymbol, color = Primary), point.padding = 0.2) +
  scale_colour_manual(values = c("blue", "red")) +
  theme_classic()

ggplot(ln_primary) +
  geom_point(aes(x = logFC, y = MinusLog10, size = abs(t)), color = "grey80", alpha = 0.5) +
  geom_point(data = ln_primary[1:20,], aes(x = logFC, y = MinusLog10, color = Primary, size = abs(t)), alpha = 0.5) +
  geom_text_repel(data = ln_primary[1:20,], aes(x = logFC, y = MinusLog10, label = GeneSymbol, color = Primary), point.padding = 0.2) +
  scale_colour_manual(values = c("blue", "red")) +
  theme_classic()

ggplot(lung_primary) +
  geom_point(aes(x = logFC, y = MinusLog10, size = abs(t)), color = "grey80", alpha = 0.5) +
  geom_point(data = lung_primary[1:20,], aes(x = logFC, y = MinusLog10, color = Primary, size = abs(t)), alpha = 0.5) +
  geom_text_repel(data = lung_primary[1:20,], aes(x = logFC, y = MinusLog10, label = GeneSymbol, color = Primary), point.padding = 0.2) +
  scale_colour_manual(values = c("blue", "red")) +
  theme_classic()

