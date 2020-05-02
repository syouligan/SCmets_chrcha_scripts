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

# Find genes detected in all replicates of a given tissue type
active_liver <- rowSums(exprs[,1:4] >= 1) == 4
active_LN <- rowSums(exprs[,5:8] >= 1) == 4
active_lung <- rowSums(exprs[,9:12] >= 1) == 4
active_primary <- rowSums(exprs[,13:16] >= 1) == 4
geneActivity <- data.frame("Liver" = active_liver, "LN" = active_LN, "Lung" = active_lung, "Primary" = active_primary)

# Plot samples based on 2 PCs of pseudo-bulk expression measurements
# --------------------------------------------------------------------------

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

# Perform principal component analysis
zscore.pca <- prcomp(t(exprs_no_batch), scale = TRUE)

# Summary of Principle Components
PCA_summary <- summary(zscore.pca)
PCA_stats <- data.frame(t(PCA_summary$importance))
PCA_stats$PC <- factor(as.character(rownames(PCA_stats)), levels = rownames(PCA_stats))
colnames(PCA_stats) <- c("StDev", "Proportion", "Cumulative", "PC")

# Label with sample and replicate information
PCA_df <- as.data.frame(zscore.pca$x)
Group <- as.character(gsub( "\\_[0-9]*$", "", rownames(PCA_df)))
Names <- as.character(rownames(PCA_df))
PCA_df$Group <- factor(Group, levels = unique(Group))
PCA_df$Names <- Names

# Make axis labels with proportion of variance explained
variance_exp  <- paste0(rownames(PCA_stats), " (", round(PCA_stats$Proportion, 2),")")

# Plot proportion of variance explained and cumulative proportion of variance explained
ggplot(PCA_stats) +
  geom_line(aes(x = PC, y = Proportion, group = 1), size = 0.5, colour = "grey", linetype = "twodash") +
  geom_line(aes(x = PC, y = Cumulative, group = 1), size = 0.5, colour = "black") +
  geom_point(aes(x = PC, y = Proportion), size = 3, colour = "grey") +
  geom_point(aes(x = PC, y = Cumulative), size = 3, colour = "black") +
  geom_hline(aes(yintercept=0.5), color="#525252", linetype="dashed") +
  labs(x = "Principle component", y = "Variance explained (proportion)") +
  theme_classic() +
  theme(plot.title = element_blank(), axis.text = element_text(color = "black"), axis.line = element_line(color = "black"), axis.ticks = element_line(color = "black"), aspect.ratio = 0.2) +
  ggsave("Pseudo-bulk_scree_all.pdf", useDingbats = FALSE)

# Make plots of the first two principle components
ggplot(PCA_df, aes(x = PCA_df[,1], y = PCA_df[,2])) +
  stat_ellipse(data = PCA_df, aes(x = PCA_df[,1], y = PCA_df[,2], fill = Group), alpha = 0.3, geom = "polygon", type = "norm", level = 0.5, color = "black") +
  geom_point(aes(fill = Group), alpha = 0.8, size = 6, shape = 21, colour = "black") +
  geom_text(aes(label = Names), alpha = 0.8, size = 4, colour = "black") +
  scale_fill_manual("Sample type", values = c("#fdae61", "#f46d43", "#d53e4f", "#3288bd")) +
  scale_x_continuous(position = "top") +
  xlab(variance_exp[1]) +
  scale_y_continuous(position = "left") +
  ylab(variance_exp[2]) +
  theme_classic() +
  theme(plot.title = element_blank(), axis.text = element_text(color = "black"), panel.background = element_rect(fill = "#f0f0f0"), axis.line = element_blank(), aspect.ratio = 1) +
  ggsave("Pseudo-bulk_PC1_PC2_all.pdf")

# Identify differentially expressed genes and enriched hallmark pathways in MDA cells from different tissues
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
dir.create("msigdb")

# Make table of comparisons across tissue
comparisons <- combinations(4, 2, unique(as.character(sampleType)))
colnames(comparisons) <- c("Tissue_2", "Tissue_1")
rownames(comparisons) <- 1:nrow(comparisons)

# Find DEGs for each comparison
for(i in 1:nrow(comparisons)) {
  
  # Make vectors for each tissue in the comparison
  tissue1 <- comparisons[i, 1]
  tissue2 <- comparisons[i, 2]
  print(tissue1)
  print(tissue2)
  
  # Subset expression data to those genes actively expressed in either sample type
  expression <- exprs[geneActivity[ ,tissue1] | geneActivity[ ,tissue2], c(sampleType == tissue1 | sampleType == tissue2)]
  genes <- data.frame(rowData_summed[is.element(rownames(rowData_summed), rownames(expression)), ])
  
  # Make DGEList object and normalise using TMM
  allDGEList <- DGEList(counts = expression, group = sampleType[sampleType == tissue1 | sampleType == tissue2], genes = genes)
  allDGEList <- calcNormFactors(allDGEList, method = "TMM")
  
  # Make design matrix, fit linear models, list comparisons
  type <- factor(sampleTable$sampleType, levels = c(tissue1, tissue2))
  rep <- factor(sampleTable$sampleRep, levels = c(1:4))
  design <- model.matrix(~0 + rep + type)
  rownames(design) <- sampleName[sampleType == tissue1 | sampleType == tissue2]
  allDGEList <- voom(allDGEList, design, plot = TRUE)

  # Identify significantly altered mSigDB genesets
  for(sig in names(mdbsig)) {
    camera.out <- camera(allDGEList$E, mdbsig[[sig]], design, contrast = 5, inter.gene.cor=0.01, trend.var = TRUE)
    write.csv(camera.out, paste0("msigdb/", tissue1, "_", tissue2, "_", sig, "_camera.csv"))
  }
  
  # Identify differentially expressed genes at 0 LFC threshold.
  lbFit <- lmFit(allDGEList, design)
  lbFit3 <- treat(lbFit, lfc = 0)
  DEG3 <- decideTests(lbFit3)
  print(summary(DEG3))
  results <- data.frame(topTreat(lbFit3, coef = 5, n = Inf, p.value = Inf))
  write.csv(results, paste0(tissue1, "_", tissue2, "_DEG_0LFC.csv"))
  
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
    fwrite(fgseaRes, file=paste0("msigdb/", tissue1, "_", tissue2, "_", sig, "_fgsea.csv"), sep=",", sep2=c("", " ", ""))
  } 
}
