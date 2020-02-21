#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Determine cell filtering thresholds using SAVER standard error, library size, mito content, etc
# --------------------------------------------------------------------------

# Set working directory, load data and libraries
# --------------------------------------------------------------------------

# Working directory
# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  location <- "/Users/mac/cloudstor/"
  place <- "local"
} else {
  location <- "/share/ScratchGeneral/scoyou/"
  place <- "wolfpack"
}

setwd(paste0(location, "sarah_projects/SCMDA231mets_chrcha/project_results/prefiltering/"))

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('EnsDb.Hsapiens.v86')
library('grid')
library('ggplot2')
library('ggridges')
library('viridis')
library('Matrix')
library('SAVER')
library('stringr')

# Load all data into single object
# --------------------------------------------------------------------------

# Load sample metadata
sample_IDs <- read.csv(paste0(location, "sarah_projects/SCMDA231mets_chrcha/sample_metadata/Sample_ids.csv"))

# Load all samples and save before editing
data_directory <- paste0(location, "sarah_projects/SCMDA231mets_chrcha/raw_data/data")

if(file.exists("all_data/Raw_experiment_all_samples_old_labels.rds")) {
  raw_experiment <- readRDS("all_data/Raw_experiment_all_samples_old_labels.rds")
  } else {
  raw_experiment <- read10xCounts(paste0(list.files(data_directory, full.names = TRUE), '/outs/filtered_feature_bc_matrix'), col.names = TRUE, sample.names = list.files(data_directory, full.names = FALSE))
  if (place == "wolfpack") {
    saveRDS(raw_experiment, "all_data/Raw_experiment_all_samples_old_labels.rds")
  } else {
    print("Working local")
  }
}

# Rename colData Sample slot. Add Tissue, Met, Replicate information.
raw_experiment$old_Sample <- raw_experiment$Sample

idx <- match(raw_experiment$Sample, sample_IDs$Old_sample_ID)
raw_experiment$Sample <- as.character(sample_IDs$New_sample_ID) [idx]
raw_experiment$Tissue <- as.character(sample_IDs$Tissue) [idx]
raw_experiment$Replicate <- as.character(sample_IDs$Replicate) [idx]
raw_experiment$cellIDs <- rownames((colData(raw_experiment)))

# Add gene labels and organism info.
ensembl_ids <- data.frame(rowData(raw_experiment))
ensembl_ids$ID <- str_replace_all(ensembl_ids$ID, "___", "_")
rowData(raw_experiment)$GeneSymbol <- str_replace_all(ensembl_ids$ID, "___", "_")
ensembl_ids <- ensembl_ids %>%
  separate(ID, c("Organism", "Ensembl"), sep = "([\\_+])", remove = TRUE)
rowData(raw_experiment)$Organism <- as.character(ensembl_ids$Organism)
rowData(raw_experiment)$Ensembl <- as.character(ensembl_ids$Ensembl)
rowData(raw_experiment)$GeneSymbol <- as.character(gsub(".*_", "", rowData(raw_experiment)$Symbol))
rownames(raw_experiment) <- uniquifyFeatureNames(rowData(raw_experiment)$Ensembl, rowData(raw_experiment)$GeneSymbol)

# Remove mouse cells and mouse genes
# --------------------------------------------------------------------------

# Remove mouse cells
stats <- perCellQCMetrics(raw_experiment, subsets=list(Human=which(rowData(raw_experiment)$Organism == "GRCh38")))
raw_experiment$Human_percent <- stats$subsets_Human_percent
raw_experiment$Human_cells <- raw_experiment$Human_percent >= 99
ggplot(data.frame(colData(raw_experiment)), aes(x = Human_percent, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Human_percent_ridge_raw.pdf", useDingbats = FALSE)
raw_experiment <- raw_experiment[,which(raw_experiment$Human_cells)]

# Remove mouse genes
raw_experiment <- raw_experiment[which(rowData(raw_experiment)$Organism == "GRCh38"),]

# Calculate SAVER stats of experiment. Remove genes with zero counts to speed up process.
# --------------------------------------------------------------------------

# Remove genes that have zero counts across all cells
undetectedGenes <- rowSums(counts(raw_experiment)) == 0
raw_experiment <- raw_experiment[!undetectedGenes, ] # Remove undetected genes

# Calculate SAVER standard error estimates and plot relative to library size. Likely observed two distinct populations, one with low library size and high standard error, the other with higher library size and lower standard error (set filter parameters to capture this population)
row_names <- split(rownames(raw_experiment), ceiling(seq_along(rownames(raw_experiment))/2500))
names(row_names)

for(i in names(row_names)){
  saver.new <- saver(counts(raw_experiment), pred.genes = which(rownames(raw_experiment) %in% row_names[[i]]), pred.genes.only = TRUE, do.fast = FALSE, ncores = 32)
  if(exists(saver.out)){
    saver.out <- combine.saver(list(saver.out, saver.new))
  } else {
    saver.out <- saver.new
    }
  } # Split genes into multiple groups of 2500 to avaoid memory issues and then recombine into single object
saveRDS(saver.out, "SAVER_run_raw_experiment.rds")

pdf("Unfiltered_SAVER_SEvsLibSize.pdf")
plot(colSums(counts(raw_experiment)), colMeans(saver.out$se), log="x")
dev.off()

raw_experiment$meanSAVERse <- colMeans(saver.out$se) # Use to filter later after viewing plots

# Plot cells based on distribution of Library size, gene number and mito content
# --------------------------------------------------------------------------

# Idenitify cells to discard based on 3MAD outlier in either number of detected genes, library size or >20% mitochondrial content
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(raw_experiment)$Ensembl, column="SEQNAME", keytype="GENEID")
stats <- perCellQCMetrics(raw_experiment, subsets=list(Mito=which(location=="MT")))
raw_experiment$Lib_size <- stats$sum
raw_experiment$Genes_detected <- stats$detected
raw_experiment$Mito_percent <- stats$subsets_Mito_percent

# Plot QC stats
ggplot(data.frame(colData(raw_experiment)), aes(x = Lib_size, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Library_size_ridge_raw.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(raw_experiment)), aes(x = Genes_detected, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Number_of_genes_ridge_raw.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(raw_experiment)), aes(x = Mito_percent, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Mito_percent_ridge_raw.pdf", useDingbats = FALSE)

# Save practice dataset (5% of cells from each sample)
raw_experiment$cellIDs <- rownames((colData(raw_experiment)))
fe_subset <- data.frame(data.frame(colData(raw_experiment)) %>%
                          group_by(Sample) %>%
                          sample_frac(0.05))
raw_experiment$Practice_subset <- is.element(rownames(colData(raw_experiment)), fe_subset$cellIDs)
practice_exp <- raw_experiment[,which(raw_experiment$Practice_subset)]

if (place == "wolfpack") {
  saveRDS(practice_exp, "practice_all_data/Raw_experiment_all_samples_labeled_practice.rds")
} else {
  print("Working local")
}

# Save total raw dataset
if (place == "wolfpack") {
  saveRDS(raw_experiment, "all_data/Raw_experiment_all_samples_labeled.rds")
} else {
  print("Working local")
}