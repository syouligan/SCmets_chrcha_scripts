#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Collate all data into SingleCellExperiment object, and filter cells and genes based on available data. Mouse cells and genes are removed.
# --------------------------------------------------------------------------

# Set working directory, load data and libraries
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  location <- "/Users/mac/cloudstor/"
  place <- "local"
} else {
  location <- "/share/ScratchGeneral/scoyou/"
  place <- "wolfpack"
}

setwd(paste0(location, "sarah_projects/SCMDA231mets_chrcha/project_results/prefiltered/"))

# Libraries
library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')
library('EnsDb.Hsapiens.v75')
library('grid')
library('ggplot2')
library('ggridges')
library('viridis')
library('Matrix')

# Load all data into single object
# --------------------------------------------------------------------------

if (place == "local") {
  
  raw_experiment <- readRDS("practice_all_data/Raw_experiment_all_samples_labeled_practice.rds") # Load practice dataset generated after the following steps
  
} else {

  raw_experiment <- readRDS("all_data/Raw_experiment_all_samples_labeled.rds") # Load practice dataset generated after the following steps

}

# Filter genes and cells based on SE, library size thresholds and mitochondrial content (>20%).
# --------------------------------------------------------------------------

# Idenitify cells to discard based on 3MAD outlier in either number of detected genes, library size or >20% mitochondrial content
discard <- quickPerCellQC(stats, batch=raw_experiment$Sample)
raw_experiment$discard_NumGenes <- discard$low_n_features
raw_experiment$discard_SE <- raw_experiment$meanSAVERse < XXX
raw_experiment$discard_LibSize <- raw_experiment$Lib_size < XXX
raw_experiment$discard_Mito <- raw_experiment$Mito_percent > 20
raw_experiment$discard <- raw_experiment$discard_NumGenes | raw_experiment$discard_SE | raw_experiment$discard_LibSize | raw_experiment$discard_Mito

sum(raw_experiment$discard_NumGenes)
sum(raw_experiment$discard_SE)
sum(raw_experiment$discard_LibSize)
sum(raw_experiment$discard_Mito)
sum(raw_experiment$discard)

filtered_exp <- raw_experiment[ ,which(!raw_experiment$discard)]

# Plot QC stats
ggplot(data.frame(colData(filtered_exp)), aes(x = Human_percent, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Human_percent_ridge_filtered.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(filtered_exp)), aes(x = Lib_size, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Library_size_ridge_filtered.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(filtered_exp)), aes(x = Genes_detected, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Number_of_genes_ridge_filtered.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(filtered_exp)), aes(x = Mito_percent, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Mito_percent_ridg_filterede.pdf", useDingbats = FALSE)

# Remove genes without counts in at least 3 cells in each samples for any tissue
tmp_structure <- data.frame(unique(colData(filtered_exp)[ ,c("Tissue", "Sample")]))
tmp_table <- data.frame(table(tmp_structure$Tissue))
colnames(tmp_table) <- c("Tissue", "Reps")

GOI <- c()
for(t in tmp_table$Tissue) {
  tissue_exp <- filtered_exp[,filtered_exp$Tissue == t]
  sampleGT3 <- c()
  for(i in unique(tissue_exp$Sample)) {
    tmp_sample <- tissue_exp[,tissue_exp$Sample == i]
    sampleGT3 <- cbind(sampleGT3, rowSums(counts(tmp_sample) > 0) > 3) # Greater than 0 counts in at least 3 cells in sample
  }
  GOI <- cbind(GOI, rowSums(sampleGT3) == tmp_table[tmp_table$Tissue == t, "Reps", drop = TRUE]) # Above 3 in each replicate of each tissue type
}
colnames(GOI) <- paste0(tmp_table$Tissue, "_active")
GOI <- data.frame(GOI)
GOI$Any_Active <- rowSums(GOI) > 0 # Above three in every replicate of any tissue type
rowData(filtered_exp) <- cbind(rowData(filtered_exp), GOI) # Add to rowData
filtered_exp <- filtered_exp[which(GOI$Any_Active), ] # Subset to active genes

colSums(GOI) # Number of genes "active" in each tumour location

# Number of cells remaining per sample
cells_remaining <- data.frame(table(filtered_exp$Sample))
write.csv(cells_remaining, "Cells_remaining_number.csv", row.names = FALSE)

ggplot(data=cells_remaining, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color="black") +
  ggtitle("Cells remaining") +
  scale_color_viridis(discrete=TRUE) +
  coord_flip() +
  theme_classic() +
  ggsave("Cells_remaining_barplot.pdf", useDingbats = FALSE)

# Save datasets.
# --------------------------------------------------------------------------

# Save practice dataset (5% of cells from each sample)
filtered_exp$cellIDs <- rownames((colData(filtered_exp)))
fe_subset <- data.frame(data.frame(colData(filtered_exp)) %>%
                          group_by(Sample) %>%
                          sample_frac(0.05))
filtered_exp$Practice_subset <- is.element(rownames(colData(filtered_exp)), fe_subset$cellIDs)
practice_exp <- filtered_exp[,which(filtered_exp$Practice_subset)]

if (place == "wolfpack") {
  saveRDS(practice_exp, "practice_all_data/Prefiltered_experiment_practice.rds")
} else {
  print("Working local")
}

# Save total filtered dataset
if (place == "wolfpack") {
  saveRDS(filtered_exp, "all_data/Prefiltered_experiment_all.rds")
} else {
  print("Working local")
}
