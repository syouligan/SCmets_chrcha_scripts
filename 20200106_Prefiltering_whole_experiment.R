#!/usr/bin/Rscript

# --------------------------------------------------------------------------
#! Collate all data into SingleCellExperiment object, and filter cells and genes based on available data. Mouse cells and genes are removed.
# --------------------------------------------------------------------------

# Set working directory, load data and libraries
# --------------------------------------------------------------------------

# Working directory
if(dir.exists("/Users/mac/cloudstor/")) {
  location <- "/Users/mac/cloudstor/"
} else {
  location <- "/share/ScratchGeneral/scoyou/"
}

setwd(paste0(location, "sarah_projects/SCmets_chrcha/project_results/prefiltered/"))

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

# Load all samples
data_directory <- paste0(location, "sarah_projects/SCmets_chrcha/raw_data/data")
raw_experiment <- read10xCounts(paste0(list.files(data_directory, full.names = TRUE), '/outs/filtered_feature_bc_matrix'), col.names = TRUE)
saveRDS(raw_experiment, "Raw_experiment_all_samples_old_labels.rds")

# Create simple informative sample ids
old_ids <- c(paste0(list.files(data_directory, full.names = TRUE), '/outs/filtered_feature_bc_matrix'))
new_ids <- c("Liver_1", "LN_1", "Lung_1", "Liver_2", "Lung_2", "LN_2", "Liver_3", "Liver_3", "LN_3", "LN_3", "Lung_3", "Lung_3", "Liver_4", "Lung_4", "LN_4", "Primary_1", "Primary_2", "Primary_3", "Primary_4")
new_ids <- data.frame(new_ids) %>%
  separate(new_ids, c("Tissue", "Replicate"), "_", remove = FALSE)
ids <- data.frame(cbind(old_ids, new_ids))
print(ids) # check id mapping
write.csv(ids, paste0(location, "sarah_projects/SCmets_chrcha/sample_metadata/20200106_Sample_ids.csv"), row.names = FALSE)

# Rename colData Sample slot. Add Tissue, Met, Replicate information.
raw_experiment$old_Sample <- raw_experiment$Sample
idx <- match(raw_experiment$Sample, ids$old_ids)
raw_experiment$Sample <- as.character(ids$new_ids) [idx]
raw_experiment$Tissue <- as.character(ids$Tissue) [idx]
raw_experiment$Met <- as.character(ids$Met) [idx]
raw_experiment$Replicate <- as.character(ids$Replicate) [idx]
raw_experiment$cellIDs <- rownames((colData(raw_experiment)))

# Add gene labels and organism info.
ensembl_ids <- data.frame(rowData(raw_experiment)) %>%
  separate(ID, c("Organism", "Ensembl"), "_", remove = TRUE)
rowData(raw_experiment)$Organism <- as.character(ensembl_ids$Organism)
rowData(raw_experiment)$Ensembl <- as.character(ensembl_ids$Ensembl)
rowData(raw_experiment)$GeneSymbol <- as.character(gsub(".*_", "", rowData(raw_experiment)$Symbol))
rownames(raw_experiment) <- uniquifyFeatureNames(rowData(raw_experiment)$Ensembl, rowData(raw_experiment)$GeneSymbol)
saveRDS(raw_experiment, "Raw_experiment_all_samples_labeled.rds")

# Filter genes and cells based on mouse or human, counts outliers, mitochondrial outliers.
# --------------------------------------------------------------------------

# Remove mouse cells
stats <- perCellQCMetrics(raw_experiment, subsets=list(Human=which(rowData(raw_experiment)$Organism == "hg19")))
raw_experiment$Human_percent <- stats$subsets_Human_percent
raw_experiment$Human_cells <- raw_experiment$Human_percent >= 99
plotColData(raw_experiment, x="Sample", y="Human_percent", colour_by="Human_cells", other_fields="Tissue") +
  facet_wrap(~Tissue) +
  ggtitle("Total count") +
  ggsave("Human_percent_distribution.pdf", useDingbats = FALSE)

filtered_exp <- raw_experiment[,which(raw_experiment$Human_cells)]

# Remove mouse genes
filtered_exp <- filtered_exp[which(rowData(filtered_exp)$Organism == "hg19"),]

# Idenitify cells to discard based on 3MAD outlier in either mito-content, number of detected genes, library size
location <- mapIds(EnsDb.Hsapiens.v75, keys=rowData(filtered_exp)$Ensembl, column="SEQNAME", keytype="GENEID")
stats <- perCellQCMetrics(filtered_exp, subsets=list(Mito=which(location=="MT")))
filtered_exp$Lib_size <- stats$sum
filtered_exp$Genes_detected <- stats$detected
filtered_exp$Mito_percent <- stats$subsets_Mito_percent

discard <- quickPerCellQC(stats, percent_subsets=c("subsets_Mito_percent"), batch=filtered_exp$Sample)
filtered_exp$discard_LibGenes <- discard$low_lib_size | discard$low_n_features
filtered_exp$discard_Mito <- discard$high_subsets_Mito_percent
filtered_exp$discard <- discard$discard
discard_stats <- DataFrame(colSums(as.matrix(discard)))
colnames(discard_stats) <- "Cell#"
write.csv(discard_stats, "Discard_stats.csv", row.names = FALSE)

# Plot QC stats
ggplot(data.frame(colData(filtered_exp)), aes(x = Lib_size, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Library_size_ridge.pdf", useDingbats = FALSE)

plotColData(filtered_exp, x="Sample", y="Lib_size", colour_by="discard", other_fields="Tissue") +
  facet_wrap(~Tissue) +
  scale_y_log10() +
  ggtitle("Total count") +
  ggsave("Library_size_violin.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(filtered_exp)), aes(x = Genes_detected, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Number_of_genes_ridge.pdf", useDingbats = FALSE)
  
plotColData(filtered_exp, x="Sample", y="Genes_detected", colour_by="discard", other_fields="Tissue") +
  facet_wrap(~Tissue) +
  scale_y_log10() +
  ggtitle("Detected features") +
  ggsave("Number_of_genes_violin.pdf", useDingbats = FALSE)

ggplot(data.frame(colData(filtered_exp)), aes(x = Mito_percent, y = Sample, fill = Tissue)) +
  geom_density_ridges() +
  theme_minimal() +
  ggsave("Mito_percent_ridge.pdf", useDingbats = FALSE)
  
plotColData(filtered_exp, x="Sample", y="Mito_percent", colour_by="discard", other_fields="Tissue") + 
  facet_wrap(~Tissue) +
  ggtitle("Mito percent") +
  ggsave("Mito_percent_violin.pdf", useDingbats = FALSE)

# Remove "discard" cells
filtered_exp <- filtered_exp[ ,which(!filtered_exp$discard)]

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
GOI$Any_Active <- rowSums(GOI) > 0 # Above three in any replicate of any tissue type
rowData(filtered_exp) <- cbind(rowData(filtered_exp), GOI) # Add to rowData
filtered_exp <- filtered_exp[which(GOI$Any_Active), ] # Subset to active genes

colSums(GOI) # Number of genes "active" in each tumour location

# Number of cells remaining per sample
cells_remaining <- data.frame(table(filtered_exp$Sample))
write.csv(cells_remaining, "Number_cells_remaining.csv", row.names = FALSE)

ggplot(data=cells_remaining, aes(x=Var1, y=Freq)) +
  geom_bar(stat="identity", color="black") +
  ggtitle("Cells remaining") +
  scale_color_viridis(discrete=TRUE) +
  coord_flip() +
  theme_classic() +
  ggsave("Cells_remaining.pdf", useDingbats = FALSE)

# Save datasets.
# --------------------------------------------------------------------------

# Save practice dataset (5% of cells from each sample)
filtered_exp$cellIDs <- rownames((colData(filtered_exp)))
fe_subset <- data.frame(data.frame(colData(filtered_exp)) %>%
                        group_by(Sample) %>%
                        sample_frac(0.05))
filtered_exp$Practice_subset <- is.element(rownames(colData(filtered_exp)), fe_subset$cellIDs)
practice_exp <- filtered_exp[,which(filtered_exp$Practice_subset)]
saveRDS(practice_exp, "practice_all_data/Prefiltered_experiment_Practice.rds")

raw_experiment$cellIDs <- rownames((colData(raw_experiment)))
fe_subset <- data.frame(data.frame(colData(raw_experiment)) %>%
                          group_by(Sample) %>%
                          sample_frac(0.05))
raw_experiment$Practice_subset <- is.element(rownames(colData(raw_experiment)), fe_subset$cellIDs)
practice_exp <- raw_experiment[,which(raw_experiment$Practice_subset)]
saveRDS(practice_exp, "practice_all_data/Raw_experiment_Practice.rds")


# Save total filtered dataset
saveRDS(filtered_exp, "all_data/Prefiltered_experiment_All.rds")

# Save individual samples
for(i in unique(colData(filtered_exp)$Sample)) {
  saveRDS(filtered_exp[,filtered_exp$Sample == i], paste0("individual/", i,"/Prefiltered_experiment_", i, ".rds"))
}

