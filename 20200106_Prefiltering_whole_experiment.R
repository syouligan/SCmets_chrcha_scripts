# --------------------------------------------------------------------------
#! Collate all data into SingleCellExperiment object, and filter cells and genes based on available data. Mouse cells and genes are removed.
# --------------------------------------------------------------------------

# Set working directory, load data and libraries
# --------------------------------------------------------------------------

#location <- "/Users/mac/cloudstor/" # if local
location <- "/share/ScratchGeneral/scoyou/" # if wolfpack
setwd(paste0(location, "sarah_projects/SCmets_chrcha/project_results/prefiltered/"))

library('DropletUtils')
library('dplyr')
library('tidyr')
library('scater')


# Load all data into single object
# --------------------------------------------------------------------------

# Load all samples
data_directory <- paste0(location, "sarah_projects/SCmets_chrcha/raw_data/data")
raw_experiment <- read10xCounts(paste0(list.files(data_directory, full.names = TRUE), '/outs/filtered_feature_bc_matrix'), col.names = TRUE)

# Create simple informative sample ids
old_ids <- c(paste0(list.files(data_directory, full.names = TRUE), '/outs/filtered_feature_bc_matrix'))
new_ids <- c("Liver_NA_1", "LN_NA_1", "Lung_NA_1", "Liver_NA_2", "Lung_NA_2", "Lymph_NA_2", "Liver_A_3", "Liver_B_3", "LN_A_3", "LN_B_3", "Lung_A_3", "Lung_B_3", "Liver_NA_4", "Lung_NA_4", "Lymph_NA_4", "Primary_NA_1", "Primary_NA_2", "Primary_NA_3", "Primary_NA_4")
new_ids <- data.frame(new_ids) %>%
  separate(new_ids, c("Tissue", "Met", "Replicate"), "_", remove = FALSE)
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

# Add gene labels and organism info.
ensembl_ids <- data.frame(rowData(raw_experiment)) %>%
  separate(ID, c("Organism", "Ensembl"), "_", remove = TRUE)
rowData(raw_experiment)$Organism <- as.character(ensembl_ids$Organism)
rowData(raw_experiment)$Ensembl <- as.character(ensembl_ids$Ensembl)
rowData(raw_experiment)$GeneSymbol <- as.character(gsub(".*_", "", rowData(raw_experiment)$Symbol))
rownames(raw_experiment) <- uniquifyFeatureNames(rowData(raw_experiment)$Ensembl, rowData(raw_experiment)$Symbol)
saveRDS(raw_experiment, "Raw_experiment_all_samples_labeled.rds")

# Filter cells based on mouse or human, mitochondrial counts outliers.
# --------------------------------------------------------------------------


# Filter genes based on mouse or human, counts across cells 
# --------------------------------------------------------------------------

