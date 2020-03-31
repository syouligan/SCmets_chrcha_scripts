#!/bin/bash

# Load R 3.6.0
module load briglo/R/3.6.0

# Set paths
projectID='SCmets_chrcha'
resultsPath='/share/ScratchGeneral/scoyou/sarah_projects/'
scriptsPath='/home/scoyou/sarah_projects/'

# Make log path
logDir=$scriptsPath$projectID'/'$projectID'_scripts/logs'
mkdir -p $logDir

#qsub -P OsteoporosisandTranslationalResearch -N 'Determine_thresholds_'$projectID -b y -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Determine_prefiltering_thresholds.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Prefiltering_'$projectID -b y -hold_jid 'Determine_thresholds_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Prefiltering_whole_experiment.R'
qsub -P OsteoporosisandTranslationalResearch -N 'PseudoDGE_tissue_'$projectID -b y -hold_jid 'Prefiltering_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Pseudo-bulk_DGE_by_tissue.R'
qsub -P OsteoporosisandTranslationalResearch -N 'NormyCluster_'$projectID -b y -hold_jid 'Prefiltering_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Whole_experiment_batch_correction_and_clustering.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Marker_genes_'$projectID -b y -hold_jid 'NormyCluster_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Marker_gene_detection.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Differential_abundance_'$projectID -b y -hold_jid 'NormyCluster_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Differential_abundance.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Cell_Cycle_'$projectID -b y -hold_jid 'NormyCluster_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Whole_experiment_cell_cycle_annotation.R'
qsub -P OsteoporosisandTranslationalResearch -N 'SeuratCluster_'$projectID -b y -hold_jid 'Cell_Cycle_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Whole_experiment_seurat_merge_and_clustering.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Marker_genes_seurat_'$projectID -b y -hold_jid 'SeuratCluster_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Marker_gene_detection_seurat.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Differential_abundance_seurat_'$projectID -b y -hold_jid 'SeuratCluster_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Differential_abundance_seurat.R'

