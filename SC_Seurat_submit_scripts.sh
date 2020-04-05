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

qsub -P OsteoporosisandTranslationalResearch -N 'Determine_thresholds_'$projectID -b y -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Determine_prefiltering_thresholds.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Prefiltering_'$projectID -b y -hold_jid 'Determine_thresholds_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Prefiltering_whole_experiment.R'
qsub -P OsteoporosisandTranslationalResearch -N 'PseudoDGE_tissue_'$projectID -b y -hold_jid 'Prefiltering_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 5 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Pseudo-bulk_DGE_by_tissue.R'
qsub -P OsteoporosisandTranslationalResearch -N 'SeuratCluster_'$projectID -b y -hold_jid 'Prefiltering_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Whole_experiment_seurat_merge_and_clustering.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Within_Tissue_SeuratCluster_'$projectID -b y -hold_jid 'Prefiltering_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Within_tissue_seurat_merge_and_clustering.R'

