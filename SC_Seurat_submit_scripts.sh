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

# QC and prefiltering
qsub -P OsteoporosisandTranslationalResearch -N 'Determine_thresholds_'$projectID -b y -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Determine_prefiltering_thresholds.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Prefiltering_'$projectID -b y -hold_jid 'Determine_thresholds_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Prefiltering_whole_experiment.R'

# Pseudo-bulk tissue level differences
qsub -P OsteoporosisandTranslationalResearch -N 'PseudoDGE_tissue_'$projectID -b y -hold_jid 'Prefiltering_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Pseudo-bulk_DGE_by_tissue.R'
qsub -P OsteoporosisandTranslationalResearch -N 'PseudoDGE_GO_'$projectID -b y -hold_jid 'PseudoDGE_tissue_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Pseudo-bulk_GO_by_tissue.R'
qsub -P OsteoporosisandTranslationalResearch -N 'PseudoDGE_GO_venn_'$projectID -b y -hold_jid 'PseudoDGE_tissue_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Pseudo-bulk_GO_by_tissue_venn_partitioning.R'
qsub -P OsteoporosisandTranslationalResearch -N 'PseudoDGE_hallmark_heatmaps_'$projectID -b y -hold_jid 'SeuratCluster_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Pseudo-bulk_by_tissue_hallmark_heatmaps.R'

# Whole experiment clustering and marker identification
qsub -P OsteoporosisandTranslationalResearch -N 'SeuratCluster_'$projectID -b y -hold_jid 'Prefiltering_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Whole_experiment_seurat_merge_and_clustering.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Markers_SeuratCluster_'$projectID -b y -hold_jid 'SeuratCluster_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Marker_enrichment_analysis_Seurat.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Markers_TissuevsTissue_'$projectID -b y -hold_jid 'Markers_SeuratCluster_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/TissuevsTissue_marker_enrichment_analysis_Seurat.R'
qsub -P OsteoporosisandTranslationalResearch -N 'Markers_TissuevsAll_'$projectID -b y -hold_jid 'Markers_TissuevsTissue_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/TissuevsAll_marker_enrichment_analysis_Seurat.R'

# Whole experiment QC visualisations
# qsub -P OsteoporosisandTranslationalResearch -N 'SeuratClusternoStress_'$projectID -b y -hold_jid 'Markers_TissuevsTissue_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Whole_experiment_seurat_merge_and_clustering_without_regress_stress.R'
# qsub -P OsteoporosisandTranslationalResearch -N 'SeuratClusternoCC_'$projectID -b y -hold_jid 'Markers_TissuevsTissue_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Whole_experiment_seurat_merge_and_clustering_without_regress_CC.R'

# Within tissue clustering and marker identification
# qsub -P OsteoporosisandTranslationalResearch -N 'Within_Tissue_SeuratCluster_'$projectID -b y -hold_jid 'Prefiltering_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Within_tissue_seurat_merge_and_clustering.R'
# qsub -P OsteoporosisandTranslationalResearch -N 'Markers_Within_Tissue_SeuratCluster_'$projectID -b y -hold_jid 'Within_Tissue_SeuratCluster_'$projectID -wd $logDir -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH $scriptsPath$projectID'/'$projectID'_scripts/Within_tissue_marker_enrichment_analysis_Seurat.R'

