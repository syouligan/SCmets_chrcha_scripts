#!/bin/bash

# Load R 3.6.0
module load briglo/R/3.6.0

qsub -P OsteoporosisandTranslationalResearch -N "Prefiltering" -b y -wd /share/ScratchGeneral/scoyou/sarah_projects/SCmets_chrcha/logs -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH /home/scoyou/sarah_projects/SCmets_chrcha/SCmets_chrcha_scripts/20200106_Prefiltering_whole_experiment.R
qsub -P OsteoporosisandTranslationalResearch -N "NormyCluster" -b y -hold_jid Prefiltering -wd /share/ScratchGeneral/scoyou/sarah_projects/SCmets_chrcha/logs -j y -R y -l mem_requested=8G -pe smp 64 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH /home/scoyou/sarah_projects/SCmets_chrcha/SCmets_chrcha_scripts/20200113_Whole_experiment_batch_correction_and_clustering.R
qsub -P OsteoporosisandTranslationalResearch -N "Cell_Cycle" -b y -hold_jid NormyCluster -wd /share/ScratchGeneral/scoyou/sarah_projects/SCmets_chrcha/logs -j y -R y -l mem_requested=8G -pe smp 64 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH /home/scoyou/sarah_projects/SCmets_chrcha/SCmets_chrcha_scripts/20200125_Whole_experiment_cell_cycle_annotation.R
qsub -P OsteoporosisandTranslationalResearch -N "WithinTissue" -b y -hold_jid Cell_Cycle -wd /share/ScratchGeneral/scoyou/sarah_projects/SCmets_chrcha/logs -j y -R y -l mem_requested=8G -pe smp 32 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH /home/scoyou/sarah_projects/SCmets_chrcha/SCmets_chrcha_scripts/20200125_Within_tissue_batch_correction_and_clustering.R
qsub -P OsteoporosisandTranslationalResearch -N "MAGIC" -b y -hold_jid Cell_Cycle -wd /share/ScratchGeneral/scoyou/sarah_projects/SCmets_chrcha/logs -j y -R y -l mem_requested=8G -pe smp 64 -V -m bea -M s.youlten@garvan.org.au R CMD BATCH /home/scoyou/sarah_projects/SCmets_chrcha/SCmets_chrcha_scripts/20200125_MAGIC_imputation.R
