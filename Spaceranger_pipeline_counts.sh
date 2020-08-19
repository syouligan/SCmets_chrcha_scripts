#!/bin/bash

## Run spaceranger pipeline to get counts
#######

# Number of cores
ncores=20

# Experimental info (update for single or paired reads, and readlength). Set paired to "yes" if appropriate. NOTE: reads are assumed to be stranded.
project="SCMDA231spatial_chrcha"

homedir="/share/ScratchGeneral/scoyou/sarah_projects"
data="data/fastqs"
images="data/slides"
tool="spaceranger"
results="project_results"
QCDir="logs"
reference="/directflow/GWCCGPipeline/projects/reference/refdata-cellranger-GRCh38-and-mm10-3.1.0"

# Path for log files
logDir=$homedir/$project/$QCDir
mkdir -p $logDir
echo "logDir $logDir"

# Path to folder containing sample folders
sample_Path=$homedir/$project/$data
echo "sample_Path "$sample_Path
image_Path=$homedir/$project/$images
echo "image_Path "$image_Path


# Make an array containing names of directories containing samples
sample_arr=( $(ls $sample_Path) )
echo "# in samples array ${#sample_arr[@]}"
echo "names in samples array ${sample_arr[@]}"

# Submit command for each sample in array
for sample in ${sample_arr[@]}; do

# Runs loop for only the first sample in array (used for development)
# for sample in ${sample_arr[0]}; do

# Define input directory, define and make output and log directories
inPath=$sample_Path/$sample
echo "inPath $inPath"

outDir=$homedir/$project/$results/$tool/$sample
mkdir -p $outDir
echo "outDir $outDir"

# Command to be executed.
Command="spaceranger count \
  --id="$sample" \
  --fastqs="$sample_Path" \
  --transcriptome="$reference" \
  --sample="$sample" \
  --image="$image_Path/$sample.jpg" \
  --unknown-slide \
  --jobmode=local \
  --localcores=16 \
  --localmem=256"

# Submit to queue
echo "Command "$Command
#$Command
qsub -P OsteoporosisandTranslationalResearch -N $tool$sample -b y -wd $logDir -j y -R y -l mem_requested=8G -pe smp $ncores -V -m bea -M s.youlten@garvan.org.au $Command

done
