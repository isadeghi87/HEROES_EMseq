#!/bin/bash

# Load necessary modules
module load anaconda3/2021.05
module load parallel/20210622

# Activate the conda environment
source activate cfdna

# Define directories
base_dir="/home/i439h/projects/Emseq_temp/results/nextflow/bismark/deduplicated/"
results="/home/i439h/projects/Emseq_temp/results/cnv_calling/cfdna"

# List all deduplicated BAM files in the base directory
bam_files=($(find "$base_dir" -name "*.deduplicated.bam"))

# Define common parameters for the Python command
GENOME="hg38"
BIN_SIZE=50000

# Iterate over each BAM file
for bam_file in "${bam_files[@]}"; do
    echo "Processing $bam_file..."
    python -m cfdna callCNVs --bam "$bam_file" --segs --genome "$GENOME" --bin_size "$BIN_SIZE"
done

# Move the output files to the results directory
mv *.pdf *.segs "$results"

# Deactivate the conda environment
conda deactivate
