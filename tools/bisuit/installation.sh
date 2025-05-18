#!/bin/bash
module load miniconda/4.9.2

# Replace /path/to/your/directory with your desired directory path
conda create --prefix /home/i439h/projects/Emseq_temp/tools/bisuit python=3.8 -y

# Activate the environment using the full path
source activate /home/i439h/projects/Emseq_temp/tools/bisuit

# Install BISCUIT
conda install -c bioconda biscuit -y
# Install Snakemake
conda install -c bioconda -c conda-forge snakemake -y
