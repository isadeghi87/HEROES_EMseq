
#!/bin/bash
cd /home/i439h/projects/EMseq/tools/hatchet

module load anaconda3/2021.05
# Create a new conda environment named 'hatchet_env' in the specified directory
conda create -n hatchet hatchet

# Activate the new environment
source activate hatchet
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# Install HATCHet and its dependencies using pip or conda
conda install hatchet
conda install -c dranew shapeit