#!/bin/bash
#PBS -N SNAKEMASTER
#PBS -o logs/workflows/workflow_output-${PBS_JOBID}.log
#PBS -l nodes=1:ppn=25
#PBS -l walltime=500:00:00
#PBS -l mem=150g

WORKDIR=/omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/results/biscuit_pipeline/
cd $WORKDIR

mkdir -p logs/workflows
CONFIG_FILE=/omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/codes/biscuit_pipeline/config.yaml

# Add snakemake to PATH here
# if [[ $(which snakemake 2>/dev/null) ]]; then
#     snakemake_module="bbc2/snakemake/snakemake-7.25.0"

#     module load $snakemake_module
# fi

# Save DAG job file with a time stamp
TIME=$(date "+%Y-%m-%d_%H.%M.%S")

## load modules 
module load snakemake/5.26.1-foss-2020a-Python-3.8.2
module load  Miniconda3/4.5.11
#snakemake --configfile ${CONFIG_FILE} --dry-run         > logs/workflows/workflow_${TIME}.txt
#snakemake --configfile ${CONFIG_FILE} --dag | dot -Tpng > logs/workflows/workflow_${TIME}.png

# Run snakemake with PBS
snakemake --use-conda --cores 16 --configfile $CONFIG_FILE --rerun-incomplete

# snakemake \
#     --printshellcmds \
#     --latency-wait 20 \
#     --keep-going \
#     --use-conda \
#     --jobs 20 \
#     --configfile ${CONFIG_FILE} \
#     --cluster "mkdir -p logs/{rule}; qsub \
#         -l nodes=1:ppn={threads} \
#         -l mem={resources.mem_gb}GB \
#         -l walltime={resources.time} \
#         --rerun-incomplete \
#         -o logs/{rule}/{rule}-%j.log \
#         -j oe"
