#!/bin/bash

## running nexflow nf-core/methylseq pipeline for HEREOES-AYA OE0290 whole genome biosulfite sequencing data 
# load nextflow module 


module load nextflow/23.10.1

##Directories and files
## codes 
codes=/home/i439h/projects/Emseq_temp/codes/nextflow
OUTDIR=/omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/results/nextflow
INPUT=/home/i439h/projects/Emseq_temp/datasets/sample_sheets/sampleSheet_undone.csv
export NXF_SINGULARITY_CACHEDIR=/home/i439h/projects/Emseq_temp/codes/nextflow/work/singularity
export NXF_OPTS="-Dleveldb.mmap=false"
export NXF_OPTS='-Xms1g -Xmx4g'

cd $codes
nextflow run nf-core/methylseq --input $INPUT --outdir $OUTDIR --genome GRCh38 --aligner bismark -profile singularity --comprehensive --em_seq true  --max_cpus 20 --max_memory 180.GB  -resume

