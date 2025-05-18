#!/bin/bash
module load anaconda3/2021.05
source activate hatchet

# Path to reference genome - make sure you have also generated the reference dictionary as /path/to/reference.dict
reference = "/home/i439h/projects/pools/AG_Thongjuea/Code/10x_multiomics/Sadeghi/Demultiplexing/references/hg38/genome.fa"
reference_version = "hg38"
normal = "/home/i439h/projects/EMseq/results/nextflow/bismark/deduplicated/OE0290-PED_0LB-060_plasma-03-02_1_val_1_bismark_bt2_pe.deduplicated.bam.sorted.bam"
# Space-delimited list of tumor BAM locations
bams = "/home/i439h/projects/EMseq/results/nextflow/bismark/deduplicated/OE0290-PED_I070-032_plasma-02-01_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam"
# Space-delimited list of tumor names
samples = "I070-032"

mincov = 2
maxcov = 300

#genotype-snps

python3 -m hatchet genotype-snps -N ${normal} -r ${reference} -j ${J} -c ${mincov} -C ${maxcov} \
                            -o ${SNP} |& tee ${BAF}bafs.log

