#!/bin/bash
module load anaconda3/2021.05
source activate /home/i439h/projects/EMseq/tools/biscuit/conda
module load samtools/1.9

#files and dirs
sample='I070-032_plasma-02-01'
outdir=/home/i439h/projects/EMseq/results/snv_biscuit/${sample}
bam=/home/i439h/projects/EMseq/results/biscuit_alignment/${sample}/$sample.bam


mkdir -p $outdir

output=${outdir}/${sample}.vcf
reference=/home/i439h/projects/pools/AG_Thongjuea/Dataset/10x_multiomics/references/human/refdata-gex-GRCh38-2024-A/fasta/genome.fa

cd $outdir

## run biscuit to call variants
#biscuit index $reference
# biscuit pileup -@ 10 -o $output $reference $bam
 bgzip -@ NTHREADS $output
 tabix -p vcf $output.gz

# Extract DNA methylation into BED format
# Also compresses and indexes the BED
biscuit vcf2bed $output.gz > ${sample}_methylation.bed
bgzip ${sample}_methylation.bed
tabix -p bed ${sample}_methylation.bed.gz

