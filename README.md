# HEROES_EMseq
Codes for analysis of EMSeq data from liquid biopsy

1. This repository contain codes and results for EMseq data analysis
The main directory is at /omics/odcf/analysis/OE0290_projects/pediatric_tumor/whole_genome_biosulfite_sequencing/
A temporary directory was created by ODCF at /omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/ which contains most updated codes and results. 

A. codes folder contains the following:
1. nextflow: for running nf-core/methylseq pipeline
2. dmr: for differentially methylated region analysis
3. qc_summary: Summarizing QC of EMseq samples
4. sample_prep: preparation of input for running the pipeline
5. cnv_calling: Running CNV calling. In this folder cfdna is the main tool that works better for liquid biopsy (i.e. other folders or tools should be ignored)

B. data folder contains normalized loaded methylation calls for DMR analysis

C. datasets folder contains samples_sheets used as input for pipeline.

D. The results folder contains all results corresponding to each analysis:
The output from nextflow pipeline are stored in nextflow folder. The most important outputs are in bismark/deduplicated and bismark/methylation_calls. The bam files from bismark/deduplicated can be used for CNV calling. The output from bismark/methylation_calls can be used for DMR analysis.

Figures folder contains figures from dmr analysis. 

biscuit_pipeline folder contains test scripts for running biscuit tool. This should be modified  and rerun to replace nextflow pipeline.



