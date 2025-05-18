## CNV calling for EMseq using cn.mops

setwd("~/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/home/projects/EMseq/")

## libraries
library(cn.mops)
bamfiles= c( "results/nextflow/bismark/deduplicated/OE0290-PED_0LB-064_plasma-01-01_1_val_1_bismark_bt2_pe.deduplicated.bam.sorted.bam",
            "results/nextflow/bismark/deduplicated/OE0290-PED_I070-032_plasma-02-01_1_val_1_bismark_bt2_pe.deduplicated.sorted.bam"  )

## get read counts
chrs = paste0('chr',1:22)
bamDataRanges <- getReadCountsFromBAM(bamfiles, sampleNames=c("normal",'tumor'),refSeqNames = chrs,parallel = 4)
