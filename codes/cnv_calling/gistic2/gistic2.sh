#!/bin/bash
module load miniconda/4.9.2
source activate /home/i439h/projects/Emseq_temp/tools/gistic2/gistic_env

output=/omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/results/cnv_calling/gistic2
segfiles=/omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/results/cnv_calling/cfdna/*segs
refgene=/omics/odcf/analysis/OE0290_projects_temp/pediatric_tumor/EMseq/tools/gistic2/support/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

## run gistic
gistic2 -b $output \
        -seg $segfiles \
        -refgene $refgene \
        -genegistic 1 \
        -rx 0 \
        -cap 1.5 \
        -conf 0.90 \
        -brlen 0.98 \
        -armpeel 1 \
        -savegene 1 \
        -gcm extreme
