#!/bin/bash -euo pipefail
bismark \
    -1 OE0290-PED_2LB-323_plasma-01-01_1_val_1.fq.gz -2 OE0290-PED_2LB-323_plasma-01-01_2_val_2.fq.gz \
    --genome BismarkIndex \
    --bam \
    --bowtie2         --maxins 1000 --multicore 5

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_ALIGN":
    bismark: $(echo $(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*$//')
END_VERSIONS
