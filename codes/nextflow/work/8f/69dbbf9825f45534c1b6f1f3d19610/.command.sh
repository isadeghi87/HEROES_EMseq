#!/bin/bash -euo pipefail
samtools sort \
     \
    -@ 6 \
    -m 6144M \
    -o OE0290-PED_5LB-003_plasma-02-01.sorted.bam \
    -T OE0290-PED_5LB-003_plasma-02-01.sorted \
    OE0290-PED_5LB-003_plasma-02-01_1_val_1_bismark_bt2_pe.bam

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:SAMTOOLS_SORT_ALIGNED":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
