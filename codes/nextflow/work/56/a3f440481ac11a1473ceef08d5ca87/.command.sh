#!/bin/bash -euo pipefail
samtools sort \
     \
    -@ 6 \
    -m 6144M \
    -o OE0290-PED_2LB-323_plasma-01-01.sorted.bam \
    -T OE0290-PED_2LB-323_plasma-01-01.sorted \
    OE0290-PED_2LB-323_plasma-01-01_1_val_1_bismark_bt2_pe.deduplicated.bam

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:SAMTOOLS_SORT_DEDUPLICATED":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
