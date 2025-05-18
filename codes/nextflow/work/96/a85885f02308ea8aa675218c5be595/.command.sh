#!/bin/bash -euo pipefail
samtools sort \
     \
    -@ 6 \
    -m 6144M \
    -o 5LB-003_plasma-01-01.sorted.bam \
    -T 5LB-003_plasma-01-01.sorted \
    5LB-003_plasma-01-01_1_val_1_bismark_bt2_pe.bam

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:SAMTOOLS_SORT_ALIGNED":
    samtools: $(echo $(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*$//')
END_VERSIONS
