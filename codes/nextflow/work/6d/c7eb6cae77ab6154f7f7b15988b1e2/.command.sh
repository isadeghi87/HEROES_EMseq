#!/bin/bash -euo pipefail
bismark_methylation_extractor \
    OE0290-PED_2LB-197_plasma-02-01_1_val_1_bismark_bt2_pe.deduplicated.bam \
    --bedGraph \
    --counts \
    --gzip \
    --report \
    -p \
    --comprehensive --merge_non_CpG      --no_overlap --ignore_r2 2 --ignore_3prime_r2 2 --multicore 5 --buffer_size 98G

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_METHYLATIONEXTRACTOR":
    bismark: $(echo $(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*$//')
END_VERSIONS
