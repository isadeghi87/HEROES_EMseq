#!/bin/bash -euo pipefail
deduplicate_bismark \
     \
    -p \
    --bam OE0290-PED_5LB-049_plasma-01-01_1_val_1_bismark_bt2_pe.bam

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_DEDUPLICATE":
    bismark: $(echo $(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*$//')
END_VERSIONS
