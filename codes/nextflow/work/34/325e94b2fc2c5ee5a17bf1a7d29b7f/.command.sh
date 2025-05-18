#!/bin/bash -euo pipefail
deduplicate_bismark \
     \
    -p \
    --bam H021-E7U1N6_plasma-03-01_1_val_1_bismark_bt2_pe.bam

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_DEDUPLICATE":
    bismark: $(echo $(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*$//')
END_VERSIONS
