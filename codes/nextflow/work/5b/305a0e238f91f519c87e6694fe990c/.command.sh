#!/bin/bash -euo pipefail
bismark2summary H021-E7U1N6_plasma-03-01_1_val_1_bismark_bt2_pe.bam 2LB-269_plasma-01-01_1_val_1_bismark_bt2_pe.bam 2LB-269_plasma-02-01_1_val_1_bismark_bt2_pe.bam 5LB-003_plasma-01-01_1_val_1_bismark_bt2_pe.bam 2LB-304_plasma-01-01_1_val_1_bismark_bt2_pe.bam 2LB-257_plasma-01-01_1_val_1_bismark_bt2_pe.bam 2LB-077_plasma-01-01_1_val_1_bismark_bt2_pe.bam 2LB-197_plasma-03-01_1_val_1_bismark_bt2_pe.bam H021-PJCDQK_plasma-03-01_1_val_1_bismark_bt2_pe.bam 5LB-017_plasma-01-01_1_val_1_bismark_bt2_pe.bam

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:BISMARK:BISMARK_SUMMARY":
    bismark: $(echo $(bismark -v 2>&1) | sed 's/^.*Bismark Version: v//; s/Copyright.*$//')
END_VERSIONS
