#!/bin/bash -euo pipefail
[ ! -f  OE0290-PED_I034-042_plasma-15-01_1.fastq.gz ] && ln -s AS-1052158-LR-69347_R1.fastq.gz OE0290-PED_I034-042_plasma-15-01_1.fastq.gz
[ ! -f  OE0290-PED_I034-042_plasma-15-01_2.fastq.gz ] && ln -s AS-1052158-LR-69347_R2.fastq.gz OE0290-PED_I034-042_plasma-15-01_2.fastq.gz
trim_galore \
    --fastqc   --clip_r1 10 --clip_r2 10 --three_prime_clip_r1 10 --three_prime_clip_r2 10 \
    --cores 8 \
    --paired \
    --gzip \
    OE0290-PED_I034-042_plasma-15-01_1.fastq.gz \
    OE0290-PED_I034-042_plasma-15-01_2.fastq.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE":
    trimgalore: $(echo $(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*$//')
    cutadapt: $(cutadapt --version)
END_VERSIONS
