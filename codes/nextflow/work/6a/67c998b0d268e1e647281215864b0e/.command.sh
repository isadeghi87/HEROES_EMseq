#!/bin/bash -euo pipefail
[ ! -f  H021-E7U1N6_plasma-03-01_1.fastq.gz ] && ln -s AS-1383049-LR-78290_R1.fastq.gz H021-E7U1N6_plasma-03-01_1.fastq.gz
[ ! -f  H021-E7U1N6_plasma-03-01_2.fastq.gz ] && ln -s AS-1383049-LR-78290_R2.fastq.gz H021-E7U1N6_plasma-03-01_2.fastq.gz
trim_galore \
    --fastqc   --clip_r1 10 --clip_r2 10 --three_prime_clip_r1 10 --three_prime_clip_r2 10 \
    --cores 8 \
    --paired \
    --gzip \
    H021-E7U1N6_plasma-03-01_1.fastq.gz \
    H021-E7U1N6_plasma-03-01_2.fastq.gz

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:TRIMGALORE":
    trimgalore: $(echo $(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*$//')
    cutadapt: $(cutadapt --version)
END_VERSIONS
