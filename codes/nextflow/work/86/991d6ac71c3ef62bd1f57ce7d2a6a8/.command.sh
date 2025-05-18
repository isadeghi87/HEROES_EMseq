#!/bin/bash -euo pipefail
preseq \
    lc_extrap \
     -verbose -bam \
    -pe \
    -output 5LB-017_plasma-01-01.lc_extrap.txt \
    5LB-017_plasma-01-01.sorted.bam
cp .command.err 5LB-017_plasma-01-01.command.log

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP":
    preseq: $(echo $(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*$//')
END_VERSIONS
