#!/bin/bash -euo pipefail
preseq \
    lc_extrap \
     -verbose -bam \
    -pe \
    -output 2LB-269_plasma-01-01.lc_extrap.txt \
    2LB-269_plasma-01-01.sorted.bam
cp .command.err 2LB-269_plasma-01-01.command.log

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP":
    preseq: $(echo $(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*$//')
END_VERSIONS
