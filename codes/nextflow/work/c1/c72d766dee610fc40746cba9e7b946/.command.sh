#!/bin/bash -euo pipefail
preseq \
    lc_extrap \
     -verbose -bam \
    -pe \
    -output H021-E7U1N6_plasma-03-01.lc_extrap.txt \
    H021-E7U1N6_plasma-03-01.sorted.bam
cp .command.err H021-E7U1N6_plasma-03-01.command.log

cat <<-END_VERSIONS > versions.yml
"NFCORE_METHYLSEQ:METHYLSEQ:PRESEQ_LCEXTRAP":
    preseq: $(echo $(preseq 2>&1) | sed 's/^.*Version: //; s/Usage:.*$//')
END_VERSIONS
