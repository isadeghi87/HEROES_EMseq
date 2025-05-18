#!/bin/bash
bedGraphToBigWig 2> /dev/null || [[ "$?" == 255 ]]


set -ex



which bedGraphToBigWig
exit 0
