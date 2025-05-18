

set -ex



biscuit 2>&1 | grep "Usage"
QC.sh -h 2>&1 | grep "Usage"
build_biscuit_QC_assets.pl -h 2>&1 | grep "Usage"
exit 0
