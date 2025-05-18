

set -ex



echo "exit [catch {package require -exact Tktable 2.10 }]" | xvfb-run tclsh
exit 0
