

set -ex



echo "if {[catch {package require -exact BWidget 1.9.14; exit 0}]} {exit 1}" | xvfb-run tclsh
exit 0
