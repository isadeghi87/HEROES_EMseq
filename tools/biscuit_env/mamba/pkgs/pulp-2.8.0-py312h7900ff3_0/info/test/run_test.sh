

set -ex



pip check
python -c "import pulp; assert pulp.LpSolverDefault, 'default solver not available'"
exit 0
