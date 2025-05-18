

set -ex



pip check
pytest -vv --color=yes --pyargs patsy --cov=patsy --cov-branch --cov-report=term-missing:skip-covered --cov-fail-under=47
exit 0
