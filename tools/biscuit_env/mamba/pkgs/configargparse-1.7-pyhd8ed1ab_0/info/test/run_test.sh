

set -ex



pip check
pytest -vv --cov=configargparse --cov-branch --cov-report=term-missing:skip-covered --no-cov-on-fail --cov-fail-under=74
exit 0
