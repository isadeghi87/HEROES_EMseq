

set -ex



pip check
typer --version
typer --help
cd typer && coverage run --source=typer --branch -m pytest -vv --color=yes --tb=long -k "not ((multiple_values and main) or completion)"
coverage report --show-missing --skip-covered --fail-under=66
exit 0
