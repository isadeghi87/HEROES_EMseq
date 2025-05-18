

set -ex



pip check
python -c "from contourpy.util import build_config; from pprint import pprint; pprint(build_config())"
python -c "import platform, sys; print(sys.version, platform.uname())"
python -m pytest -v
exit 0
