if [ -z ${CONDA_BUILD+x} ]; then
    source /opt/conda/conda-bld/bismark_1698247648149/work/build_env_setup.sh
fi
#!/bin/bash

mkdir -p $PREFIX/bin/
rm -rf Docs
cp -r * $PREFIX/bin/
