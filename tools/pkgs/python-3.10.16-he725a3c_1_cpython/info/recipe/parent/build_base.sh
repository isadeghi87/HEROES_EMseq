#!/bin/bash
set -ex

# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/libtool/build-aux/config.* .

# The LTO/PGO information was sourced from @pitrou and the Debian rules file in:
# http://http.debian.net/debian/pool/main/p/python3.6/python3.6_3.6.2-2.debian.tar.xz
# https://packages.debian.org/source/sid/python3.6
# or:
# http://bazaar.launchpad.net/~doko/python/pkg3.5-debian/view/head:/rules#L255
# .. but upstream regrtest.py now has --pgo (since >= 3.6) and skips tests that are:
# "not helpful for PGO".

VERFULL=${PKG_VERSION}
VER=${PKG_VERSION%.*}
VERNODOTS=${VER//./}
TCLTK_VER=${tk}
# Disables some PGO/LTO
QUICK_BUILD=no

_buildd_static=build-static
_buildd_shared=build-shared
_ENABLE_SHARED=--enable-shared
# We *still* build a shared lib here for non-static embedded use cases
_DISABLE_SHARED=--disable-shared
# Hack to allow easily comparing static vs shared interpreter performance
# .. hack because we just build it shared in both the build-static and
# build-shared directories.
# Yes this hack is a bit confusing, sorry about that.
if [[ ${PY_INTERP_LINKAGE_NATURE} == shared ]]; then
  _DISABLE_SHARED=--enable-shared
  _ENABLE_SHARED=--enable-shared
fi

# For debugging builds, set this to no to disable profile-guided optimization
if [[ ${DEBUG_C} == yes ]]; then
  _OPTIMIZED=no
else
  _OPTIMIZED=yes
fi

# Since these take very long to build in our emulated ci, disable for now
if [[ ${target_platform} == linux-aarch64 ]]; then
  _OPTIMIZED=no
fi
if [[ ${target_platform} == linux-ppc64le ]]; then
  _OPTIMIZED=no
fi

declare -a _dbg_opts
if [[ ${DEBUG_PY} == yes ]]; then
  # This Python will not be usable with non-debug Python modules.
  _dbg_opts+=(--with-pydebug)
  DBG=d
else
  DBG=
fi

ABIFLAGS=${DBG}
VERABI=${VER}${DBG}

# Make sure the "python" value in conda_build_config.yaml is up to date.
test "${PY_VER}" = "${VER}"

# This is the mechanism by which we fall back to default gcc, but having it defined here
# would probably break the build by using incorrect settings and/or importing files that
# do not yet exist.
unset _PYTHON_SYSCONFIGDATA_NAME
unset _CONDA_PYTHON_SYSCONFIGDATA_NAME

# Prevent lib/python${VER}/_sysconfigdata_*.py from ending up with full paths to these things
# in _build_env because _build_env will not get found during prefix replacement, only _h_env_placeh ...
AR=$(basename "${AR}")

# CC must contain the string 'gcc' or else distutils thinks it is on macOS and uses '-R' to set rpaths.
if [[ ${target_platform} == osx-* ]]; then
  CC=$(basename "${CC}")
else
  CC=$(basename "${GCC}")
fi
CXX=$(basename "${CXX}")
RANLIB=$(basename "${RANLIB}")
READELF=$(basename "${READELF}")

if [[ ${HOST} =~ .*darwin.* ]] && [[ -n ${CONDA_BUILD_SYSROOT} ]]; then
  # Python's setup.py will figure out that this is a macOS sysroot.
  CFLAGS="-isysroot ${CONDA_BUILD_SYSROOT} "${CFLAGS}
  LDFLAGS="-isysroot ${CONDA_BUILD_SYSROOT} "${LDFLAGS}
  CPPFLAGS="-isysroot ${CONDA_BUILD_SYSROOT} "${CPPFLAGS}
fi

# Debian uses -O3 then resets it at the end to -O2 in _sysconfigdata.py
if [[ ${_OPTIMIZED} = yes ]]; then
  CPPFLAGS=$(echo "${CPPFLAGS}" | sed "s/-O2/-O3/g")
  CFLAGS=$(echo "${CFLAGS}" | sed "s/-O2/-O3/g")
  CXXFLAGS=$(echo "${CXXFLAGS}" | sed "s/-O2/-O3/g")
fi

if [[ ${CONDA_FORGE} == yes ]]; then
  ${SYS_PYTHON} ${RECIPE_DIR}/brand_python.py
fi

if [[ "$target_platform" == linux-* ]]; then
  cp ${PREFIX}/include/uuid/uuid.h ${PREFIX}/include/uuid.h
fi

declare -a LTO_CFLAGS=()

# Following is needed for building extensions like zlib
CPPFLAGS=${CPPFLAGS}" -I${PREFIX}/include"

re='^(.*)(-I[^ ]*)(.*)$'
if [[ ${CFLAGS} =~ $re ]]; then
  CFLAGS="${BASH_REMATCH[1]}${BASH_REMATCH[3]}"
fi

# Force rebuild to avoid:
# ../work/Modules/unicodename_db.h:24118:30: note: (near initialization for 'code_hash')
# ../work/Modules/unicodename_db.h:24118:33: warning: excess elements in scalar initializer
#      0, 0, 12018, 0, 0, 0, 0, 0, 4422, 4708, 3799, 119358, 119357, 0, 120510,
#                                  ^~~~
# This should have been fixed by https://github.com/python/cpython/commit/7c69c1c0fba8c1c8ff3969bce4c1135736a4cc58
# .. but that appears incomplete. In particular, the generated files contain:
# /* this file was generated by Tools/unicode/makeunicodedata.py 3.2 */
# .. yet the PR updated to version of makeunicodedata.py to 3.3
# rm -f Modules/unicodedata_db.h Modules/unicodename_db.h
# ${SYS_PYTHON} ${SRC_DIR}/Tools/unicode/makeunicodedata.py
# .. instead we revert this commit for now.

export CPPFLAGS CFLAGS CXXFLAGS LDFLAGS

if [[ ${target_platform} == osx-* ]]; then
  sed -i -e "s/@OSX_ARCH@/$ARCH/g" Lib/distutils/unixccompiler.py
fi

if [[ "${CONDA_BUILD_CROSS_COMPILATION}" == "1" ]]; then
  # Build the exact same Python for the build machine. It would be nice (and might be
  # possible already?) to be able to make this just an 'exact' pinned build dependency
  # of a split-package?
  BUILD_PYTHON_PREFIX=${PWD}/build-python-install
  mkdir build-python-build
  pushd build-python-build
    (unset CPPFLAGS LDFLAGS;
     export CC=${CC_FOR_BUILD} \
            CXX=${CXX_FOR_BUILD} \
            CPP="${CC_FOR_BUILD} -E" \
            CFLAGS="-O2" \
            AR="$(${CC_FOR_BUILD} --print-prog-name=ar)" \
            RANLIB="$(${CC_FOR_BUILD} --print-prog-name=ranlib)" \
            LD="$(${CC_FOR_BUILD} --print-prog-name=ld)" && \
      ${SRC_DIR}/configure --build=${BUILD} \
                           --host=${BUILD} \
                           --prefix=${BUILD_PYTHON_PREFIX} \
                           --with-ensurepip=no \
                           --with-tzpath=${PREFIX}/share/zoneinfo \
                           --with-platlibdir=lib && \
      make -j${CPU_COUNT} && \
      make install)
    export PATH=${BUILD_PYTHON_PREFIX}/bin:${PATH}
    ln -s ${BUILD_PYTHON_PREFIX}/bin/python${VER} ${BUILD_PYTHON_PREFIX}/bin/python
  popd
  echo "ac_cv_file__dev_ptmx=yes"        > config.site
  echo "ac_cv_file__dev_ptc=yes"        >> config.site
  echo "ac_cv_pthread=yes"              >> config.site
  echo "ac_cv_little_endian_double=yes" >> config.site
  echo "ac_cv_aligned_required=no" >> config.site
  if [[ ${target_platform} == osx-arm64 ]]; then
      echo "ac_cv_file__dev_ptc=no" >> config.site
      echo "ac_cv_pthread_is_default=yes" >> config.site
      echo "ac_cv_working_tzset=yes" >> config.site
      echo "ac_cv_pthread_system_supported=yes" >> config.site
  fi
  export CONFIG_SITE=${PWD}/config.site
  # This is needed for libffi:
  export PKG_CONFIG_PATH=${PREFIX}/lib/pkgconfig
fi

# This causes setup.py to query the sysroot directories from the compiler, something which
# IMHO should be done by default anyway with a flag to disable it to workaround broken ones.
# Technically, setting _PYTHON_HOST_PLATFORM causes setup.py to consider it cross_compiling
if [[ -n ${HOST} ]]; then
  if [[ ${HOST} =~ .*darwin.* ]]; then
    # Even if BUILD is .*darwin.* you get better isolation by cross_compiling (no /usr/local)
    IFS='-' read -r host_arch host_os host_kernel <<<"${HOST}"
    export _PYTHON_HOST_PLATFORM=darwin-${host_arch}
  else
    IFS='-' read -r host_arch host_vendor host_os host_libc <<<"${HOST}"
    export _PYTHON_HOST_PLATFORM=${host_os}-${host_arch}
  fi
fi

if [[ ${target_platform} == osx-64 ]]; then
  export MACHDEP=darwin
  export ac_sys_system=Darwin
  export ac_sys_release=13.4.0
  export MACOSX_DEFAULT_ARCH=x86_64
  # TODO: check with LLVM 12 if the following hack is needed.
  # https://reviews.llvm.org/D76461 may have fixed the need for the following hack.
  echo '#!/bin/bash' > $BUILD_PREFIX/bin/$HOST-llvm-ar
  echo "$BUILD_PREFIX/bin/llvm-ar --format=darwin" '"$@"' >> $BUILD_PREFIX/bin/$HOST-llvm-ar
  chmod +x $BUILD_PREFIX/bin/$HOST-llvm-ar
  export ARCHFLAGS="-arch x86_64"
elif [[ ${target_platform} == osx-arm64 ]]; then
  export MACHDEP=darwin
  export ac_sys_system=Darwin
  export ac_sys_release=20.0.0
  export MACOSX_DEFAULT_ARCH=arm64
  echo '#!/bin/bash' > $BUILD_PREFIX/bin/$HOST-llvm-ar
  echo "$BUILD_PREFIX/bin/llvm-ar --format=darwin" '"$@"' >> $BUILD_PREFIX/bin/$HOST-llvm-ar
  chmod +x $BUILD_PREFIX/bin/$HOST-llvm-ar
  export ARCHFLAGS="-arch arm64"
  export CFLAGS="$CFLAGS $ARCHFLAGS"
elif [[ ${target_platform} == linux-* ]]; then
  export MACHDEP=linux
  export ac_sys_system=Linux
  export ac_sys_release=
fi

# Not used at present but we should run 'make test' and finish up TESTOPTS (see debians rules).
declare -a TEST_EXCLUDES
TEST_EXCLUDES+=(test_ensurepip test_venv)
TEST_EXCLUDES+=(test_tcl test_codecmaps_cn test_codecmaps_hk
                test_codecmaps_jp test_codecmaps_kr test_codecmaps_tw
                test_normalization test_ossaudiodev test_socket)
if [[ ! -f /dev/dsp ]]; then
  TEST_EXCLUDES+=(test_linuxaudiodev test_ossaudiodev)
fi
# hangs on Aarch64, see LP: #1264354
if [[ ${CC} =~ .*-aarch64.* ]]; then
  TEST_EXCLUDES+=(test_faulthandler)
fi
if [[ ${CC} =~ .*-arm.* ]]; then
  TEST_EXCLUDES+=(test_ctypes)
  TEST_EXCLUDES+=(test_compiler)
fi

declare -a _common_configure_args
_common_configure_args+=(--prefix=${PREFIX})
_common_configure_args+=(--build=${BUILD})
_common_configure_args+=(--host=${HOST})
_common_configure_args+=(--enable-ipv6)
_common_configure_args+=(--with-ensurepip=no)
_common_configure_args+=(--with-tzpath=${PREFIX}/share/zoneinfo)
_common_configure_args+=(--with-computed-gotos)
_common_configure_args+=(--with-system-ffi)
_common_configure_args+=(--enable-loadable-sqlite-extensions)
_common_configure_args+=(--with-tcltk-includes="-I${PREFIX}/include")
_common_configure_args+=("--with-tcltk-libs=-L${PREFIX}/lib -ltcl8.6 -ltk8.6")
_common_configure_args+=(--with-platlibdir=lib)

# Add more optimization flags for the static Python interpreter:
declare -a PROFILE_TASK=()
if [[ ${_OPTIMIZED} == yes ]]; then
  _common_configure_args+=(--with-lto)
  if [[ "$CONDA_BUILD_CROSS_COMPILATION" != "1" ]]; then
    _common_configure_args+=(--enable-optimizations)
    _MAKE_TARGET=profile-opt
    # To speed up build times during testing (1):
    if [[ ${QUICK_BUILD} == yes ]]; then
      # TODO :: It seems this is just profiling everything, on Windows, only 40 odd tests are
      #         run while on Unix, all 400+ are run, making this slower and less well curated
      _PROFILE_TASK+=(PROFILE_TASK="-m test --pgo")
    else
      # From talking to Steve Dower, who implemented pgo/pgo-extended, it is really not worth
      # it to run pgo-extended (which runs the whole test-suite). The --pgo set of tests are
      # curated specifically to be useful/appropriate for pgo instrumentation.
      # _PROFILE_TASK+=(PROFILE_TASK="-m test --pgo-extended")
      _PROFILE_TASK+=(PROFILE_TASK="-m test --pgo")
    fi
  fi
  if [[ ${CC} =~ .*gcc.* ]]; then
    LTO_CFLAGS+=(-fuse-linker-plugin)
    LTO_CFLAGS+=(-ffat-lto-objects)
    # -flto must come after -flto-partition due to the replacement code
    # TODO :: Replace the replacement code using conda-build's in-build regex replacement.
    LTO_CFLAGS+=(-flto-partition=none)
    LTO_CFLAGS+=(-flto)
  else
    # TODO :: Check if -flto=thin gives better results. It is about faster
    #         compilation rather than faster execution so probably not:
    # http://clang.llvm.org/docs/ThinLTO.html
    # http://blog.llvm.org/2016/06/thinlto-scalable-and-incremental-lto.html
    LTO_CFLAGS+=(-flto)
    # -flto breaks the check to determine whether float word ordering is bigendian
    # see:
    # https://bugs.python.org/issue28015
    # https://bugs.python.org/issue38527
    # manually specify this setting
    export ax_cv_c_float_words_bigendian=no
  fi
  export CFLAGS="${CFLAGS} ${LTO_CFLAGS[@]}"
else
  _MAKE_TARGET=
fi

mkdir -p ${_buildd_shared}
pushd ${_buildd_shared}
  ${SRC_DIR}/configure "${_common_configure_args[@]}" \
                       "${_dbg_opts[@]}" \
                       --oldincludedir=${BUILD_PREFIX}/${HOST}/sysroot/usr/include \
                       --enable-shared
popd

mkdir -p ${_buildd_static}
pushd ${_buildd_static}
  ${SRC_DIR}/configure "${_common_configure_args[@]}" \
                       "${_dbg_opts[@]}" \
                       -oldincludedir=${BUILD_PREFIX}/${HOST}/sysroot/usr/include \
                       ${_DISABLE_SHARED} "${_PROFILE_TASK[@]}"
popd

if [[ "${CI}" == "travis" ]]; then
  # Travis has issues with long logs
  make -j${CPU_COUNT} -C ${_buildd_static} \
       EXTRA_CFLAGS="${EXTRA_CFLAGS}" \
       ${_MAKE_TARGET} "${_PROFILE_TASK[@]}" 2>&1 >make-static.log
else
  make -j${CPU_COUNT} -C ${_buildd_static} \
       EXTRA_CFLAGS="${EXTRA_CFLAGS}" \
       ${_MAKE_TARGET} "${_PROFILE_TASK[@]}" 2>&1 | tee make-static.log
fi
if rg "Failed to build these modules" make-static.log; then
  echo "(static) :: Failed to build some modules, check the log"
  exit 1
fi

if [[ "${CI}" == "travis" ]]; then
  # Travis has issues with long logs
  make -j${CPU_COUNT} -C ${_buildd_shared} \
          EXTRA_CFLAGS="${EXTRA_CFLAGS}" 2>&1 >make-shared.log
else
  make -j${CPU_COUNT} -C ${_buildd_shared} \
          EXTRA_CFLAGS="${EXTRA_CFLAGS}" 2>&1 | tee make-shared.log
fi
if rg "Failed to build these modules" make-shared.log; then
  echo "(shared) :: Failed to build some modules, check the log"
  exit 1
fi

# build a static library with PIC objects and without LTO/PGO
make -j${CPU_COUNT} -C ${_buildd_shared} \
        EXTRA_CFLAGS="${EXTRA_CFLAGS}" \
        LIBRARY=libpython${VERABI}-pic.a libpython${VERABI}-pic.a

make -C ${_buildd_static} install

declare -a _FLAGS_REPLACE=()
if [[ ${_OPTIMIZED} == yes ]]; then
  _FLAGS_REPLACE+=(-O3)
  _FLAGS_REPLACE+=(-O2)
  _FLAGS_REPLACE+=("-fprofile-use")
  _FLAGS_REPLACE+=("")
  _FLAGS_REPLACE+=("-fprofile-correction")
  _FLAGS_REPLACE+=("")
  _FLAGS_REPLACE+=("-L.")
  _FLAGS_REPLACE+=("")
  for _LTO_CFLAG in "${LTO_CFLAGS[@]}"; do
    _FLAGS_REPLACE+=(${_LTO_CFLAG})
    _FLAGS_REPLACE+=("")
  done
fi
# Install the shared library (for people who embed Python only, e.g. GDB).
# Linking module extensions to this on Linux is redundant (but harmless).
# Linking module extensions to this on Darwin is harmful (multiply defined symbols).
cp -pf ${_buildd_shared}/libpython*${SHLIB_EXT}* ${PREFIX}/lib/
if [[ ${target_platform} =~ .*linux.* ]]; then
  ln -sf ${PREFIX}/lib/libpython${VERABI}${SHLIB_EXT}.1.0 ${PREFIX}/lib/libpython${VERABI}${SHLIB_EXT}
fi

SYSCONFIG=$(find ${_buildd_static}/$(cat ${_buildd_static}/pybuilddir.txt) -name "_sysconfigdata*.py" -print0)
cat ${SYSCONFIG} | ${SYS_PYTHON} "${RECIPE_DIR}"/replace-word-pairs.py \
  "${_FLAGS_REPLACE[@]}"  \
    > ${PREFIX}/lib/python${VER}/$(basename ${SYSCONFIG})
MAKEFILE=$(find ${PREFIX}/lib/python${VER}/ -path "*config-*/Makefile" -print0)
cp ${MAKEFILE} /tmp/Makefile-$$
cat /tmp/Makefile-$$ | ${SYS_PYTHON} "${RECIPE_DIR}"/replace-word-pairs.py \
  "${_FLAGS_REPLACE[@]}"  \
    > ${MAKEFILE}
# Check to see that our differences took.
# echo diff -urN ${SYSCONFIG} ${PREFIX}/lib/python${VER}/$(basename ${SYSCONFIG})
# diff -urN ${SYSCONFIG} ${PREFIX}/lib/python${VER}/$(basename ${SYSCONFIG})

# Python installs python${VER}m and python${VER}, one as a hardlink to the other. conda-build breaks these
# by copying. Since the executable may be static it may be very large so change one to be a symlink
# of the other. In this case, python${VER}m will be the symlink.
if [[ -f ${PREFIX}/bin/python${VER}m ]]; then
  rm -f ${PREFIX}/bin/python${VER}m
  ln -s ${PREFIX}/bin/python${VER} ${PREFIX}/bin/python${VER}m
fi
ln -s ${PREFIX}/bin/python${VER} ${PREFIX}/bin/python
ln -s ${PREFIX}/bin/pydoc${VER} ${PREFIX}/bin/pydoc
# Workaround for https://github.com/conda/conda/issues/10969
ln -s ${PREFIX}/bin/python3.10 ${PREFIX}/bin/python3.1

# Remove test data to save space
# Though keep `support` as some things use that.
# TODO :: Make a subpackage for this once we implement multi-level testing.
pushd ${PREFIX}/lib/python${VER}
  mkdir test_keep
  mv test/__init__.py test/support test/test_support* test/test_script_helper* test_keep/
  rm -rf test */test
  mv test_keep test
popd

# Size reductions:
pushd ${PREFIX}
  if [[ -f lib/libpython${VERABI}.a ]]; then
    chmod +w lib/libpython${VERABI}.a
    ${STRIP} -S lib/libpython${VERABI}.a
  fi
  CONFIG_LIBPYTHON=$(find lib/python${VER}/config-${VERABI}* -name "libpython${VERABI}.a")
  if [[ -f lib/libpython${VERABI}.a ]] && [[ -f ${CONFIG_LIBPYTHON} ]]; then
    chmod +w ${CONFIG_LIBPYTHON}
    rm ${CONFIG_LIBPYTHON}
  fi
popd

# OLD_HOST is with CentOS version in them. When building this recipe
# with the compilers from conda-forge OLD_HOST != HOST, but when building
# with the compilers from defaults OLD_HOST == HOST. Both cases are handled in the
# code below
case "$target_platform" in
  linux-64)
    OLD_HOST=$(echo ${HOST} | sed -e 's/-conda-/-conda_cos6-/g')
    ;;
  linux-*)
    OLD_HOST=$(echo ${HOST} | sed -e 's/-conda-/-conda_cos7-/g')
    ;;
  *)
    OLD_HOST=$HOST
    ;;
esac

# Copy sysconfig that gets recorded to a non-default name
# using the new compilers with python will require setting _PYTHON_SYSCONFIGDATA_NAME
# to the name of this file (minus the .py extension)
pushd "${PREFIX}"/lib/python${VER}
  # On Python 3.5 _sysconfigdata.py was getting copied in here and compiled for some reason.
  # This breaks our attempt to find the right one as recorded_name.
  find lib-dynload -name "_sysconfigdata*.py*" -exec rm {} \;
  recorded_name=$(find . -name "_sysconfigdata*.py")
  our_compilers_name=_sysconfigdata_$(echo ${HOST} | sed -e 's/[.-]/_/g').py
  # So we can see if anything has significantly diverged by looking in a built package.
  cp ${recorded_name} ${recorded_name}.orig
  cp ${recorded_name} sysconfigfile
  # fdebug-prefix-map for python work dir is useless for extensions
  sed -i.bak "s@-fdebug-prefix-map=$SRC_DIR=/usr/local/src/conda/python-$PKG_VERSION@@g" sysconfigfile
  sed -i.bak "s@-fdebug-prefix-map=$PREFIX=/usr/local/src/conda-prefix@@g" sysconfigfile
  # Append the conda-forge zoneinfo to the end
  sed -i.bak "s@zoneinfo'@zoneinfo:$PREFIX/share/tzinfo'@g" sysconfigfile
  # Remove osx sysroot as it depends on the build machine
  # be sure CONDA_BUILD_SYSROOT has value, as other we will remove here instead spaces
  if [[ "${target_platform}" == osx-* ]] && [[ -n ${CONDA_BUILD_SYSROOT} ]]; then
    sed -i.bak "s@-isysroot @@g" sysconfigfile
    sed -i.bak "s@$CONDA_BUILD_SYSROOT @@g" sysconfigfile
  fi
  # Remove unfilled config option
  sed -i.bak "s/@SGI_ABI@//g" sysconfigfile
  sed -i.bak "s@$BUILD_PREFIX/bin/${HOST}-llvm-ar@${HOST}-ar@g" sysconfigfile
  # Remove GNULD=yes to make sure new-dtags are not used
  sed -i.bak "s/'GNULD': 'yes'/'GNULD': 'no'/g" sysconfigfile
  cp sysconfigfile ${our_compilers_name}

  sed -i.bak "s@${HOST}@${OLD_HOST}@g" sysconfigfile
  old_compiler_name=_sysconfigdata_$(echo ${OLD_HOST} | sed -e 's/[.-]/_/g').py
  cp sysconfigfile ${old_compiler_name}

  # For system gcc remove the triple
  sed -i.bak "s@$OLD_HOST-c++@g++@g" sysconfigfile
  sed -i.bak "s@$OLD_HOST-@@g" sysconfigfile
  if [[ "$target_platform" == linux* ]]; then
    # For linux, make sure the system gcc uses our linker
    sed -i.bak "s@-pthread@-pthread -B $PREFIX/compiler_compat@g" sysconfigfile
  fi
  # Don't set -march and -mtune for system gcc
  sed -i.bak "s@-march=[^( |\\\"|\\\')]*@@g" sysconfigfile
  sed -i.bak "s@-mtune=[^( |\\\"|\\\')]*@@g" sysconfigfile
  # Remove these flags that older compilers and linkers may not know
  for flag in "-fstack-protector-strong" "-ffunction-sections" "-pipe" "-fno-plt" \
            "-ftree-vectorize" "-Wl,--sort-common" "-Wl,--as-needed" "-Wl,-z,relro" \
            "-Wl,-z,now" "-Wl,--disable-new-dtags" "-Wl,--gc-sections" "-Wl,-O2" \
            "-fPIE" "-ftree-vectorize" "-mssse3" "-Wl,-pie" "-Wl,-dead_strip_dylibs" \
            "-Wl,-headerpad_max_install_names"; do
    sed -i.bak "s@$flag@@g" sysconfigfile
  done
  # Cleanup some extra spaces from above
  sed -i.bak "s@' [ ]*@'@g" sysconfigfile
  cp sysconfigfile $recorded_name
  echo "========================sysconfig==========================="
  cat $recorded_name
  echo "============================================================"

  rm sysconfigfile
  rm sysconfigfile.bak
popd

if [[ ${HOST} =~ .*linux.* ]]; then
  mkdir -p ${PREFIX}/compiler_compat
  ln -s ${PREFIX}/bin/${HOST}-ld ${PREFIX}/compiler_compat/ld
  echo "Files in this folder are to enhance backwards compatibility of anaconda software with older compilers."   > ${PREFIX}/compiler_compat/README
  echo "See: https://github.com/conda/conda/issues/6030 for more information."                                   >> ${PREFIX}/compiler_compat/README
fi

# There are some strange distutils files around. Delete them
rm -rf ${PREFIX}/lib/python${VER}/distutils/command/*.exe

python -c "import compileall,os;compileall.compile_dir(os.environ['PREFIX'])"
rm ${PREFIX}/lib/libpython${VER}.a
if [[ "$target_platform" == linux-* ]]; then
  rm ${PREFIX}/include/uuid.h
fi

# Workaround for old conda versions which fail to install noarch packages for Python 3.10+
# https://github.com/conda/conda/issues/10969
ln -s "${PREFIX}/lib/python3.10" "${PREFIX}/lib/python3.1"
