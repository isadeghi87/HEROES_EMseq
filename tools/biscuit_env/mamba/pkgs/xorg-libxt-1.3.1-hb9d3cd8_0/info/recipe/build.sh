#! /bin/bash
# Get an updated config.sub and config.guess
cp $BUILD_PREFIX/share/gnuconfig/config.* .

set -e
IFS=$' \t\n' # workaround for conda 4.2.13+toolchain bug

# Adopt a Unix-friendly path if we're on Windows (see bld.bat).
[ -n "$PATH_OVERRIDE" ] && export PATH="$PATH_OVERRIDE"

# On Windows we want $LIBRARY_PREFIX in both "mixed" (C:/Conda/...) and Unix
# (/c/Conda) forms, but Unix form is often "/" which can cause problems.
if [ -n "$LIBRARY_PREFIX_M" ] ; then
    mprefix="$LIBRARY_PREFIX_M"
    if [ "$LIBRARY_PREFIX_U" = / ] ; then
        uprefix=""
    else
        uprefix="$LIBRARY_PREFIX_U"
    fi
else
    mprefix="$PREFIX"
    uprefix="$PREFIX"
fi

# On Windows we need to regenerate the configure scripts.
if [ -n "$CYGWIN_PREFIX" ] ; then
    am_version=1.15 # keep sync'ed with meta.yaml
    export ACLOCAL=aclocal-$am_version
    export AUTOMAKE=automake-$am_version
    autoreconf_args=(
        --force
        --install
        -I "$mprefix/share/aclocal"
        -I "$BUILD_PREFIX_M/Library/usr/share/aclocal"
    )
    autoreconf "${autoreconf_args[@]}"

    # ./configure needs this
    export CPP=x86_64-w64-mingw32-cpp.exe

    # And we need to add the search path that lets libtool find the
    # msys2 stub libraries for ws2_32.
    platlibs=$(cd $(dirname $($CC --print-prog-name=ld))/../sysroot/usr/lib && pwd -W)
    test -f $platlibs/libws2_32.a || { exit "error locating libws2_32" ; exit 1 ; }
    export LDFLAGS="$LDFLAGS -L$platlibs"
else
    # for other platforms we just need to reconf to get the correct achitecture
    echo libtoolize
    libtoolize
    echo aclocal -I $PREFIX/share/aclocal -I $BUILD_PREFIX/share/aclocal
    aclocal -I $PREFIX/share/aclocal -I $BUILD_PREFIX/share/aclocal
    echo autoheader
    autoheader
    echo autoconf
    autoconf
    echo automake --force-missing --add-missing --include-deps
    automake --force-missing --add-missing --include-deps

    export CONFIG_FLAGS="--build=${BUILD}"
fi

# Workaround for osx-arm64 copied from xorg-libx11-feedstock -- the activation
# scripts should probably provide $CPP themselves.
if [[ "$(uname)" == "Darwin" ]]; then
    export CPP=clang-cpp
    ln -s $BUILD_PREFIX/bin/clang-cpp $BUILD_PREFIX/bin/cpp
fi

export PKG_CONFIG_LIBDIR=$uprefix/lib/pkgconfig:$uprefix/share/pkgconfig
configure_args=(
    --prefix=$mprefix
    --disable-static
    --disable-dependency-tracking
    --disable-selective-werror
    --disable-silent-rules
)

if [ -n "$CYGWIN_PREFIX" ] ; then
    # The default value for this argument messes up Libtool on Windows
    configure_args+=(--with-xfile-search-path="$uprefix/etc/X11/%L/%T/%N%C%S:$uprefix/share/X11/%L/%T/%N%C%S:$uprefix/lib/X11/%L/%T/%N%C%S")
fi

if [[ "${CONDA_BUILD_CROSS_COMPILATION}" == "1" ]] ; then
    configure_args+=(--enable-malloc0returnsnull)
fi

./configure "${configure_args[@]}"
make -j$CPU_COUNT
make install
if [[ "${CONDA_BUILD_CROSS_COMPILATION:-}" != "1" || "${CROSSCOMPILING_EMULATOR}" != "" ]]; then
make check
fi
rm -rf $uprefix/share/man $uprefix/share/doc/${PKG_NAME#xorg-}

# Remove any new Libtool files we may have installed. It is intended that
# conda-build will eventually do this automatically.
find $uprefix/. -name '*.la' -delete
