#!/bin/bash
eval "$(garden-exec)"
set -e

garden env-keep-only

# shared library dependencies
PYTHON3=desres-python/3.7.7-06c7
SCONS=scons/3.1.2-01c7
BOOST=boost/1.57.0-02c7
INCHI=inchi/1.05-01c7
SQLITE=sqlite/3.24.0-02c7
PYBIND11=pybind11/2.3.0-01c7

loadmodules() {
    garden load \
        $PYTHON3/bin \
        $SCONS/bin \
        $BOOST/lib \
        $INCHI/lib \
        $SQLITE/lib \
        $PYBIND11/lib \
        enscons/0.23.0-01c7/lib-python37 \
        auditwheel/3.1.1-01c7/bin \
        patchelf/0.9-01c7/bin \

    garden prepend-path DESRES_MODULE_CXXFLAGS -fpermissive
    garden prepend-path PYTHONPATH $(readlink -f $(dirname $0))    # for sconsutils
    export PYTHONVER=37
    export MSYS_WITH_INCHI=1
}

gendocs() {
    prefix=${PREFIX:-$PWD/build}
    (cd doc && BINDIR=$prefix/bin garden with -m $PYTHON3/bin make clean genhelp html )
    cp -r doc/build/html $PREFIX/doc/
}

loadmodules
BUILD_WHEEL=1 BUILD_WHEEL_VERSION=$PYTHONVER DESRES_LOCATION= scons "$@"

version=$(src/version.py)
auditwheel repair \
    --plat manylinux2014_x86_64 \
    -w build/wheel \
    build/wheel/msys-${version}-cp${PYTHONVER}-cp${PYTHONVER}m-linux_x86_64.whl

