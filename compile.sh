#!/bin/bash

# exit immediately on error
set -e

UNAME=`uname`
HIH="${HOME}/houdini${HOUDINI_MAJOR_RELEASE}.${HOUDINI_MINOR_RELEASE}"

if [ ${UNAME} = 'Darwin' ] ; then
    # pushd 3rdparty 
    # ./build_osx.sh
    # popd
    IFLAGS="-I. -I./3rdparty/include -I ./3rdparty/osx/include"
    LFLAGS="-L ./3rdparty/osx/lib  -l fftw3f -l blitz"
    FLAGS="$IFLAGS $LFLAGS -i ${HIH}/dso"
elif [ ${UNAME} = 'Linux' ] ; then
    pushd 3rdparty 
    ./build_linux.sh
    popd
    IFLAGS="-I. -I./3rdparty/include -I ./3rdparty/linux/include"
    LFLAGS="-L ./3rdparty/linux/lib  -l fftw3f -l blitz"
    FLAGS="$IFLAGS $LFLAGS -i ${HIH}/dso"
    hcustom  -e $FLAGS SOP_Ocean.C 
    hcustom  -e $FLAGS VEX_Ocean.C 
    hcustom  SOP_Cleave.C
else
    echo "Unknown architecture, sorry."
    exit 1
fi

hcustom  -e $FLAGS SOP_Ocean.C 
hcustom  -e $FLAGS VEX_Ocean.C 
hcustom  SOP_Cleave.C

VEXDSO="${HIH}/vex"
mkdir -p ${VEXDSO}
cp VEXdso_${UNAME} ${VEXDSO}/VEXdso
