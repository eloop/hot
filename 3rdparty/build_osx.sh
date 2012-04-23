#!/bin/bash

IDIR=`pwd`/osx
echo "building in src/fftw-3.2.2"
pushd src/fftw-3.2.2
CC='gcc -m64' ./configure --with-pic --enable-threads --enable-optimize --enable-static --enable-single --enable-sse --prefix=$IDIR  
make clean
make install && make distclean
popd

echo "installing to $IDIR"
echo "building in src/blitz"
pushd src/blitz
make distclean
LDFLAGS="-m64" CXX="c++ -m64" ./configure  --enable-threadsafe --enable-optimize --enable-static --prefix=$IDIR 
# we don't want to build the docs 
make clean
(cd lib && make install)
(cd blitz && make install)
(cd random && make install)
make distclean
popd

