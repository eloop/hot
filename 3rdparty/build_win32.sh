#!/bin/bash
pushd win32
unzip ../packages/fftw-3.2.2.pl1-dll32.zip
lib /def:libfftw3-3.def
lib /def:libfftw3f-3.def
lib /def:libfftw3l-3.def
cp -r ../src/blitz/blitz .
popd
