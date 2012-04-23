#!/bin/bash
pushd win64
unzip ../packages/fftw-3.2.2-dll64.zip
lib /machine:X64 /def:libfftw3-3.def
lib /machine:X64 /def:libfftw3f-3.def
lib /machine:X64 /def:libfftw3l-3.def
cp -r ../src/blitz/blitz .
popd
