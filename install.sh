#!/bin/sh

cd src
python3 -m numpy.f2py -c qSDHC.pyf qSDHC.for CPF3.for CPF.for --f90flags="-Ofast"
python3 -m numpy.f2py -c qSDV.pyf qSDV.for CPF3.for CPF.for --f90flags="-Ofast"
mkdir build
cd build
cmake ..
make
