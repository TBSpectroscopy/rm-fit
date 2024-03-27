# RM-Fit

What is RM-Fit?
---------------

RM-Fit (Relaxation Matrix-Fit) is a multi-spectrum fitting python script using the relaxation matrix formalism (C++ and Fortran). It is capable of modeling full line mixing in combination with speed-dependent and narrowing effects. The current version does not allow for first-order line mixing, but it will be added eventually.

Compilation
-----------
### Dependencies

- gcc/g++/gfortran
- CMake
- [pybind11](https://github.com/pybind/pybind11)
- numpy
- [Eigen](https://eigen.tuxfamily.org)

### Instructions

1. Clone source code
    ```
    git clone https://github.com/TBSpectroscopy/rm-fit.git
    ```

2. Access the source code directory
    ```
    cd rm-fit/src
    ```

3. Locate pybind11 in your system (depends on installation method) and replace it in `CMakelist.txt` line 4.

    `set(pybind11_DIR /usr/include/pybind11)` -> `set(pybind11_DIR [your_pybind11_location])`

4. Compile the Fortran and C++ libraries
    ```
    python -m numpy.f2py -c qSDHC.pyf qSDHC.for CPF3.for CPF.for --f90flags="-Ofast"
    python -m numpy.f2py -c qSDV.pyf qSDV.for CPF3.for CPF.for --f90flags="-Ofast"
    mkdir build
    cd build
    cmake ../
    make
    ```

Usage
-----

Open a terminal in the main folder `rm-fit/`

```
python rm-fit.py [option] [input_file]
```

| `option` | Description |
|---|---|
| -h, --help | Display help text|
| -f | Fit |
| -c | Calculate |

### Dependencies

- python 3
- numpy
- scipy

### Input

An example is provided in `data/input/` named `ch4_co2.txt`. The input file contains all the relevant information to calculate spectra. It is separated into blocks marked by a '$' symbol followed by a name. What follows the name on the same line is ignored.

`$SPECTRUM` contains all the information specific to one spectrum. Can be as many as needed

`$LINELIST` contains the inputs for a species/isotopologue line list (diagonals and off-diagonals). Can be as many as needed

`$CALCULATION` contains the information used for all spectra and the output path. Unique

 Initial parameters are followed by an 'f' if you want to fit them, or 'c' if you want them to stay constant. The parameters inside the blocks are detailed in the example file.
