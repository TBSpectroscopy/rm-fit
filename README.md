# RM-Fit

What is RM-Fit?
---------------

RM-Fit (Relaxation Matrix-Fit) is a free (libre) multi-spectrum fitting python script using the relaxation matrix formalism (C++ and Fortran). It is capable of modeling full line mixing in combination with speed-dependent and narrowing effects, as well as first-order line-mixing.

Compilation
-----------
### Dependencies

- [gcc, g++](https://gcc.gnu.org/) / [clang](https://clang.llvm.org/) with C++ >= 14
- [gfortran](https://gcc.gnu.org/fortran/)
- [CMake](https://cmake.org/)
- [Python 3.x](https://www.python.org/)
- [pybind11](https://github.com/pybind/pybind11)
- [meson](https://mesonbuild.com/)
- [NumPy](https://numpy.org/)
- [Eigen](https://eigen.tuxfamily.org)

#### Linux
Arch Linux: Every dependencies are available through the package manager and the AUR

#### macOS
Clang is available through Xcode Command Line Tools (`sudo xcode-select --install`). Everything else is available through [Homebrew](https://brew.sh/) and [pip](https://pip.pypa.io/)

### Instructions

1. Clone source code
    ```
    git clone https://github.com/TBSpectroscopy/rm-fit.git
    ```

2. Access the source code directory
    ```
    cd rm-fit/src
    ```

3. Compile the Fortran and C++ libraries
    ```
    python3 -m numpy.f2py -c qSDHC.pyf qSDHC.for CPF3.for CPF.for --f90flags="-Ofast"
    python3 -m numpy.f2py -c qSDV.pyf qSDV.for CPF3.for CPF.for --f90flags="-Ofast"
    mkdir build
    cd build
    cmake ..
    make
    ```

Usage
-----

Open a terminal in the main folder `rm-fit/`

```
python3 rm-fit.py [option] [input_file]
```

| `option` | Description |
|---|---|
| -h, --help | Display help text|
| -f | Fit |
| -c | Calculate |

### Dependencies

- Python 3.x
- NumPy
- SciPy

### Input

An example is provided in `data/input/` named `ch4_co2.txt`. The input file contains all the relevant information to calculate spectra. It is separated into blocks marked by a '$' symbol followed by a name. What follows the name on the same line is ignored. All file/folder paths inside must be absolute or relative to the location of the input file.

`$SPECTRUM` contains all the information specific to one spectrum. Can be as many as needed. Note that this repository does not provide the spectrum to run the example input. [However, we do provide it as a download here.](https://owncloud.ulb.ac.be/index.php/s/yLinSsQ1t480Qxc/download)

`$LINELIST` contains the inputs for a species/isotopologue line list (diagonals and off-diagonals). Can be as many as needed

`$CALCULATION` contains the information used for all spectra and the output path. Unique

Initial parameters are followed by an 'f' if you want to fit them, or 'c' if you want them to stay constant. The parameters inside the blocks are detailed in the example file.

### Linelist

Lines needed to compute the example input file are provided by the linelists `diag_61.txt` and `diag_62.txt` and the off-diagonal list `offdiag_61.txt`. All files must contain a header giving details about the species, and a format section inside quotation marks allowing RM-Fit to know how to read the parameters. The format follows a structure similar to python's string `format` method, with braces containing the name of a parameter and its width separated by a colon. The width can include a letter indicating how the numbers should be formatted. Spaces are considered inside the format section and can be used to ignore sections of the lists (ignored parameters are assumed to be 0).


Citation
--------

If you use RM-Fit in your work, please cite it using the following references:

T. Bertin and J. Vander Auwera, “CO2 collision-induced line parameters for the ν3 band of 12CH4 measured using a hard-collision speed-dependent line shape and the relaxation matrix formalism,” Journal of Quantitative Spectroscopy and Radiative Transfer, vol. 324, p. 109069, Sep. 2024, doi: 10.1016/j.jqsrt.2024.109069. [Link to article](https://www.sciencedirect.com/science/article/abs/pii/S0022407324001766)

H. Tran, N. H. Ngo, and J. Hartmann, “Efficient computation of some speed-dependent isolated line profiles,” Journal of Quantitative Spectroscopy and Radiative Transfer, vol. 129, pp. 199–203, Nov. 2013, doi: 10.1016/j.jqsrt.2013.06.015. [Link to article](https://www.sciencedirect.com/science/article/abs/pii/S0022407313002598)

F. Schreier, S. Gimeno García, P. Hochstaffl, and S. Städt, “Py4CAtS—PYthon for Computational ATmospheric Spectroscopy,” Atmosphere, vol. 10, no. 5, Art. no. 5, May 2019, doi: 10.3390/atmos10050262. [Link to article](https://www.mdpi.com/2073-4433/10/5/262)
