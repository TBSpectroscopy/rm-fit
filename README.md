# RM-Fit

What is RM-Fit?
---------------

RM-Fit (Relaxation Matrix-Fit) is a free (libre) multi-spectrum fitting python script using the relaxation matrix formalism (C++ and Fortran). It is capable of modeling full line mixing in combination with speed-dependent and narrowing effects, as well as first-order line-mixing.

Citation
--------

If you use RM-Fit in your work, please cite it using the following references:

T. Bertin and J. Vander Auwera, “CO2 collision-induced line parameters for the ν3 band of 12CH4 measured using a hard-collision speed-dependent line shape and the relaxation matrix formalism,” Journal of Quantitative Spectroscopy and Radiative Transfer, vol. 324, p. 109069, Sep. 2024, doi: 10.1016/j.jqsrt.2024.109069. [Link to article](https://www.sciencedirect.com/science/article/abs/pii/S0022407324001766)

H. Tran, N. H. Ngo, and J. Hartmann, “Efficient computation of some speed-dependent isolated line profiles,” Journal of Quantitative Spectroscopy and Radiative Transfer, vol. 129, pp. 199–203, Nov. 2013, doi: 10.1016/j.jqsrt.2013.06.015. [Link to article](https://www.sciencedirect.com/science/article/abs/pii/S0022407313002598)

F. Schreier, S. Gimeno García, P. Hochstaffl, and S. Städt, “Py4CAtS—PYthon for Computational ATmospheric Spectroscopy,” Atmosphere, vol. 10, no. 5, Art. no. 5, May 2019, doi: 10.3390/atmos10050262. [Link to article](https://www.mdpi.com/2073-4433/10/5/262)

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

2. Run `install.sh`

3. Start using RM-Fit

OR

2. Access the source code directory
    ```
    cd rm-fit/src
    ```

3. Compile the Fortran and C++ libraries
    ```
    python3 -m numpy.f2py -c qSDHC.pyf qSDHC.for CPF3.for CPF.for --f90flags="-Ofast"
    python3 -m numpy.f2py -c mHT.pyf mHT.for CPF3.for CPF.for --f90flags="-Ofast"
    python3 -m numpy.f2py -c qSDV.pyf qSDV.for CPF3.for CPF.for --f90flags="-Ofast"
    python3 -m numpy.f2py -c pCqSDHC.pyf pCqSDHC.for CPF3.for CPF.for --f90flags="-Ofast"
    mkdir build
    cd build
    cmake ..
    make
    ```

4. Start using RM-Fit

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

Note: rm-fit.py can be called from any directory.

### Dependencies

- Python 3.x
- NumPy
- SciPy

### Input

An example is provided in `data/input/` named `ch4_co2.txt`. The input file contains all the relevant information to calculate spectra. It is separated into blocks marked by a '$' symbol followed by a name. What follows the name on the same line is ignored. All file/folder paths inside must be absolute or relative to the location of the input file.

`$SPECTRUM` contains all the information specific to one spectrum. Can be as many as needed. Note that this repository does not provide the spectrum to run the example input. [However, we do provide it as a download here](https://owncloud.ulb.ac.be/index.php/s/yLinSsQ1t480Qxc/download) or [here.](https://www.dropbox.com/scl/fi/bdo6sxd3lfg1p3lz7tbmt/ch4_co2_spec.txt?rlkey=iw2mq2r02skcasxabfyp2h6na&st=ox0a3jlx&dl=1)

`$LINELIST` contains the inputs for a species/isotopologue line list (diagonals and off-diagonals). Can be as many as needed

`$CALCULATION` contains the information used for all spectra and the output path. Unique

Initial parameters are followed by an 'f' if you want to fit them, or 'c' if you want them to stay constant. The parameters inside the blocks are detailed in the example file.

### Linelist

Lines needed to compute the example input file are provided by the linelists `diag_61.txt` and `diag_62.txt` and the off-diagonal list `offdiag_61.txt`. All files must contain a header giving details about the species, and a format section inside quotation marks allowing RM-Fit to know how to read the parameters. The format follows a structure similar to python's string `format` method, with braces containing the name of a parameter and its width separated by a colon. The width can include a letter indicating how the numbers should be formatted. As an example, `{foreign-broadening:10E}` indicates that $\gamma _0 ^{\mathrm{foreign}}$ is written with 10 characters in scientific E notation.

Each parameter that needs to be fitted must include a column `parameter_unc` and `parameter_fit` in the format, where `parameter` is replaced by the key of the parameter (the list of available keys is given below). `parameter_unc` is used to output the uncertainty of the fit, and `parameter_fit` is used to indicate which parameter to fit:
```
{foreign-broadening:10E} {foreign-broadening_unc:7E}{foreign-broadening_fit:1}
```
Here, 10 characters are used for the foreign-broadening parameter, 1 for an empty space, 7 for the uncertainty of the fitted broadening, and 1 to tell RM-Fit if the parameter is to be fitted. To activate the fit of a parameter in a certain line, replace the space at the position of `parameter_fit` and `fit` by an asterisk, `*`.

Spaces are considered inside the format section and can be used to ignore parts of the lists (ignored parameters are assumed to be 0). You can also use a format keyword that is not used by RM-Fit. For example: `{empty:20}` means that the next 20 characters are ignored.

#### Available format keys

##### Diagonal

| Format key | Description |
|---|---|
| fit | Activate the fit for the line |
| name | Identifier used by the off-diagonal |
| wavenumber | $\tilde{\nu }_{ij}$, transition wavenumber |
| intensity | $S_{ij}$, line intensity |
| foreign-broadening | $\gamma _0 ^{\mathrm{foreign}}$, foreign speed-averaged broadening parameter |
| foreign-broadening_temp | Temperature exponent 1 for $\gamma _0 ^{\mathrm{foreign}}$ |
| self-broadening | $\gamma _0 ^{\mathrm{self}}$, self speed-averaged broadening parameter |
| self-broadening_temp | Temperature exponent 1 for $\gamma _0 ^{\mathrm{self}}$ |
| SD-broadening | $\gamma _2 / \gamma _0$, quadratic speed-dependence ratio of the broadening |
| narrowing | $\nu _{VC}$, velocity-changing collision frequency parameter |
| narrowing_temp | Temperature exponent 1 for $\nu _{VC}$ |
| narrowing_i | $\mathrm{Im}( \nu _{\mathrm{opt}} )$, imaginary part of the Dicke narrowing parameter |
| narrowing_i_temp | Temperature exponent 1 for $\mathrm{Im}( \nu _{\mathrm{opt}} )$ |
| foreign-shift | $\delta _0 ^{\mathrm{foreign}}$, foreign speed-averaged shift parameter |
| foreign-shift_temp | Temperature exponent 1 for $\delta _0 ^{\mathrm{foreign}}$ |
| self-shift | $\delta _0 ^{\mathrm{self}}$, self speed-averaged shift parameter |
| self-shift_temp | Temperature exponent 1 for $\delta _0 ^{\mathrm{self}}$ |
| SD-shift | $\delta _2 / \delta _0$, quadratic speed-dependence ratio of the shift |
| correlation_HT | $\eta $, partial correlation parameter |
| line-mixing_fo | $y$, first-order line-mixing parameter |
| line-mixing_fo_temp | Temperature exponent 1 for $y$ |
| statistical weight | $g''$, statistical weight of the lower level |
| energy | Lower energy level |

##### Off-diagonal

| Format key | Description |
|---|---|
| name | Line identifier from the diagonal (used to indicate lines coupled) |
| line-mixing | $W_{mn}$, line-mixing parameter |
| line-mixing_temp | Temperature exponent 1 for $W_{mn}$ |

#### Not yet implemented

##### Diagonal

| Format key | Description |
|---|---|
| SD-broadening_temp | Temperature exponent 1 for $\gamma _2 / \gamma _0$ |
| SD-shift_temp | Temperature exponent 1 for $\delta _2 / \delta _0$ |
| foreign-SD-broadening | $(\gamma _2 / \gamma _0)^{\mathrm{foreign}}$, foreign quadratic speed-dependence ratio of the broadening |
| foreign-SD-broadening_temp | Temperature exponent 1 for $(\gamma _2 / \gamma _0)^{\mathrm{foreign}}$ |
| self-SD-broadening | $(\gamma _2 / \gamma _0)^{\mathrm{self}}$, self quadratic speed-dependence ratio of the broadening |
| self-SD-broadening_temp | Temperature exponent 1 for $(\gamma _2 / \gamma _0)^{\mathrm{self}}$ |
| foreign-SD-shift | $(\delta _2 / \delta _0)^{\mathrm{foreign}}$, foreign quadratic speed-dependence ratio of the shift |
| foreign-SD-shift_temp | Temperature exponent 1 for $(\delta _2 / \delta _0)^{\mathrm{foreign}}$ |
| self-SD-shift | $(\delta _2 / \delta _0)^{\mathrm{self}}$, self quadratic speed-dependence ratio of the shift |
| self-SD-shift_temp | Temperature exponent 1 for $(\delta _2 / \delta _0)^{\mathrm{self}}$ |
| foreign-broadening_2 | $\gamma _2 ^{\mathrm{foreign}}$, foreign order 2 speed-dependent broadening parameter |
| foreign-broadening_2_temp | Temperature exponent 1 for $\gamma _2 ^{\mathrm{foreign}}$ |
| self-broadening_2 | $\gamma _2 ^{\mathrm{self}}$, self order 2 speed-dependent broadening parameter |
| self-broadening_2_temp | Temperature exponent 1 for $\gamma _2 ^{\mathrm{self}}$ |
| foreign-shift_2 | $\delta _2 ^{\mathrm{foreign}}$, foreign order 2 speed-dependent shift parameter |
| foreign-shift_2_temp | Temperature exponent 1 for $\delta _2 ^{\mathrm{foreign}}$ |
| self-shift_2 | $\delta _2 ^{\mathrm{self}}$, self order 2 speed-dependent shift parameter |
| self-shift_2_temp | Temperature exponent 1 for $\delta _2 ^{\mathrm{self}}$ |
| foreign-narrowing | $\nu _{VC} ^{\mathrm{foreign}}$, foreign velocity-changing collision frequency parameter |
| foreign-narrowing_temp | Temperature exponent 1 for $\nu _{VC} ^{\mathrm{foreign}}$ |
| self-narrowing | $\nu _{VC} ^{\mathrm{self}}$, self velocity-changing collision frequency parameter |
| self-narrowing_temp | Temperature exponent 1 for $\nu _{VC} ^{\mathrm{self}}$ |
| foreign-narrowing_i | $\mathrm{Im}( \nu _{\mathrm{opt}} ^{\mathrm{foreign}} )$, imaginary part of the foreign Dicke narrowing parameter |
| foreign-narrowing_i_temp | Temperature exponent 1 for $\mathrm{Im}( \nu _{\mathrm{opt}} ^{\mathrm{foreign}} )$ |
| self-narrowing_i | $\mathrm{Im}( \nu _{\mathrm{opt}} ^{\mathrm{self}} )$, imaginary part of the self Dicke narrowing parameter |
| self-narrowing_i_temp | Temperature exponent 1 for $\mathrm{Im}( \nu _{\mathrm{opt}} ^{\mathrm{self}} )$ |
| foreign-correlation_HT | $\eta ^{\mathrm{foreign}}$, foreign partial correlation parameter |
| self-correlation_HT | $\eta ^{\mathrm{self}}$, self partial correlation parameter |

##### Off-diagonal

| Format key | Description |
|---|---|
| foreign_line-mixing | $W_{mn} ^{\mathrm{foreign}}$, foreign line-mixing parameter |
| foreign_line-mixing_temp | Temperature exponent 1 for $W_{mn} ^{\mathrm{foreign}}$ |
| self_line-mixing | $W_{mn} ^{\mathrm{self}}$, self line-mixing parameter |
| self_line-mixing_temp | Temperature exponent 1 for $W_{mn} ^{\mathrm{self}}$ |
