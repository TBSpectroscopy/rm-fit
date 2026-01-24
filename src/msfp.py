#-------------------------------------------------------------------------------------------------
#
# Thibault BERTIN
# Atomic and Molecular Physics
# Center for Astrophysics | Harvard & Smithsonian
# 60 Garden St., 02138 MA, USA
# E-mail: thibault.bertin@cfa.harvard.edu
#
#-------------------------------------------------------------------------------------------------


import os
import ctrl_pars
import linelist
import calc
import tips
import opus
import fit
import output
import numpy as np
import scipy.optimize as optimize
import sys

def msfp(control_file, option):

    # Program start-up
    if not os.path.exists(control_file) :
        print("ERROR: file \"{}\" does not exist".format(control_file))
        return

    # Get parameters from control file
    spectral_data, spectral_fit = ctrl_pars.get(control_file)
    if spectral_data == dict() :
        return

    fold_path = ''
    folders = spectral_data["calculation"]["output_path"].split(os.sep)
    for i in range(0, len(folders) - 1, 1):
        fold_path = "{}{}{}".format(fold_path, folders[i], os.sep)
        if not os.path.exists(fold_path):
            print("WARNING: creating the non-existent output directory provided (\"{}\")" \
            .format(fold_path))
            os.mkdir(fold_path)
    

    linelists, offdiags = linelist.get_blocks(spectral_data) # Separate line parameters according to their interactions with each other -> [file index][block (lines mixed) index]
    line_fit = linelist.get_fitted_parameters(linelists, offdiags)

    
    specs = [np.array([], dtype = np.double), np.array([], dtype = np.double)]  # [all wavenumbers, all irradiances]
    for i in range(0, len(spectral_data["spectra"]), 1):
        spec = opus.readSpectrum(spectral_data["spectra"][i]["spectrum"], spectral_data["calculation"]["range"])
        #spec = calc.calibrate_spectrum(spec[0], spec[1], spectral_data["calculation"]["x_calibration_factor"], spectral_data["spectra"][i]["baseline"])
        # Bad idea to calibrate the spectrum before fitting. The more abs(xcalibration) increases in the fit, the more experimental data is lost.
        # Calibration should be applied to the calculations during the fit, and to the experiment after.
        spectral_data["spectra"][i]["delimitations"] = [len(specs[0]), len(specs[0]) + len(spec[0])]
        if len(spec[0]) < 2:
            print("ERROR: spectral region in control file is out of boundaries of experiment.\nRM-Fit is a fitting software and does not currently support calculations without a spectrum file. For these use cases we recommend HAPI or HAPI2.")
            sys.exit()

        specs[0] = np.concatenate((specs[0], spec[0]), axis = None)
        specs[1] = np.concatenate((specs[1], spec[1]), axis = None)


    params_id, params = fit.group_fitted_parameters(spectral_fit, line_fit) # Concatenate all parameters into one array

    # Check if file is a directory
    if os.path.isdir(spectral_data["calculation"]["output_path"]):
        print("\nOutput set to folder.\nUsing {}rm-fit_out as output files".format(spectral_data["calculation"]["output_path"] + ("/" * (spectral_data["calculation"]["output_path"][-1] != os.sep))))
        spectral_data["calculation"]["output_path"] = spectral_data["calculation"]["output_path"] + ("/" * (spectral_data["calculation"]["output_path"][-1] != os.sep)) +"rm-fit_out"

    # Check before overwriting files
    spectral_data["calculation"]["output_path"] = os.path.splitext(spectral_data["calculation"]["output_path"])[0] # Remove extension
    spectral_data["calculation"]["log_file"] = os.path.splitext(spectral_data["calculation"]["output_path"])[0] # Remove extension
    if os.path.exists("{}.txt".format(spectral_data["calculation"]["output_path"])) or os.path.exists("{}_spectra.txt".format(spectral_data["calculation"]["output_path"])):
        print("\nOutput file exists already.")
        overwrite = input("Do you want to overwrite it? (y/n)\n") == "y"
        if not overwrite:
            count = 1
            spectral_data["calculation"]["output_path"] = "{}_{:d}".format(spectral_data["calculation"]["output_path"], count) # Add _1 to the file name
            while os.path.exists("{}.txt".format(spectral_data["calculation"]["output_path"])) or os.path.exists("{}_spectra.txt".format(spectral_data["calculation"]["output_path"])): # Increase number of _n if file exists
                count += 1
                dlen = len(spectral_data["calculation"]["output_path"].split("_")[-1])
                spectral_data["calculation"]["output_path"] = "{}{:d}".format(spectral_data["calculation"]["output_path"][:-dlen], count)

    # Either fit ("-f") or calculate ("-c") spectra depending on command parameter
    if option == "-f":

        y_resid = fit.fit_spectra(params, specs, params_id, len(spectral_fit), spectral_data, linelists, offdiags)

        y_calc = np.array([])
        spec_calib = [np.array([], dtype = np.double), np.array([], dtype = np.double)]
        for i in range(0, len(spectral_data["spectra"]), 1):
            lims = spectral_data["spectra"][i]["delimitations"]
            calib = calc.calibrate_spectrum(specs[0][lims[0] : lims[1]], specs[1][lims[0] : lims[1]], spectral_data["calculation"]["x_calibration_factor"])
            spec_calib[0] = np.concatenate((spec_calib[0], calib[0]), axis = None)
            spec_calib[1] = np.concatenate((spec_calib[1], calib[1]), axis = None)
        del specs
        for i in range(0, len(spectral_data["spectra"]), 1):
            y_calc = np.concatenate((y_calc, calc.calc_spectrum(spectral_data, i, linelists, offdiags, spec_calib[0][spectral_data["spectra"][i]["delimitations"][0] : spectral_data["spectra"][i]["delimitations"][1]], apply_xcal = False)), axis = None)

        y_resid = spec_calib[1] - y_calc

        
        spec_filename = "{}_spectra.txt".format(os.path.splitext(spectral_data["calculation"]["output_path"])[0])
        with open(spec_filename, "w") as f:
            for i in range(0, len(spec_calib[0]), 1):
                f.write("{:18.10f}{:17.7E}{:17.7E}{:17.7E}\n".format(spec_calib[0][i], spec_calib[1][i], y_calc[i], y_resid[i]))

        return

    elif option == "-c":
        y_calc = np.array([])
        spec_calib = [np.array([], dtype = np.double), np.array([], dtype = np.double)]
        for i in range(0, len(spectral_data["spectra"]), 1):
            lims = spectral_data["spectra"][i]["delimitations"]
            calib = calc.calibrate_spectrum(specs[0][lims[0] : lims[1]], specs[1][lims[0] : lims[1]], spectral_data["calculation"]["x_calibration_factor"])
            spec_calib[0] = np.concatenate((spec_calib[0], calib[0]), axis = None)
            spec_calib[1] = np.concatenate((spec_calib[1], calib[1]), axis = None)
        del specs
        for i in range(0, len(spectral_data["spectra"]), 1):
            y_calc = np.concatenate((y_calc, fit.recalc_spectrum(params, i, spec_calib, params_id, len(spectral_fit), spectral_data, linelists, offdiags, apply_xcal = False)), axis = None)

        spec_filename = "{}_spectra.txt".format(os.path.splitext(spectral_data["calculation"]["output_path"])[0])
        with open(spec_filename, "w") as f:
            for i in range(0, len(spec_calib[0]), 1):
                f.write("{:18.10f}{:17.7E}{:17.7E}{:17.7E}\n".format(spec_calib[0][i], spec_calib[1][i], y_calc[i], spec_calib[1][i] - y_calc[i]))


    return




