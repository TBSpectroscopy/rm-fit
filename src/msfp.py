#-------------------------------------------------------------------------------------------------
#
# Thibault BERTIN
# Spectroscopy, Quantum Chemistry and Atmospheric Remote Sensing (SQUARES), C.P. 160/09
# Universite Libre de Bruxelles
# 50 avenue F. D. Roosevelt, B-1050 Brussels, Belgium
# Phone: +32.2.650.24.18 - E-mail: thibault.bertin@ulb.be - Web: http://www.ulb.ac.be/cpm
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
import time

def msfp(control_file, option):

    # Program start-up
    if not os.path.exists(control_file) :
        print("ERROR - File \"{}\" does not exist.".format(control_file))
        return

    # Get parameters from control file
    spectral_data, spectral_fit = ctrl_pars.get(control_file)
    if spectral_data == dict() :
        return

    if spectral_data["calculation"]["output_path"][-1] == os.sep :
        spectral_data["calculation"]["output_path"] = spectral_data["calculation"]["output_path"][:-1]

    fold_path = ''
    folders = spectral_data["calculation"]["output_path"].split(os.sep)
    for i in range(0, len(folders) - 1, 1):
        fold_path = "{}{}/".format(fold_path, folders[i])
        if not os.path.exists(fold_path):
            print("WARNING - Creating the non-existent output directoy provided (\"{}\")." \
            .format(fold_path))
            os.mkdir(fold_path)
    

    linelists, offdiags = linelist.get_blocks(spectral_data) # Separate line parameters according to their interactions with each other -> [file index][block (lines mixed) index]
    line_fit = linelist.get_fitted_parameters(linelists, offdiags)

    
    specs = [np.array([], dtype = np.double), np.array([], dtype = np.double)]  # [all wavenumbers, all irradiances]
    nu = np.array([], dtype = np.double)
    for i in range(0, len(spectral_data["spectra"]), 1):
        spec = opus.readSpectrum(spectral_data["spectra"][i]["spectrum"], spectral_data["calculation"]["range"])
        nu = spec[0]
        specs[0] = np.concatenate((specs[0], spec[0]), axis = None)
        specs[1] = np.concatenate((specs[1], spec[1]), axis = None)


    params_id, params = fit.group_fitted_parameters(spectral_fit, line_fit) # Concatenate all parameters into one array

    # Either fit ("-f") or calculate ("-c") spectra depending on command parameter
    if option == "-f":

        if spectral_data["calculation"]["output_path"][-1] == os.sep :  # Remove folder separator
            spectral_data["calculation"]["output_path"] = spectral_data["calculation"]["output_path"][:-1]

        spectral_data["calculation"]["output_path"] = os.path.splitext(spectral_data["calculation"]["output_path"])[0] # Remove extension
        spectral_data["calculation"]["log_file"] = os.path.splitext(spectral_data["calculation"]["output_path"])[0] # Remove extension
        if os.path.exists("{}.txt".format(spectral_data["calculation"]["output_path"])):
            print("\nOutput file exists already.")
            overwrite = input("Do you want to overwrite it? (y/n)\n") == "y"
            if not overwrite:
                count = 1
                spectral_data["calculation"]["output_path"] = "{}_{:d}".format(spectral_data["calculation"]["output_path"], count)
                while os.path.exists("{}.txt".format(spectral_data["calculation"]["output_path"])):
                    count += 1
                    dlen = len(spectral_data["calculation"]["output_path"].split("_")[-1])
                    spectral_data["calculation"]["output_path"] = "{}{:d}".format(spectral_data["calculation"]["output_path"][:-dlen], count)

        y_resid = fit.fit_spectra(params, specs, params_id, len(spectral_fit), spectral_data, linelists, offdiags)
        y_calc = y_resid + specs[1]

        
        spec_filename = "{}_spectra.txt".format(os.path.splitext(spectral_data["calculation"]["output_path"])[0])
        with open(spec_filename, "w") as f:
            for i in range(0, len(specs[0]), 1):
                f.write("{:18.10f}{:17.7E}{:17.7E}{:17.7E}\n".format(specs[0][i], specs[1][i], y_calc[i], -y_resid[i]))

        return

    elif option == "-c":
        starttime = time.time()
        y_calc = np.array([])
        for i in range(0, len(spectral_data["spectra"]), 1):
            y_calc = np.concatenate((y_calc, fit.recalc_spectrum(params, i, specs, params_id, len(spectral_fit), spectral_data, linelists, offdiags)), axis = None)

        spec_filename = "{}_spectra.txt".format(os.path.splitext(spectral_data["calculation"]["output_path"])[0])
        with open(spec_filename, "w") as f:
            for i in range(0, len(specs[0]), 1):
                f.write("{:18.10f}{:17.7E}{:17.7E}{:17.7E}\n".format(specs[0][i], specs[1][i], y_calc[i], specs[1][i] - y_calc[i]))

    else:
        print("ERROR: unrecognized option \"{}\"\nTry \"python rm-fit -h\" or \"python rm-fit --help\" for more information".format(option))

    return




