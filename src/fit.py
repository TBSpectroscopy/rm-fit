#-------------------------------------------------------------------------------------------------
#
# Thibault BERTIN
# Atomic and Molecular Physics
# Center for Astrophysics | Harvard & Smithsonian
# 60 Garden St., 02138 MA, USA
# E-mail: thibault.bertin@cfa.harvard.edu
#
#-------------------------------------------------------------------------------------------------


import numpy as np
from scipy.optimize import least_squares
import math
import calc
import datetime
import sys
import os
import output
import multiprocessing


#-------------------------------------------------------------------------------------------------
def group_fitted_parameters(spectral_fit, line_fit):
    # spectral_fit:    file parameter value
    # line_fit :       file block line parameter value
    fitted_par = [() for _ in range(0, (len(spectral_fit) + len(line_fit)), 1)] # identifier
    fitted_val = [0.0] * (len(spectral_fit) + len(line_fit)) # value

    for i in range(0, len(spectral_fit), 1):
        fitted_par[i] = spectral_fit[i][:-1]
        fitted_val[i] = spectral_fit[i][-1]

    for i in range(len(spectral_fit), len(fitted_par), 1):
        fitted_par[i] = line_fit[i - len(spectral_fit)][:-1]
        if fitted_par[i][3] == "intensity":
            fitted_val[i] = line_fit[i - len(spectral_fit)][-1] * 1E20
        else:
            fitted_val[i] = line_fit[i - len(spectral_fit)][-1]

    return fitted_par, fitted_val


#-------------------------------------------------------------------------------------------------
def get_rel_y(params_id : list, params : list, linelists, offdiags, method, option = "-f"):

    y = [[[] for j in range(0, len(linelists[i]), 1)] for i in range(0, len(linelists), 1)]
    for i in range(0, len(params), 1):
        if len(params_id[i]) == 4:
            if params_id[i][3] == "line-mixing_fo":
                y[params_id[i][0]][params_id[i][1]].append(i)

    rel_y = [[() for j in range(0, len(linelists[i]), 1)] for i in range(0, len(linelists), 1)]
    del_y = []
    for i in range(0, len(y), 1):
        for j in range(0, len(y[i]), 1):
            for k in range(0, len(y[i][j]), 1):
                if method == "first-order" or offdiags[params_id[y[i][j][k]][0]][params_id[y[i][j][k]][1]] == {} or offdiags[params_id[y[i][j][k]][0]][params_id[y[i][j][k]][1]] == []:
                    if k == int(len(y[i][j])/2):
                        rel_y[i][j] = params_id[y[i][j][k]]
                        del_y.append(y[i][j][k])
                        if option == "-f" or option == "--fit":
                            update_y(params_id[y[i][j][k]], rel_y, linelists) # Update the selected parameter
                        break
                else: # Keep deleting if the full line-mixing is used
                    del_y.append(y[i][j][k])

    print("Parameters set to enforce sum rule")
    print(rel_y)
    for i in range(len(del_y) - 1, -1, -1):
        del params[del_y[i]]
        del params_id[del_y[i]]


    return rel_y


#-------------------------------------------------------------------------------------------------
def get_hessian_inv(jacobian, params_id, params, par_prob):
    hessian = np.matmul(np.transpose(jacobian), jacobian)   # A Gauss Newton approximation of the Hessian of the cost function

    try:
        hessian_inv = np.linalg.inv(hessian)
    except:
        print("\nWARNING: Fit could not converge: unable to get uncertainties")
        jac_norm = np.linalg.norm(jacobian, axis = 0)
        recurs = False
        for i in range(0, len(jac_norm), 1):
            if jac_norm[i] < 1e-16: # Detection limit for the fit
                recurs = True
                print("{} {:e} could not converge".format(params_id[i], params[i]))
                par_prob.append(i)
                jacobian[:,i] += 1e-10
        if recurs:
            print("\nAttempting uncertainty calculation without problematic parameters")
            hessian_inv = get_hessian_inv(jacobian, params_id, params, par_prob)
            return hessian_inv
        print("\nUnable to find the problem. Check your initial parameters")
        hessian_inv = np.identity(len(hessian))
    return hessian_inv


#-------------------------------------------------------------------------------------------------
def fit_spectra(params, specs, params_id, ns, spectral_data, linelists, offdiags, rel_y : list):
    # ns: number of spectrum specific parameters

    results_filename = "{}.txt".format(os.path.splitext(spectral_data["calculation"]["output_path"])[0])
    log_filename = "{}_log.txt".format(os.path.splitext(spectral_data["calculation"]["log_file"])[0])
    params = np.array(params, dtype = np.double)

    bounds = get_bounds(params_id, params)

    for i in range(0, len(params), 1):
        print("{}  {}  {} - {}".format(params_id[i], params[i], bounds[0][i], bounds[1][i]))

    noend = False
    if os.path.exists(log_filename):
        noend = True
        with open(log_filename) as f:
            for line in f:
                if line.startswith("END"):
                    noend = False
    if noend:
        print("\nPrevious fit did not end properly.")
        cont = input("Continue from last results? (y/n)\n") == "y"
        if cont:
            with open(log_filename) as f:
                for x, line in enumerate(f):
                    if x >= 2:
                        spl = line.split()
                        if line.startswith("["):
                            params = [float(spl[i]) for i in range(1, len(spl), 1)]
                        elif line.endswith("]\n"):
                            for i in range(0, len(spl) - 1, 1):
                                params.append(float(spl[i]))
                            params.append(float(spl[-1][:-1]))
                        else:
                            for i in range(0, len(spl), 1):
                                params.append(float(spl[i]))
        else:
            with open(log_filename, "w") as f:
                f.write("Run time: {}\n\n".format(datetime.datetime.now().isoformat(sep = " ", timespec="minutes")))
    else:
        with open(log_filename, "w") as f:
            f.write("Run time: {}\n\n".format(datetime.datetime.now().isoformat(sep = " ", timespec="minutes")))



    with open(results_filename, "w") as f:

        f.write("RM-Fit\n"\
                "Run time: {}\n\n".format(datetime.datetime.now().isoformat(sep = " ", timespec="minutes")))

        results = least_squares(leastsq_fun, params, jac = calc_jac, bounds = bounds, args = (specs, params_id, ns, spectral_data, linelists, offdiags, rel_y), method = "trf", verbose = 2, max_nfev = 12, gtol = None, ftol = None, xtol = 1e-12)

        # Estimate the uncertainty with the rescaled inverse Hessian
        stdev = math.sqrt(np.sum((results.fun)**2)/float((results.fun).size - len(params)))
        jacobian = results.jac

        par_prob = []
        hessian_inv = get_hessian_inv(jacobian, params_id, results.x, par_prob)

        main_diag = np.diagonal(hessian_inv)
        sqrt_main_diag = np.sqrt(main_diag)
        unc = stdev*sqrt_main_diag

        update_parameters(results.x, specs, params_id, ns, spectral_data, linelists, offdiags, rel_y, unc = unc)
        for i in range(0, len(results.x), 1):
            f.write("{:13.5E} +- {:7.1E}..........{}\n".format(results.x[i], unc[i], params_id[i]))


        f.write("\nEnd time: {}".format(datetime.datetime.now().isoformat(sep = " ", timespec="minutes")))
        linelist_inputs = [i["linelist_in"] for i in spectral_data["linelists"]]
        linelist_outputs = [i["linelist_out"] for i in spectral_data["linelists"]]
        offdiag_inputs = [i["off_diagonal_in"] for i in spectral_data["linelists"]]
        offdiag_outputs = [i["off_diagonal_out"] for i in spectral_data["linelists"]]
        output.write_linelists(linelist_inputs, linelist_outputs, linelists, spectral_data["calculation"]["range"], [i["format"] for i in spectral_data["linelists"]])
        output.write_offdiags(offdiag_inputs, offdiag_outputs, offdiags, [i["format_offdiag"] for i in spectral_data["linelists"]])

        hessian = np.matmul(np.transpose(jacobian), jacobian)   # A Gauss Newton approximation of the Hessian of the cost function
        with open("{}_matrices.txt".format(os.path.splitext(results_filename)[0]), "w") as f2:
            f2.write("Hessian & Covariance matrices\n\n")
            f2.write("{:7}".format(""))
            for j in range(0, len(results.x), 1):
                f2.write(" {:4d}      ".format(j + 1))
                if j == len(results.x) - 1:
                    f2.write("\n")
            for i in range(0, len(results.x), 1):
                for j in range(0, len(results.x), 1):
                    if j == 0:
                        f2.write("{:4d}{:3}".format(i + 1, ""))
                    f2.write("{:8.1E}   ".format(hessian[i][j]))
                    if j == len(results.x) - 1:
                        f2.write("\n")
            f2.write("\n\n")

            f2.write("{:7}".format(""))
            for j in range(0, len(results.x), 1):
                f2.write(" {:4d}      ".format(j + 1))
                if j == len(results.x) - 1:
                    f2.write("\n")
            for i in range(0, len(results.x), 1):
                for j in range(0, len(results.x), 1):
                    if j == 0:
                        f2.write("{:4d}{:3}".format(i + 1, ""))
                    f2.write("{:8.1E}   ".format(stdev * hessian_inv[i][j]))
                    if j == len(results.x) - 1:
                        f2.write("\n")

    with open(log_filename, "a") as f:
        f.write("\n\nEND")


    return results.fun


#-------------------------------------------------------------------------------------------------
def update_parameters(params, *args, unc = []):
    nu = args[0][0]
    params_id = args[1]
    ns = args[2]
    spectral_data = args[3]
    linelists = args[4]
    offdiags = args[5]
    rel_y = args[6]
    if len(unc) == 0:
        unc = [0.0] * len(params)


    # Update parameters
    for i in range(0, ns, 1):
        if params_id[i][1].split()[0] == "baseline":
            spectral_data["spectra"][params_id[i][0]][params_id[i][1].split()[0]][int(params_id[i][1].split()[1])] = params[i]
        elif params_id[i][1] == "x_calibration_factor":
            spectral_data["calculation"][params_id[i][1]] = params[i]
        elif params_id[i][1].split()[0] == "mole_fraction" or params_id[i][1].split()[0] == "intensity_factor":
            spectral_data["spectra"][params_id[i][0]][params_id[i][1].split()[0]][params_id[i][2]] = (params[i], spectral_data["spectra"][params_id[i][0]][params_id[i][1].split()[0]][params_id[i][2]][1])
        else:
            spectral_data["spectra"][params_id[i][0]][params_id[i][1]] = params[i]
    for i in range(ns, len(params_id), 1):
        if params_id[i][3] == "intensity":
            linelists[params_id[i][0]][params_id[i][1]][params_id[i][3]][params_id[i][2]] = (params[i] / 1E20, True, unc[i] / 1E20)
        elif params_id[i][3] != "line-mixing":
            linelists[params_id[i][0]][params_id[i][1]][params_id[i][3]][params_id[i][2]] = (params[i], True, unc[i])
            if params_id[i][3] == "line-mixing_fo":
                update_y(params_id[i], rel_y, linelists)
        else:
            offdiags[params_id[i][0]][params_id[i][1]][params_id[i][3]][params_id[i][2]] = (params[i], True, unc[i])
    
    return


#-------------------------------------------------------------------------------------------------
def update_y(param_id, rel_y, linelists):
    y_sum = 0.0
    unc_sum = 0.0
    len_sum = 0
    rel_id = rel_y[param_id[0]][param_id[1]]
    for i in range(0, len(linelists[param_id[0]][param_id[1]][param_id[3]]), 1):
        if i != rel_id[2]:
            y_sum += linelists[param_id[0]][param_id[1]]["intensity"][i][0] * linelists[param_id[0]][param_id[1]]["line-mixing_fo"][i][0]
            if linelists[param_id[0]][param_id[1]]["line-mixing_fo"][i][2] != 0.0:
                unc_sum += ((linelists[param_id[0]][param_id[1]]["intensity"][i][0] / linelists[rel_id[0]][rel_id[1]]["intensity"][rel_id[2]][0])  * linelists[param_id[0]][param_id[1]]["line-mixing_fo"][i][2]) ** 2 # unc_a = sqrt(Sum((S_i/S_a) unc_i)^2) / sqrt(n)
                len_sum += 1

    linelists[rel_id[0]][rel_id[1]]["line-mixing_fo"][rel_id[2]] = (-(y_sum / linelists[rel_id[0]][rel_id[1]]["intensity"][rel_id[2]][0]), True, np.sqrt(unc_sum / max(1, len_sum)))
    return


#-------------------------------------------------------------------------------------------------
def get_bounds(params_id, params):

    bound_min = []
    bound_max = []
    for i in range(0, len(params_id), 1):
        if len(params_id[i]) != 4:
            if params_id[i][1].split()[0] == "mole_fraction":
                bound_min.append(0.0)
                bound_max.append(1.0)
            elif params_id[i][1].split()[0] == "x_calibration_factor" or params_id[i][1].split()[0] == "baseline":    # Also acccept baseline 0 to be negative just in case
                bound_min.append(-np.inf)
                bound_max.append(np.inf)
            else:
                bound_min.append(0.0)
                bound_max.append(np.inf)
        else:
            if params_id[i][3] in ["foreign-shift", "self-shift", "narrowing_i"]:
                bound_min.append(-np.inf)
                bound_max.append(np.inf)
            elif params_id[i][3] in ["line-mixing", "line-mixing_fo"]:
                lims = [0.0, params[i] * np.inf]
                bound_min.append(min(lims))
                bound_max.append(max(lims))
            else:
                bound_min.append(0.0)
                bound_max.append(np.inf)

    return (bound_min, bound_max)

#-------------------------------------------------------------------------------------------------
def recalc_spectrum(params, spec_index, *args, apply_xcal : bool = True):

    # Update parameters
    update_parameters(params, *args)
    # Note: lists in Python are passed by reference so args itself is modified

    nu = args[0][0]
    params_id = args[1]
    ns = args[2]
    spectral_data = args[3]
    linelists = args[4]
    offdiags = args[5]

    # Get n_points for single spectrum
    lims = spectral_data["spectra"][spec_index]["delimitations"]
    
    return calc.calc_spectrum(spectral_data, spec_index, linelists, offdiags, nu[lims[0] : lims[1]], apply_xcal = apply_xcal)


#-------------------------------------------------------------------------------------------------
# Get jacobian column for line parameters
def get_jac_column_specs(params, i, col, y0, lims, transmittance, *args):

    # Get arguments
    nu = args[0][0]
    params_id = args[1]
    spectral_data = args[3]

    # Estimate dy/dx
    dx = np.abs(0.001 * params[i])
    params[i] += dx
    print("{} {}".format(params_id[i], params[i]))
    if params_id[i][1].split()[0] == "baseline":
        x = nu[lims[params_id[i][0]][0] : lims[params_id[i][0]][1]] - ((nu[lims[params_id[i][0]][1] - 1] + nu[lims[params_id[i][0]][0]])/2)
        if spectral_data["spectra"][params_id[i][0]]["form"] == "transmittance":
            col[lims[params_id[i][0]][0] : lims[params_id[i][0]][1]] = transmittance[params_id[i][0]] * (x ** int(params_id[i][1].split()[1]))
        elif spectral_data["spectra"][params_id[i][0]]["form"] == "absorption_coefficient":
            col[lims[params_id[i][0]][0] : lims[params_id[i][0]][1]] = x ** int(params_id[i][1].split()[1])
    elif params_id[i][0] == None:
        if dx != 0.0:
            for j in range(0, len(spectral_data["spectra"]), 1):
                dy = recalc_spectrum(params, j, *args) - y0[lims[j][0] : lims[j][1]]
                col[lims[j][0] : lims[j][1]] = dy/dx
        else:
            col = np.zeros(len(nu), dtype = np.double)
    else:
        if dx != 0.0:
            dy = recalc_spectrum(params, params_id[i][0], *args) - y0[lims[params_id[i][0]][0] : lims[params_id[i][0]][1]]
            col[lims[params_id[i][0]][0] : lims[params_id[i][0]][1]] = dy/dx
        else:
            col[lims[params_id[i][0]][0] : lims[params_id[i][0]][1]] = np.zeros(lims[params_id[i][0]][1] - lims[params_id[i][0]][0], dtype = np.double)


    params[i] -= dx

    return col


#-------------------------------------------------------------------------------------------------
# Get jacobian column for line parameters
def get_jac_column_line(params, i, col, y0, lims, *args):
    # Get arguments
    nu = args[0][0]
    params_id = args[1]
    spectral_data = args[3]
    linelists = args[4]
    rel_y = args[6]

    dx = 0.001 * params[i]
    if params_id[i][-1] == "wavenumber":
        dx = 0.0001
    params[i] += dx
    print("{} {}".format(params_id[i], params[i]))
    if dx != 0.0:
        for j in range(0, len(spectral_data["spectra"]), 1):
            dy = recalc_spectrum(params, j, *args) - y0[lims[j][0] : lims[j][1]]
            col[lims[j][0] : lims[j][1]] = dy/dx
    else:
        col = np.zeros(len(nu), dtype = np.double)

    params[i] -= dx
    if params_id[i][3] == "line-mixing_fo":
        update_y(params_id[i], rel_y, linelists)

    return col


#-------------------------------------------------------------------------------------------------
# Calculate jacobian for (ycalc - yobs) (equal to jacobian for ycalc)
def calc_jac(params, *args):

    print("Calculating jacobian")
    
    # Get arguments
    nu = args[0][0]
    params_id = args[1]
    ns = args[2]
    spectral_data = args[3]


    # Calculate initial spectrum
    y0 = np.array([])
    for i in range(0, len(spectral_data["spectra"]), 1):
        y0 = np.concatenate((y0, recalc_spectrum(params, i, *args)), axis = None)

    jac = np.zeros((len(nu), len(params_id)), dtype = np.double)

    lims = [[spectral_data["spectra"][i]["delimitations"][0], spectral_data["spectra"][i]["delimitations"][1]] for i in range(0, len(spectral_data["spectra"]), 1)]

    # Calculate signal without baseline from calculated spectrum (used for baseline derivative)
    transmittance = []
    for i in range(0, len(spectral_data["spectra"]), 1):
        if spectral_data["spectra"][i]["form"] == "transmittance":
            transmittance.append(y0[lims[i][0] : lims[i][1]] / calc.calc_baseline(spectral_data["spectra"][params_id[i][0]]["baseline"], nu[lims[i][0] : lims[i][1]]))

    params = [i for i in params]

    # Estimate dy/dx
    pool = multiprocessing.Pool()

    inputs = [(params, i, jac[:,i], y0, lims, transmittance, *args) for i in range(0, ns, 1)]
    outputs = pool.starmap(get_jac_column_specs, inputs)
    for i in range(0, ns, 1):   # Parameters inside input file
        jac[:,i] = outputs[i]

    inputs = [(params, i, jac[:,i], y0, lims, *args) for i in range(ns, len(params), 1)]
    outputs = pool.starmap(get_jac_column_line, inputs)
    for i in range(ns, len(params), 1): # Parameters inside linelists
        jac[:,i] = outputs[i - ns]


    return jac


#-------------------------------------------------------------------------------------------------
it_count = 1
def leastsq_fun(params, *args):
    global it_count
    print("Leastsq fun\n"\
            "Iteration number {:d}".format(it_count))
    it_count += 1
    print(params)
    y_obs = args[0][1]
    spectral_data = args[3]
    y_calc = np.array([])
    for i in range(0, len(spectral_data["spectra"]), 1):
        y_calc = np.concatenate((y_calc, recalc_spectrum(params, i, *args)), axis = None)

    log_filename = "{}_log.txt".format(os.path.splitext(spectral_data["calculation"]["log_file"])[0])
    with open(log_filename, "a") as f:
        f.write("{}\n".format(params))

    
    return y_calc - y_obs


