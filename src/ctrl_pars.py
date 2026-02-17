#-------------------------------------------------------------------------------------------------
#
# Thibault BERTIN
# Atomic and Molecular Physics
# Center for Astrophysics | Harvard & Smithsonian
# 60 Garden St., 02138 MA, USA
# E-mail: thibault.bertin@cfa.harvard.edu
#
#-------------------------------------------------------------------------------------------------


import sys


def get(control_file):

    spec = []
    linelist = []
    blocpos = []
    control = []
    with open(control_file) as f:
        control = f.readlines()
        for x, line in enumerate(control):
            if line.startswith("$SPECTRUM "):
                blocpos.append({"pos" : x, "type" : "spectrum"})
                spec.append(dict())
            elif line.startswith("$LINELIST "):
                blocpos.append({"pos" : x, "type" : "linelist"})
                linelist.append(dict())
            elif line.startswith("$CALCULATION"):
                blocpos.append({"pos" : x, "type" : "calculation"})
            elif line.startswith("$FIT"):
                blocpos.append({"pos" : x, "type" : "fit"})
        blocpos.append({"pos" : x, "type" : "end"})
    calculation = dict()
    fit = dict()

    count_spectra = -1
    count_linelists = -1
    fitted_par = []
    for i in range(0, len(blocpos) - 1, 1):
        if blocpos[i]["type"] == "spectrum":
            count_spectra += 1
        elif blocpos[i]["type"] == "linelist":
            count_linelists += 1
        for j in range(blocpos[i]["pos"] + 1, blocpos[i+1]["pos"] + 1, 1):
            if blocpos[i]["type"] == "spectrum":
                if control[j].startswith("file ="):
                    spec[count_spectra]["spectrum"] = get_parameter(j, control[j])
                    spec[count_spectra]["delimitations"] = [0, 0]
                elif control[j].startswith("form ="):
                    spec[count_spectra]["form"] = get_parameter(j, control[j])
                elif control[j].startswith("baseline ="):
                    spec[count_spectra]["baseline"] = [parse_fit(i)[0] for i in (get_parameter(j, control[j]).split())]
                    truth = [parse_fit(i)[1] for i in (get_parameter(j, control[j]).split())]
                    for x, k in enumerate(truth):
                        if k:
                            fitted_par.append([count_spectra, "baseline {:d}".format(x), spec[count_spectra]["baseline"][x]])
                elif control[j].startswith("lowest_value ="):
                    spec[count_spectra]["lowest_value"] = float(get_parameter(j, control[j]))
                elif control[j].startswith("mole_fraction ="):
                    spec[count_spectra]["mole_fraction"] = [(parse_fit(k.split(":")[0])[0],) + ([int(l) for l in (k.split(":")[1]).split()],) for k in (get_parameter(j, control[j]).split("/"))]
                    truth = [parse_fit(k.split(":")[0])[1] for k in (get_parameter(j, control[j]).split("/"))]
                    for x, k in enumerate(truth):
                        if k:
                            fitted_par.append([count_spectra, "mole_fraction", x, spec[count_spectra]["mole_fraction"][x][0]])
                elif control[j].startswith("intensity_factor ="):
                    spec[count_spectra]["intensity_factor"] = [(parse_fit(k.split(":")[0])[0],) + ([int(l) for l in (k.split(":")[1]).split()],) for k in (get_parameter(j, control[j]).split("/"))]
                    truth = [parse_fit(k.split(":")[0])[1] for k in (get_parameter(j, control[j]).split("/"))]
                    for x, k in enumerate(truth):
                        if k:
                            fitted_par.append([count_spectra, "intensity_factor", x, spec[count_spectra]["intensity_factor"][x][0]])
                elif control[j].startswith("path_length ="):
                    spec[count_spectra]["path_length"], truth = parse_fit(get_parameter(j, control[j]))
                    if truth:
                        fitted_par.append([count_spectra, "path_length", spec[count_spectra]["path_length"]])
                elif control[j].startswith("total_pressure ="):
                    spec[count_spectra]["total_pressure"], truth = parse_fit(get_parameter(j, control[j]))
                    if truth:
                        fitted_par.append([count_spectra, "total_pressure", spec[count_spectra]["total_pressure"]])
                elif control[j].startswith("temperature ="):
                    spec[count_spectra]["temperature"], truth = parse_fit(get_parameter(j, control[j]))
                    if truth:
                        fitted_par.append([count_spectra, "temperature", spec[count_spectra]["temperature"]])
                elif control[j].startswith("include_ils ="):
                    spec[count_spectra]["include_ils"] = get_parameter(j, control[j]) == "yes"
                elif control[j].startswith("ils_type ="):
                    spec[count_spectra]["ils_type"] = (get_parameter(j, control[j])).split()
                elif control[j].startswith("snr_ils ="):
                    spec[count_spectra]["snr_ils"] = float(get_parameter(j, control[j]))
                elif control[j].startswith("downsample_factor ="):
                    spec[count_spectra]["downsample_factor"] = int(get_parameter(j, control[j]))
            elif blocpos[i]["type"] == "linelist":
                if control[j].startswith("include ="):
                    linelist[count_linelists]["include"] = get_parameter(j, control[j])
                elif control[j].startswith("linelist_in ="):
                    linelist[count_linelists]["linelist_in"] = get_parameter(j, control[j])
                elif control[j].startswith("linelist_out ="):
                    linelist[count_linelists]["linelist_out"] = get_parameter(j, control[j])
                elif control[j].startswith("off_diagonal_in ="):
                    if len((control[j].split("%")[0]).split("=")) == 1:
                        linelist[count_linelists]["off_diagonal_in"] = ""
                    else:
                        linelist[count_linelists]["off_diagonal_in"] = get_parameter(j, control[j])
                elif control[j].startswith("off_diagonal_out ="):
                    if len((control[j].split("%")[0]).split("=")) == 1:
                        linelist[count_linelists]["off_diagonal_out"] = ""
                    else:
                        linelist[count_linelists]["off_diagonal_out"] = get_parameter(j, control[j])
                elif control[j].startswith("tips ="):
                    linelist[count_linelists]["tips"] = get_parameter(j, control[j])
            elif blocpos[i]["type"] == "calculation":
                if control[j].startswith("spectra_considered ="):
                    if control[j].split("=")[1] == "all":
                        calculation["spectra_considered"] = "all"
                    else:
                        calculation["spectra_considered"] = [(int(i) - 1) for i in get_parameter(j, control[j]).split()]
                elif control[j].startswith("range ="):
                    calculation["range"] = [float(i) for i in get_parameter(j, control[j]).split()]
                elif control[j].startswith("line_profile ="):
                    calculation["line_profile"] = get_parameter(j, control[j])
                elif control[j].startswith("method ="):
                    calculation["method"] = get_parameter(j, control[j])
                elif control[j].startswith("x_calibration_factor ="):
                    calculation["x_calibration_factor"], truth = parse_fit(get_parameter(j, control[j]))
                    if truth:
                        fitted_par.append((None, "x_calibration_factor", calculation["x_calibration_factor"]))
                elif control[j].startswith("output_path ="):
                    calculation["output_path"] = get_parameter(j, control[j])
                    calculation["log_file"] = get_parameter(j, control[j])  # Keeps the original output path
            elif blocpos[i]["type"] == "fit":
                if control[j].startswith("chi2_minimum_threshold ="):
                    fit["chi2_minimum_threshold"] = float(get_parameter(j, control[j]))
                elif control[j].startswith("epsilon ="):
                    fit["epsilon"] = get_parameter(j, control[j])
                elif control[j].startswith("max_iter ="):
                    fit["max_iter"] = get_parameter(j, control[j])

    # Remove spectra not considered spectra
    spec_list = [i for i in range(0, len(spec), 1) if i in calculation["spectra_considered"]]
    spec_dict = {spec_list[i] : i for i in range(0, len(spec_list), 1)}
    if calculation["spectra_considered"] != "all":
        for i in range(len(spec)-1, -1, -1):
            if not i in calculation["spectra_considered"]:
                del spec[i]
                for j in range(len(fitted_par)-1, -1, -1):
                    if fitted_par[j][0] == i:
                        del fitted_par[j]

    for i in range(0, len(fitted_par), 1):
        if type(fitted_par[i][0]) is int:
            fitted_par[i][0] = spec_dict[fitted_par[i][0]]
        fitted_par[i] = tuple(fitted_par[i])


    spectral_data = {"spectra": spec, "linelists": linelist, "calculation": calculation, "fit": fit}

    for i in linelist:
        if i["off_diagonal_in"] != "" and i["off_diagonal_out"] == "":
             print("ERROR: \"off_diagonal_out\" not set for \"{}\"".format(i["off_diagonal_in"]))
             sys.exit()

    for x, i in enumerate(spec):
        linelists_list = {j : False for j in range(1, len(linelist)+1, 1)}
        for j in i["mole_fraction"]:
            for key in linelists_list:
                if key in j[1]:
                    linelists_list[key] = True
        for key in linelists_list:
            if not linelists_list[key]:
                print("WARNING: mole_fraction not set for linelist #{:d} in spectrum #{:d}".format(key, x + 1))
        for j in i["intensity_factor"]:
            for key in linelists_list:
                if key in j[1]:
                    linelists_list[key] = True


    return spectral_data, fitted_par



def get_parameter(i, line):
    if "%" in line:
        line = line.split("%")[0]
    list = line.split('=')
    if len(list) != 2:
        print("\nERROR: incomplete field '{}' at line no. {}".format(line, i + 1))
        sys.exit()
    parameter = list[1].rstrip('\n').strip()
    return(parameter)


def parse_fit(par_str):
    par_str = par_str.strip()
    fit = (par_str[-1] == "f")
    par = float(par_str[:-1])
    return par, fit

