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
import sys
import os


def set_default_par(par, profile):
    prof_comb = {"voigt": ["lorentz"], "rautian": ["lorentz", "narrowing"], "qsd_voigt": ["lorentz", "speed-dependence"], "qsd_rautian": ["lorentz", "speed-dependence", "narrowing"]}
    par_prof = {"essentials": ["wavenumber", "intensity"], "general": ["name", "statistical weight", "energy"], "lorentz": ["self-broadening", "foreign-broadening", "self-shift", "foreign-shift", "line-mixing_fo", "self-broadening_temp", "foreign-broadening_temp"], "narrowing": ["narrowing"], "speed-dependence": ["SD-broadening", "SD-shift"]}

    for i in par_prof["essentials"]:
        if not i in par:
            print("ERROR: column \"{}\" not provided".format(i))
            sys.exit()

    if "name" in par and (not "statistical weight" in par or not "energy" in par):
        print("ERROR: \"statistical weight\" and \"energy\" columns must be provided if \"name\" column exists\nEither remove \"name\" or provide the missing column(s)")
        sys.exit()
    elif not "name" in par:
        par["name"] = ["" for k in par["wavenumber"]]

    for i in prof_comb[profile]:
        for j in par_prof[i]:
            if not j in par:
                par[j] = [(0.0, False, 0.0) for k in par["wavenumber"]]

    return


def get_lines(linelist_file, spec_range, profile):

    mass = 0.0
    mass_perturbing = 0.0
    parpos = dict()
    par = dict()


    startx = 0
    with open(linelist_file) as f:
        for x, line in enumerate(f):
            if line.startswith("$LINE_PARAMETERS"):
                startx = x+1
                break

    if startx == 0:
        print("ERROR: linelist \"{}\" does not contain \"$LINE_PARAMETERS\" block")
        sys.exit()
    with open(linelist_file) as f:
        for x, line in enumerate(f):
            if x < startx:
                if "Mass.....................:" in line:
                    mass = float(line.split()[1])
                elif "Mass Perturbing..........:" in line:
                    mass_perturbing = float(line.split()[2])
                elif "Format...................:" in line:
                    form = line.split("..:")[1]
                    form_st = form.find("\"")
                    form_en = form.find("\"", form_st + 1)
                    form = form[form_st + 1 : form_en]
                    parpos = get_format(form)

                    par = {key : [] for key in parpos.keys() if key != "fit"}
                    par["line_position"] = []
            elif x >= startx:
                if not (line.strip()).startswith("F") and not (line.strip()).startswith("%") and not line.strip() == "":
                    if spec_range[0] < float(line[parpos["wavenumber"][0] : parpos["wavenumber"][1]]) < spec_range[1]:
                        if parpos == dict():
                            print("ERROR: Format not provided in {}".format(linelist_file))
                            sys.exit()
                        read_linelist(par, line, parpos, x)
    if mass <= 0.0:
        print("ERROR: absorber mass not provided in {}".format(linelist_file))
        sys.exit()

    set_default_par(par, profile)

    return par, mass, mass_perturbing, parpos


def get_format(form, name_width = None):   # Convert "Format" field string to absolute positions
    form = form.rstrip()
    parpos = dict()
    end = 0
    pos = 0

    count = 1

    while end < len(form) - 1:
        start = form.find("{", end)
        space = start - end - (end != 0) # Number of columns between two parameters
        end = form.find("}", start)
        par = form[start + 1:end]

        split = par.split(":")
        if not name_width or split[0] != "name":
            var_type = split[1]
            try:
                width = int(split[1]) # Number of columns for the parameter
            except:
                width = int(split[1][:-1]) # Number of columns for the parameter
        else:
            split[0] += "_{:d}".format(count)
            count += 1
            width = name_width
            var_type = "{:d}".format(name_width)
        parpos[split[0]] = [pos + space, pos + width + space, var_type]

        pos += width + space    # Position of the last character
    parname = list(set([i.removesuffix("_unc").removesuffix("_fit") for i in parpos.keys()]))
    for key in parname:
        if "{}_fit".format(key) in parpos:
            try:
                parpos[key] += parpos["{}_unc".format(key)] + parpos["{}_fit".format(key)]
                del parpos["{}_unc".format(key)]
                del parpos["{}_fit".format(key)]
            except:
                print("ERROR: No uncertainty column was found for fit parameter {}\n".format(key))
                sys.exit()
        elif "{}_unc".format(key) in parpos and not "{}_fit".format(key) in parpos:
            parpos[key] += parpos["{}_unc".format(key)]
            del parpos["{}_unc".format(key)]

    return parpos



def read_linelist(par, line, parpos, x):

    fit = True
    if "fit" in parpos:
        fit = (line[parpos["fit"][0] : parpos["fit"][1]] == "*")

    for key in par.keys():
        if key != "line_position":
            if key in ["wavenumber", "intensity", "self-broadening", "foreign-broadening", "self-shift", "foreign-shift", "line-mixing_fo", "narrowing", "SD-broadening", "SD-shift", "foreign-broadening_temp", "self-broadening_temp"]:
                if len(parpos[key]) == 3:
                    par[key].append((float(line[parpos[key][0] : parpos[key][1]]), False, 0.0))
                elif len(parpos[key]) == 6:
                    par[key].append((float(line[parpos[key][0] : parpos[key][1]]), False, float(line[parpos[key][3] : parpos[key][4]])))
                else:
                    par[key].append((float(line[parpos[key][0] : parpos[key][1]]), (line[parpos[key][6] : parpos[key][7]] == "*") and fit, float(line[parpos[key][3] : parpos[key][4]])))
            else:
                try:
                    par[key].append(float(line[parpos[key][0] : parpos[key][1]]))
                except:
                    par[key].append(line[parpos[key][0] : parpos[key][1]])
    par["line_position"].append(x)
    return


def get_offdiag(offdiag_list, names):
    offdiag = {"names": [], "line-mixing": []}
    names_used = []
    parpos = dict()

    if os.path.exists(offdiag_list):
        for name in names:
            if len(name.strip()) != 0:
                names_used.append(name)
        with open(offdiag_list) as f:
            start = False
            for line in f:
                if start:
                    for name in names_used:
                        if len(parpos) != 0:
                            if name == line[parpos["name_1"][0] : parpos["name_1"][1]]:
                                offdiag["names"].append("{} {}".format(line[parpos["name_1"][0] : parpos["name_1"][1]], line[parpos["name_2"][0] : parpos["name_2"][1]]))
                                offdiag["line-mixing"].append((float(line[parpos["line-mixing"][0] : parpos["line-mixing"][1]]), line[parpos["line-mixing"][6] : parpos["line-mixing"][7]] == "*", float(line[parpos["line-mixing"][3] : parpos["line-mixing"][4]])))
                elif "Format...................:" in line and len(names_used) != 0:
                    form = line.split("..:")[1]
                    form_st = form.find("\"")
                    form_en = form.find("\"", form_st + 1)
                    form = form[form_st + 1 : form_en]
                    parpos = get_format(form, len(names[0]))
                    if not "name_1" in parpos or not "name_2" in parpos:
                        parpos = dict()
                        print("WARNING: \"name_1\" and \"name_2\" columns must be provided for the off-diagonal elements to be used")
                elif line.startswith("$OFF-DIAGONAL_PARAMETERS"):
                    start = True

    return offdiag, parpos


def separate_lines(offdiag, linelist, offdiag_parpos):
    names = [offdiag["names"][i] for i in range(len(offdiag["names"]) - 1, -1, -1)]
    value = [offdiag["line-mixing"][i] for i in range(len(offdiag["line-mixing"]) - 1, -1, -1)]
    namepos = []
    if len(names) != 0:
        namepos = [0, offdiag_parpos["name_1"][1] - offdiag_parpos["name_1"][0], offdiag_parpos["name_1"][1] - offdiag_parpos["name_1"][0] + 1, (2 * (offdiag_parpos["name_1"][1] - offdiag_parpos["name_1"][0])) + 1]
    blocknames = []
    offdiag_blocks = []
    linelist_blocks = []
    added = False
    while len(names) != 0:
        if not added: # If no line is matched, all the combinations of the block were found
            blocknames.append([names[-1][namepos[0] : namepos[1]], names[-1][namepos[2] : namepos[3]]])
            offdiag_blocks.append({"names": [], "line-mixing": []})
        added = False
        for i in range(len(names)-1, -1, -1):
            if names[i][namepos[0] : namepos[1]] in blocknames[-1]:
                offdiag_blocks[-1]["names"].append(names[i])
                offdiag_blocks[-1]["line-mixing"].append(value[i])
                if not names[i][namepos[2] : namepos[3]] in blocknames[-1]:
                    blocknames[-1].append(names[i][namepos[2] : namepos[3]])
                    added = True
                del names[i]
                del value[i]
            elif names[i][namepos[2] : namepos[3]] in blocknames[-1]:
                offdiag_blocks[-1]["names"].append(names[i])
                offdiag_blocks[-1]["line-mixing"].append(value[i])
                if not names[i][namepos[0] : namepos[1]] in blocknames[-1]:
                    blocknames[-1].append(names[i][namepos[0] : namepos[1]])
                    added = True
                del names[i]
                del value[i]

    linelist_temp = {key:value for key, value in linelist.items()}
    for i in range(0, len(blocknames), 1):
        linelist_blocks.append({key:[] for key, value in linelist_temp.items()})
        for j in range(len(linelist_temp["name"])-1, -1, -1):
            if linelist_temp["name"][j] in blocknames[i]:
                for key, value in linelist_temp.items():
                    linelist_blocks[i][key].append(value[j])
                    del linelist_temp[key][j]
        for key, value in linelist_blocks[i].items():
            linelist_blocks[i][key].reverse()
    if len(linelist_temp["name"]) != 0:
        linelist_blocks.append(linelist_temp)
    
    return linelist_blocks, offdiag_blocks


def get_blocks(spectral_data):

    linelists = []
    offdiags = []
    blocks = [[] for i in range(0, len(spectral_data["linelists"]), 1)]
    for i in range(0, len(spectral_data["linelists"]), 1):
        linelist_temp, spectral_data["linelists"][i]["mass"], spectral_data["linelists"][i]["mass_perturbing"], spectral_data["linelists"][i]["format"] = get_lines(spectral_data["linelists"][i]["linelist_in"], spectral_data["calculation"]["range"], spectral_data["calculation"]["line_profile"])
        offdiag_temp, spectral_data["linelists"][i]["format_offdiag"] = get_offdiag(spectral_data["linelists"][i]["off_diagonal_in"], linelist_temp["name"])

        linelist_blocks, offdiag_blocks = separate_lines(offdiag_temp, linelist_temp, spectral_data["linelists"][i]["format_offdiag"])
        linelists.append(linelist_blocks)
        offdiags.append(offdiag_blocks)

    return linelists, offdiags


def get_fitted_parameters(linelists, offdiags):
    fitted_par = []
    for ii, i in enumerate(linelists): # loop over files
        for jj, j in enumerate(i): # loop over blocks
            for key, value in j.items(): # loop over parameters
                if not key in ["name", "statistical weight", "energy", "line_position"]:
                    for kk, k in enumerate(value): # loop over lines
                        if type(k) is tuple and k[1]:
                            fitted_par.append((ii, jj, kk, key, k[0]))

    for ii, i in enumerate(offdiags): # loop over files
        for jj, j in enumerate(i): # loop over blocks
            for key, value in j.items():
                if key != "names":
                    for kk, k in enumerate(value):
                        if type(k) is tuple and k[1]:
                            fitted_par.append((ii, jj, kk, key, k[0]))
    return fitted_par


