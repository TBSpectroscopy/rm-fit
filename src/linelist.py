#-------------------------------------------------------------------------------------------------
#
# Thibault BERTIN
# Spectroscopy, Quantum Chemistry and Atmospheric Remote Sensing (SQUARES), C.P. 160/09
# Universite Libre de Bruxelles
# 50 avenue F. D. Roosevelt, B-1050 Brussels, Belgium
# Phone: +32.2.650.24.18 - E-mail: thibault.bertin@ulb.be - Web: http://www.ulb.ac.be/cpm
#
#-------------------------------------------------------------------------------------------------


import numpy as np
import sys


def get_lines(linelist_file, spec_range):

    par = {\
            "name": [],\
            "wavenumber": [],\
            "intensity": [],\
            "foreign-broadening": [],\
            "self-broadening": [],\
            "SD-broadening": [],\
            "narrowing": [],\
            "foreign-shift": [],\
            "self-shift": [],\
            "SD-shift": [],\
            "statistical weight": [],\
            "energy": []\
            }

    mass = 0.0
    mass_perturbing = 0.0

    startx = 0
    with open(linelist_file) as f:
        for x, line in enumerate(f):
            if line.startswith("$LINE_PARAMETERS"):
                startx = x+2
                break

    if startx == 0:
        print("ERROR: linelist \"{}\" does not contain \"$LINE_PARAMETERS\" block")
        sys.exit()
    with open(linelist_file) as f:
        mix = False
        for x, line in enumerate(f):
            if x < startx:
                if "Mass.....................:" in line:
                    mass = float(line.split()[1])
                elif "Sample...................:" in line:
                    mix = (line.split()[1] == "mixture")
                elif "Mass Perturbing..........:" in line:
                    mass_perturbing = float(line.split()[2])
            elif x >= startx:
                if spec_range[0] < float(line[13:26]) < spec_range[1]:
                    read_linelist(par, line, mix)
    if mass <= 0.0:
        print("ERROR: absorber mass not provided")
        sys.exit()
    if mix and mass_perturbing <= 0.0:
        print("ERROR: perturber mass not provided")
        sys.exit()

    return par, mass, mass_perturbing


def read_linelist(par, line, mix):
    fit = (line[0] == "*")
    if mix:
        fp = 2
        lp = 13
        par["name"].append(line[fp:lp])
        fp = lp + 1
        lp = fp + 12
        par["wavenumber"].append((float(line[fp:lp]), (line[lp+8] == "*") and fit, float(line[lp+1:lp+8])))
        fp = lp + 9
        lp = fp + 11
        par["intensity"].append((float(line[fp:lp]), line[lp+8] == "*" and fit, float(line[lp+1:lp+8])))
        fp = lp + 9
        lp = fp + 10
        par["foreign-broadening"].append((float(line[fp:lp]), line[lp+8] == "*" and fit, float(line[lp+1:lp+8])))
        fp = lp + 9
        lp = fp + 10
        par["self-broadening"].append((float(line[fp:lp]), line[lp+8] == "*" and fit, float(line[lp+1:lp+8])))
        fp = lp + 9
        lp = fp + 10
        par["SD-broadening"].append((float(line[fp:lp]), line[lp+8] == "*" and fit, float(line[lp+1:lp+8])))
        fp = lp + 9
        lp = fp + 10
        par["narrowing"].append((float(line[fp:lp]), line[lp+8] == "*" and fit, float(line[lp+1:lp+8])))
        fp = lp + 9
        lp = fp + 10
        par["foreign-shift"].append((float(line[fp:lp]), line[lp+8] == "*" and fit, float(line[lp+1:lp+8])))
        fp = lp + 9
        lp = fp + 10
        par["self-shift"].append((float(line[fp:lp]), line[lp+8] == "*" and fit, float(line[lp+1:lp+8])))
        fp = lp + 9
        lp = fp + 10
        par["SD-shift"].append((float(line[fp:lp]), line[lp+8] == "*" and fit, float(line[lp+1:lp+8])))
        fp = lp + 9
        lp = fp + 9
        par["statistical weight"].append(float(line[fp:lp]))
        fp = lp
        lp = fp + 13
        par["energy"].append(float(line[fp:lp]))

    return


def get_offdiag(offdiag_list, names):
    offdiag = {"names": [], "value": []}
    names_used = []
    for name in names:
        if name != "     ":
            names_used.append(name)
    with open(offdiag_list) as f:
        for line in f:
            for name in names_used:
                if name == line[0:11]:
                    offdiag["names"].append(line[0:23])
                    offdiag["value"].append((float(line[24:34]), line[42] == "*", 0.0))
    return offdiag


def separate_lines(offdiag, linelist):
    names = [offdiag["names"][i] for i in range(len(offdiag["names"]) - 1, -1, -1)]
    value = [offdiag["value"][i] for i in range(len(offdiag["value"]) - 1, -1, -1)]
    blocknames = []
    count = 0
    offdiag_blocks = []
    linelist_blocks = []
    while len(names) != 0:
        if count % 2 == 0:  # One line is associated to every other lines it is coupled with. Therefore a single pass gives all the mixed lines, and a second gives all the combinations.
            blocknames.append([names[-1][0:11], names[-1][12:23]])
            offdiag_blocks.append({"names": [], "value": []})
        for i in range(len(names)-1, -1, -1):
            if names[i][0:11] in blocknames[-1]:
                offdiag_blocks[-1]["names"].append(names[i])
                offdiag_blocks[-1]["value"].append(value[i])
                if not names[i][12:23] in blocknames[-1]:
                    blocknames[-1].append(names[i][12:23])
                del names[i]
                del value[i]
            elif names[i][12:23] in blocknames[-1]:
                offdiag_blocks[-1]["names"].append(names[i])
                offdiag_blocks[-1]["value"].append(value[i])
                if not names[i][0:11] in blocknames[-1]:
                    blocknames[-1].append(names[i][0:11])
                del names[i]
                del value[i]
        count += 1

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
        linelist_temp, spectral_data["linelists"][i]["mass"], spectral_data["linelists"][i]["mass_perturbing"] = get_lines(spectral_data["linelists"][i]["linelist_in"], spectral_data["calculation"]["range"])
        offdiag_temp = (get_offdiag(spectral_data["linelists"][i]["off_diagonal_in"], linelist_temp["name"]))

        linelist_blocks, offdiag_blocks = separate_lines(offdiag_temp, linelist_temp)
        linelists.append(linelist_blocks)
        offdiags.append(offdiag_blocks)

    return linelists, offdiags


def get_fitted_parameters(linelists, offdiags):
    fitted_par = []
    for ii, i in enumerate(linelists): # loop over files
        for jj, j in enumerate(i): # loop over blocks
            for key, value in j.items(): # loop over parameters
                if key != "name" and key != "statistical weight" and key != "energy":
                    for kk, k in enumerate(value): # loop over lines
                        if k[1]:
                            fitted_par.append((ii, jj, kk, key, k[0]))

    for ii, i in enumerate(offdiags): # loop over files
        for jj, j in enumerate(i): # loop over blocks
            for key, value in j.items():
                if key != "names":
                    for kk, k in enumerate(value):
                        if k[1]:
                            fitted_par.append((ii, jj, kk, key, k[0]))
    return fitted_par


