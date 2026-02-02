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
import os


# Change line to fit parameter
def set_line_par(line : str, parpos, par):
    tup_pos = [0, 2, 1]
    for i in range(2, len(parpos) - 3, 3):
        dot_pos = parpos[i][1].find(".")
        if dot_pos == -1:   # Get number of decimals from input line
            dot_pos = line.find(".", parpos[i - 2], parpos[i - 1])
            e_pos = parpos[i][-1]
            if parpos[i][-1] in ["e", "E"]:
                e_pos = line.find(parpos[i][-1], dot_pos, parpos[i - 1])
            dec = e_pos - dot_pos - 1
            if dot_pos == -1:
                dec = 0
            line = line[:parpos[i - 2]] + "{{:{:d}.{:d}{:1}}}".format(int(parpos[i][:-1]), dec, parpos[i][-1]).format(par[tup_pos[int(i/3)]]) + line[parpos[i - 1]:]
        else:
            line = line[:parpos[i - 2]] + "{{:{}}}".format(parpos[i]).format(par[tup_pos[int(i/3)]]) + line[parpos[i - 1]:]
    return line


def write_linelists(linelist_inputs, linelist_outputs, linelists, spec_range, linelist_formats):

    for i in range(0, len(linelist_inputs), 1):

        # Linelist is separated in blocks and needs to be flattened back
        linelist = dict()
        if len(linelists[i]) != 0:
            for key, item in linelists[i][0].items():
                linelist[key] = []

            for j in range(0, len(linelists[i]), 1):
                for key, item in linelists[i][j].items():
                    linelist[key] += item

            order = np.argsort(np.array(linelist["line_position"]), kind = "stable")

            for key, item in linelist.items():
                linelist[key] = np.array(item)[order]
                if linelist[key].ndim == 2:     # Numpy changes tuple into array and convert booleans to numbers
                    linelist[key] = linelist[key].tolist()
                    for j in linelist[key]:
                        j[1] = int(j[1])    # Make sure that converted booleans are integers

    
        with open(linelist_inputs[i], "r") as f, open("{}.temp".format(linelist_outputs[i]), "w") as nf:
            count = 0
            for x, line in enumerate(f):
                if linelist["line_position"][0] <= x <= linelist["line_position"][-1]:
                    for key in linelist_formats[i].keys():
                        if len(linelist_formats[i][key]) == 9 and linelist[key][count][1]:
                            line = set_line_par(line, linelist_formats[i][key], linelist[key][count])
                    count += 1

                nf.write(line)

        os.replace("{}.temp".format(linelist_outputs[i]), linelist_outputs[i])

    return


def write_offdiags(offdiag_inputs, offdiag_outputs, offdiags, offdiag_format):

    for i in range(len(offdiags)-1, -1, -1):
        if len(offdiag_format[i]) != 0:

            # Offdiag is separated in blocks and needs to be flattened back
            offdiag = dict()
            if len(offdiags[i]) != 0:
                for key, item in offdiags[i][0].items():
                    offdiag[key] = []

                for j in range(0, len(offdiags[i]), 1):
                    for key, item in offdiags[i][j].items():
                        offdiag[key] += item


            with open(offdiag_inputs[i]) as f, open("{}.temp".format(offdiag_outputs[i]), "w") as nf:
                count = 0
                for line in f:
                    if "names" in offdiag:
                        if "{} {}".format(line[offdiag_format[i]["name_1"][0] : offdiag_format[i]["name_1"][1]], line[offdiag_format[i]["name_2"][0] : offdiag_format[i]["name_2"][1]]) in offdiag["names"]:
                            if offdiag["line-mixing"][count][1]:
                                line = set_line_par(line, offdiag_format[i]["line-mixing"], offdiag["line-mixing"][count])
                            count += 1
                    nf.write(line)

            os.replace("{}.temp".format(offdiag_outputs[i]), offdiag_outputs[i])

    return




