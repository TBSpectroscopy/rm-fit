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
import os

def write_linelists(linelist_inputs, linelist_outputs, linelists, spec_range):

    for i in range(0, len(linelist_inputs), 1):

        linelist = dict()
        if len(linelists[i]) != 0:
            for key, item in linelists[i][0].items():
                linelist[key] = []

            for j in range(0, len(linelists[i]), 1):
                for key, item in linelists[i][j].items():
                    linelist[key] += item

            order = np.argsort(np.array(linelist["wavenumber"])[:,0], kind = "stable")

            for key, item in linelist.items():
                linelist[key] = np.array(item)[order]
                if linelist[key].ndim == 2:
                    linelist[key] = linelist[key].tolist()
                    for j in linelist[key]:
                        j[1] = int(j[1])

    
        with open(linelist_inputs[i], "r") as f, open("{}.temp".format(linelist_outputs[i]), "w") as nf:
            start = False
            count = 0
            for line in f:
                if line.startswith("$LINE_PARAMETERS"):
                    start = True
                    nf.write(line)
                elif start and not line.startswith("%") and not line.startswith("F"):
                    if float(line[13:26]) < spec_range[0] or float(line[13:26]) > spec_range[1]:
                        nf.write(line)
                    else:
                        nf.write("{:1} {:11} {:12.6f}{:8.1E}{:1}{:11.4E}{:8.1E}{:1}{:10.3E}{:8.1E}{:1}{:10.3E}{:8.1E}{:1}{:10.3E}{:8.1E}{:1}{:10.3E}{:8.1E}{:1}{:10.3E}{:8.1E}{:1}{:10.3E}{:8.1E}{:1}{:10.3E}{:8.1E}{:1}{:9.1f}{:13.5f}{}".format(line[0], linelist["name"][count], linelist["wavenumber"][count][0], linelist["wavenumber"][count][2], "*" * linelist["wavenumber"][count][1], linelist["intensity"][count][0], linelist["intensity"][count][2], "*" * linelist["intensity"][count][1], linelist["foreign-broadening"][count][0], linelist["foreign-broadening"][count][2], "*" * linelist["foreign-broadening"][count][1], linelist["self-broadening"][count][0], linelist["self-broadening"][count][2], "*" * linelist["self-broadening"][count][1], linelist["SD-broadening"][count][0], linelist["SD-broadening"][count][2], "*" * linelist["SD-broadening"][count][1], linelist["narrowing"][count][0], linelist["narrowing"][count][2], "*" * linelist["narrowing"][count][1], linelist["foreign-shift"][count][0], linelist["foreign-shift"][count][2], "*" * linelist["foreign-shift"][count][1], linelist["self-shift"][count][0], linelist["self-shift"][count][2], "*" * linelist["self-shift"][count][1], linelist["SD-shift"][count][0], linelist["SD-shift"][count][2], "*" * linelist["SD-shift"][count][1], linelist["statistical weight"][count], linelist["energy"][count], line[210:]))
                        count += 1
                else:
                    nf.write(line)
        os.replace("{}.temp".format(linelist_outputs[i]), linelist_outputs[i])


    return


def write_offdiags(offdiag_inputs, offdiag_outputs, offdiags):

    for i in range(len(offdiags)-1, -1, -1):

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
                if not line[0:23] in offdiag["names"]:
                    nf.write(line)
                else:
                    nf.write("{:23} {:10.3E}{:8.1E}{:1}\n".format(offdiag["names"][count], offdiag["value"][count][0], offdiag["value"][count][2], "*" * offdiag["value"][count][1]))
                    count += 1

        os.replace("{}.temp".format(offdiag_outputs[i]), offdiag_outputs[i])

    return




