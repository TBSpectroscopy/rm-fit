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
import os
import struct
import numpy as np
import math


def readSpectrum(spec_file, nu_range):
    with open(spec_file, mode = "rb") as f:
        sign = f.read(4)
    if sign == b"\x0A\x0A\xFE\xFE":
        return readOpus(spec_file, nu_range)
    else:
        with open(spec_file) as f:
            split = f.readline().split()
            try:
                float(split[0])
                float(split[1])
            except:
                print("ERROR: \"{}\" unrecognized spectrum file type\nOnly Opus and text files are allowed\n\nText files must contain 1 column for the wavenumber, and one for the signal. The columns must be separated by spaces/tab's".format(spec_file))
                sys.exit(1)

        return readAscii(spec_file, nu_range)



def readOpus(opus_file, nu_range):
    offset = 24
    with open(opus_file, mode = "rb") as f:
        f.seek(offset, 1)
        header = f.read(504 - offset)

        dat_pos = 0
        dat_size = 0

        par_pos = 0
        par_size = 0

        for i in range(0, len(header), 12):
            data_type = struct.unpack("<B", header[i:i+1])[0]
            if data_type in [7, 11, 15]:
                dat_pos = struct.unpack("<I", header[i + 8 : i + 12])[0]
                dat_size = struct.unpack("<I", header[i + 4 : i + 8])[0]
            elif data_type in [23, 31]:
                par_pos = struct.unpack("<I", header[i + 8 : i + 12])[0]
                par_size = struct.unpack("<I", header[i + 4 : i + 8])[0]

        del header
        f.seek(par_pos, 0)
        par = f.read(par_size * 4)
        npt_pos = 0
        fxv_pos = 0
        lxv_pos = 0
        for i in range(0, len(par), 4):
            if par[i: i+4] == b"NPT\x00":
                npt_pos = i
            elif par[i: i+4] == b"FXV\x00":
                fxv_pos = i
            elif par[i: i+4] == b"LXV\x00":
                lxv_pos = i

        n = struct.unpack("=1I", par[npt_pos + 8:npt_pos + 12])[0]

        fxv = struct.unpack("=1d", par[fxv_pos + 8:fxv_pos + 16])[0]
        lxv = struct.unpack("=1d", par[lxv_pos + 8:lxv_pos + 16])[0]
        if nu_range[0] > lxv or nu_range[1] < fxv:
            return []

        inter = (lxv - fxv) / (n - 1)
        new_nu_range = [0, 0]
        nu_pos = [math.ceil((nu_range[0] - fxv)/inter) * (nu_range[0] >= fxv)]
        nu_pos.append((math.floor((nu_range[1] - fxv) / inter) * (nu_range[1] <= lxv)) + (n * (nu_range[1] > lxv)))
        new_nu_range[0] = fxv + (inter * nu_pos[0])
        new_nu_range[1] = fxv + (inter * nu_pos[1])
        n = nu_pos[1] - nu_pos[0]

        del par
        f.seek(dat_pos + (nu_pos[0] * 4))
        dat = f.read(n * 4)
    
    y = np.array(struct.unpack("={:d}f".format(n), dat))
    x = np.linspace(new_nu_range[0], new_nu_range[1], n)

    arr = [x, y]

    return arr



def readAscii(ascii_file, nu_range):
    x = []
    y = []
    with open(ascii_file) as f:
        for line in f:
            split = line.split()
            if nu_range[0] <= float(split[0]) <= nu_range[1]:
                x.append(float(split[0]))
                y.append(float(split[1]))

    return [np.array(x), np.array(y)]

