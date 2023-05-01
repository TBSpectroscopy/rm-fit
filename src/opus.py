#-------------------------------------------------------------------------------------------------
#
# Thibault BERTIN
# Spectroscopy, Quantum Chemistry and Atmospheric Remote Sensing (SQUARES), C.P. 160/09
# Universite Libre de Bruxelles
# 50 avenue F. D. Roosevelt, B-1050 Brussels, Belgium
# Phone: +32.2.650.24.18 - E-mail: thibault.bertin@ulb.be - Web: http://www.ulb.ac.be/cpm
#
#-------------------------------------------------------------------------------------------------


import sys
import os
import struct
import numpy as np
import math


def readOpus(opus_file, nu_range):
    f = open(opus_file, mode = "rb")
    opus = f.read()
    npt_pos = opus.find(b"NPT\x00\x00\x00")
    end_pos = []

    # Append all END positions
    end_char = b"END\x00\x00\x00"
    offset = 0
    while offset < len(opus):
        offset = opus.find(end_char, offset)
        if offset == -1:
            break
        end_pos.append(offset)
        offset += len(end_char)
    del offset

    fxv_pos = npt_pos + 12
    lxv_pos = fxv_pos + 16
    n = struct.unpack("=1I", opus[npt_pos + 8:npt_pos + 12])[0]

    dat_pos = 0
    for i in range(0, len(end_pos)-1, 1):
        if end_pos[i+1] - end_pos[i] > 4*n:
            dat_pos = end_pos[i] + 8

    fxv = struct.unpack("=1d", opus[fxv_pos + 8:fxv_pos + 16])[0]
    lxv = struct.unpack("=1d", opus[lxv_pos + 8:lxv_pos + 16])[0]
    if nu_range[0] < fxv or nu_range[1] > lxv:
        return []

    inter = (lxv - fxv) / (n - 1)
    new_nu_range = [0, 0]
    a = math.ceil((nu_range[0] - fxv)/inter)
    new_nu_range[0] = fxv + (inter * a)
    new_nu_range[1] = new_nu_range[0] + (inter * math.floor((nu_range[1] - new_nu_range[0]) / inter))
    n = round(((new_nu_range[1] - new_nu_range[0]) / (inter)) + 1)
    
    y = np.array(struct.unpack("={:d}f".format(n), opus[dat_pos + (a*4):dat_pos + (a*4) + (n*4)]))
    x = np.linspace(new_nu_range[0], new_nu_range[1], n)

    arr = [x, y]

    f.close()

    return arr


