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
import racef


# x = sqrt(ln(2))(nu-nu0)/yD
# yD = sqrt(ln(2))y'D
# y'D = Doppler half width at max/e
# y = sqrt(ln(2))yL/yD
# z = x + iy

def calc_rautian(nu, nu0, doppler_width, lorentz_width, narrowing, shift, lm):
    x = np.sqrt(np.log(2)) * (nu - nu0 - shift) / doppler_width
    y = np.sqrt(np.log(2)) * lorentz_width / doppler_width
    b = np.sqrt(np.log(2)) * narrowing / doppler_width
    w = racef.weideman40a(x, (y+b))
    if b:
        shape = (np.sqrt(np.log(2) / np.pi) / doppler_width) * (((1 - (1j * lm)) * w) / (1 - (np.sqrt(np.pi) * b * w))).real
        return shape
    return ((np.sqrt(np.log(2) / np.pi) / doppler_width) * ((1 - (1j * lm)) * w).real )

