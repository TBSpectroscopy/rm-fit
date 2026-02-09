#-------------------------------------------------------------------------------------------------
#
# Atomic and Molecular Physics
# Center for Astrophysics | Harvard & Smithsonian
# 60 Garden St., 02138 MA, USA
# E-mail: thibault.bertin@cfa.harvard.edu
#
#-------------------------------------------------------------------------------------------------


import numpy as np
import scipy.constants as constants

def doppler(nu, nu0, T, M):
    nu = np.array(nu)
    gamma = np.sqrt(np.log(2)) * (np.sqrt(2 * constants.Boltzmann * T / (M / (1000*constants.Avogadro))) / constants.c) * (nu0)
    shape = (np.sqrt(np.log(2) / np.pi) / gamma) * np.exp( - np.log(2) * (((nu - nu0)/gamma)**2))   
    return shape

def hwhm(nu0, T, M):
    gamma = (np.sqrt(2 * np.log(2) * constants.Boltzmann * T * constants.Avogadro / (M / (1000))) / constants.c) * (nu0)
    return gamma

def calc_alpha(spectrum_data, mass, linelist, linelist_index, tips, nu):

    mole_fraction = 1.0
    for i in range(0, len(spectrum_data["mole_fraction"]), 1):
        if (linelist_index + 1) in spectrum_data["mole_fraction"][i][2]:
            mole_fraction = spectrum_data["mole_fraction"][i][0]

    alpha = np.zeros(len(nu), dtype = np.double)
    for i in range(0, len(linelist["wavenumber"]), 1):
        alpha += doppler(nu, linelist["wavenumber"][i][0], spectrum_data["temperature"][0], mass) * (spectrum_data["total_pressure"][0] / 1013.25) * linelist["intensity"][i][0] * mole_fraction * constants.physical_constants["Loschmidt constant (273.15 K, 101.325 kPa)"][0] * 1E-6 * (273.15 / spectrum_data["temperature"][0])

    return alpha


