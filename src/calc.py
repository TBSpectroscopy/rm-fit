#---------------------------------------------------------------------------------------------------
#
# Thibault BERTIN
# Atomic and Molecular Physics
# Center for Astrophysics | Harvard & Smithsonian
# 60 Garden St., 02138 MA, USA
# E-mail: thibault.bertin@cfa.harvard.edu
#
#
#---------------------------------------------------------------------------------------------------

import numpy as np
import scipy.constants as constants
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import scipy.signal as signal
import sd
import math
import calc_g
import qSDV
import qSDHC
import pCqSDHC
import mHT
import ctrl_pars
import doppler
import rautian
import tips
import sys

#---------------------------------------------------------------------------------------------------
# Calculate limits up to which the line is calculated in the wavenumber axis -----------------------

def calc_lims(lowest_value, doppler_width, lorentz_width, shift, nu0, vnu):

    numaxd = 0
    numaxl = 0
    dcheck = False
    lcheck = False
    if lowest_value < (math.sqrt(np.log(2) / math.pi) / doppler_width): # If lowest line profile point is lower than maximum of Doppler and Lorentzian
        numaxd = doppler_width * math.sqrt( - np.log(lowest_value / (math.sqrt(np.log(2) / math.pi) / doppler_width)) / math.log(2))
        dcheck = True
    if lorentz_width and lowest_value < (1 / (math.pi * lorentz_width)):
        numaxl = math.sqrt((lorentz_width / (math.pi * lowest_value)) - (lorentz_width ** 2))
        lcheck = True

    numax = max(numaxd, numaxl)
    lims = np.searchsorted(vnu, [nu0 + shift - numax, nu0 + shift + numax])

    return lims[0], min(lims[1] + 1, len(vnu)), dcheck + lcheck > 0

#---------------------------------------------------------------------------------------------------

# Get parameters for the different profiles ========================================================
#---------------------------------------------------------------------------------------------------
def get_intensity(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i):
    intensity = linelist["intensity"][i][0] * intensity_factor * tips_ratio * np.exp(-1.4387769 * linelist["energy"][i] / temperature) * (1 - np.exp(-1.4387769 * linelist["wavenumber"][i][0] / temperature)) / (np.exp(-1.4387769 * linelist["energy"][i] / 296.0) * (1 - np.exp(-1.4387769 * linelist["wavenumber"][i][0] / 296.0)))

    return intensity
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i):
    narrowing = (linelist["narrowing"][i][0] * ((296.0/temperature)**linelist["narrowing_temp"][i][0])) * pressure

    return narrowing
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_narrowing_i(pressure, temperature, mass, mole_fraction, linelist, i):
    narrowing_i = (linelist["narrowing_i"][i][0] * ((296.0/temperature)**linelist["narrowing_i_temp"][i][0])) * pressure

    return narrowing_i
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_line_mixing_fo(pressure, temperature, mass, mole_fraction, linelist, i):
    fo_lm = (linelist["line-mixing_fo"][i][0] * ((296.0/temperature)**linelist["line-mixing_fo_temp"][i][0])) * pressure

    return fo_lm
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_line_mixing(pressure, temperature, mass, mole_fraction, offdiag, i):
    lm = (offdiag["line-mixing"][i][0] * ((296.0/temperature)**offdiag["line-mixing_temp"][i][0])) * pressure

    return lm
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_lorentz_par(pressure, temperature, mass, mole_fraction, linelist, i):
    lorentz_width = ((linelist["self-broadening"][i][0] * ((296.0/temperature)**linelist["self-broadening_temp"][i][0]) * mole_fraction) + (linelist["foreign-broadening"][i][0] * ((296.0/temperature)**linelist["foreign-broadening_temp"][i][0]) * (1 - mole_fraction))) * pressure
    lorentz_shift = ((linelist["self-shift"][i][0] * ((296.0/temperature)**linelist["self-shift_temp"][i][0]) * mole_fraction) + (linelist["foreign-shift"][i][0] * ((296.0/temperature)**linelist["foreign-shift_temp"][i][0]) * (1 - mole_fraction))) * pressure
    return lorentz_width, lorentz_shift
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_voigt_par(pressure, temperature, mass, mole_fraction, linelist, i):
    lorentz_width, lorentz_shift = get_lorentz_par(pressure, temperature, mass, mole_fraction, linelist, i)
    doppler_width = doppler.hwhm(linelist["wavenumber"][i][0], temperature, mass)

    return lorentz_width, lorentz_shift, doppler_width
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_rautian_par(pressure, temperature, mass, mole_fraction, linelist, i):
    lorentz_width, lorentz_shift, doppler_width = get_voigt_par(pressure, temperature, mass, mole_fraction, linelist, i)
    narrowing = get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i)

    return lorentz_width, lorentz_shift, doppler_width, narrowing
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_qsd_lorentz_par(pressure, temperature, mass, mole_fraction, linelist, i):
    n_v, u_max = 21, 4.0                      # u = v/vp where vp is the most probable speed
    boltzmann_constant, avogadro_constant = 1.380649E-16, constants.Avogadro
    vp = math.sqrt(2.0*boltzmann_constant * temperature / (mass/avogadro_constant))
    u = np.linspace(0, u_max*vp, n_v, dtype = np.double) / vp

    b_self = sd.sdpar(linelist["self-broadening"][i][0], linelist["self-broadening_temp"][i][0], temperature, u, mass, mass, "quadratic", dpp0 = linelist["SD-broadening"][i][0])
    d_self = sd.sdpar(linelist["self-shift"][i][0], linelist["self-shift_temp"][i][0], temperature, u, mass, mass, "quadratic", dpp0 = linelist["SD-shift"][i][0])
    b_other = sd.sdpar(linelist["foreign-broadening"][i][0], linelist["foreign-broadening_temp"][i][0], temperature, u, mass, mass, "quadratic", dpp0 = linelist["SD-broadening"][i][0])
    d_other = sd.sdpar(linelist["foreign-shift"][i][0], linelist["foreign-shift_temp"][i][0], temperature, u, mass, mass, "quadratic", dpp0 = linelist["SD-shift"][i][0])
    sd_width = ((b_self * mole_fraction) + (b_other *(1 - mole_fraction))) * pressure
    sd_shift = ((d_self * mole_fraction) + (d_other * (1 - mole_fraction))) * pressure

    return sd_width, sd_shift
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_qsd_lorentz_par_only(pressure, temperature, mass, mole_fraction, linelist, i):
    lorentz_width, lorentz_shift = get_lorentz_par(pressure, temperature, mass, mole_fraction, linelist, i)

    sd_width = linelist["SD-broadening"][i][0] * lorentz_width
    sd_shift = linelist["SD-shift"][i][0] * lorentz_shift

    return sd_width, sd_shift
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_qsd_voigt_par(pressure, temperature, mass, mole_fraction, linelist, i):
    sd_width, sd_shift = get_qsd_lorentz_par(pressure, temperature, mass, mole_fraction, linelist, i)
    doppler_width = doppler.hwhm(linelist["wavenumber"][i][0], temperature, mass)

    return sd_width, sd_shift, doppler_width
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_qsd_rautian_par(pressure, temperature, mass, mole_fraction, linelist, i):
    sd_width, sd_shift, doppler_width = get_qsd_rautian_par(pressure, temperature, mass, mole_fraction, linelist, i)
    narrowing = get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i)

    return sd_width, sd_shift, doppler_width, narrowing
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_hartmann_tran_par(pressure, temperature, mass, mole_fraction, linelist, i):
    sd_width, sd_shift, doppler_width = get_qsd_rautian_par(pressure, temperature, mass, mole_fraction, linelist, i)
    narrowing = get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i)
    narrowing_i = get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i)

    return sd_width, sd_shift, doppler_width, narrowing, linelist["correlation_HT"][i][0]
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
def get_mod_hartmann_tran_par(pressure, temperature, mass, mole_fraction, linelist, i):
    sd_width, sd_shift, doppler_width = get_qsd_rautian_par(pressure, temperature, mass, mole_fraction, linelist, i)
    narrowing = get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i)
    narrowing_i = get_narrowing_i(pressure, temperature, mass, mole_fraction, linelist, i)

    return sd_width, sd_shift, doppler_width, narrowing, narrowing_i
#---------------------------------------------------------------------------------------------------
# Get intensity, width, shift, ... -----------------------------------------------------------------
def get_params(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i, get_profile_par : callable):
    intensity = get_intensity(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i)
    fo_lm = get_line_mixing_fo(pressure, temperature, mass, mole_fraction, linelist, i)

    return intensity, get_profile_par(pressure, temperature, mass, mole_fraction, linelist, i), fo_lm
#---------------------------------------------------------------------------------------------------
#===================================================================================================

# Calculate absorption cross section ---------------------------------------------------------------

def calc_alpha(profile, spectrum_data, linelist_data, linelist, linelist_index, vnu, offdiag = dict(), method = ""):
    mass = linelist_data["mass"]
    mass_p = linelist_data["mass_perturbing"]
    pressure = spectrum_data["total_pressure"] / 1013.25

    mole_fraction = 0.0
    for i in spectrum_data["mole_fraction"]:
        if linelist_index + 1 in i[1]:
            mole_fraction = i[0]
    intensity_factor = 1.0
    for i in spectrum_data["intensity_factor"]:
        if linelist_index + 1 in i[1]:
            intensity_factor = i[0]

    if mole_fraction == 0.0 or intensity_factor == 0.0:
        return np.zeros(len(vnu), dtype = np.double)

    lowest_value =  - np.log(1 - spectrum_data["lowest_value"]) / ((constants.physical_constants["Loschmidt constant (273.15 K, 101.325 kPa)"][0] * 1E-6 * (273.15 / spectrum_data["temperature"]) * (mole_fraction * pressure)) * spectrum_data["path_length"])   # - ln(1 - Transmittance) divided by density * path length

    tips_ = tips.get_tips(linelist_data["tips"], spectrum_data["temperature"])
    tips_ratio = tips.get_tips(linelist_data["tips"], 296.0) / tips_

    alpha = np.zeros(len(vnu), dtype = np.double)
    if profile in ["qsd_voigt", "qsd_rautian", "qsd_nelkin-ghatak", "hartmann-tran", "modified_hartmann-tran"]:
        if offdiag == dict() or method == "first-order" or offdiag["line-mixing"] == []:
            alpha = calc_qsdprofile(profile, pressure, spectrum_data["temperature"], mass, mole_fraction, linelist, vnu, lowest_value, tips_ratio, intensity_factor)
        else:
            alpha = calc_sd_mat(profile, pressure, spectrum_data["temperature"], mass, mass_p, mole_fraction, linelist, offdiag, tips_, vnu, lowest_value, method, tips_ratio, intensity_factor)
    elif profile == "rautian":
        if offdiag == dict() or method == "first-order" or offdiag["line-mixing"] == []:
            alpha = calc_rautian(pressure, spectrum_data["temperature"], mass, mole_fraction, linelist, vnu, lowest_value, tips_ratio, intensity_factor)
        else:
            alpha = calc_nosd_mat(pressure, spectrum_data["temperature"], mass, mole_fraction, linelist, offdiag, tips_, vnu, lowest_value, profile, tips_ratio, intensity_factor)
    elif profile == "voigt":
        if offdiag == dict() or method == "first-order" or offdiag["line-mixing"] == []:
            alpha = calc_voigt(pressure, spectrum_data["temperature"], mass, mole_fraction, linelist, vnu, lowest_value, tips_ratio, intensity_factor)
        else:
            alpha = calc_nosd_mat(pressure, spectrum_data["temperature"], mass, mole_fraction, linelist, offdiag, tips_, vnu, lowest_value, profile, tips_ratio, intensity_factor)

    alpha *= constants.physical_constants["Loschmidt constant (273.15 K, 101.325 kPa)"][0] * 1E-6 * (273.15 / spectrum_data["temperature"]) * (mole_fraction * pressure)

    return alpha
#---------------------------------------------------------------------------------------------------

# Calculate Voigt absorption cross section using Weideman algorithm from Py4CAtS -------------------

def calc_voigt(pressure, temperature, mass, mole_fraction, linelist, vnu, lowest_value, tips_ratio, intensity_factor):

    sigma = np.zeros(len(vnu), dtype = np.double)

    dnu = (vnu[-1] - vnu[0]) / (len(vnu) - 1)
    for i in range(0, len(linelist["wavenumber"]), 1):

        intensity, (lorentz_width, lorentz_shift, doppler_width), fo_lm = get_params(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i, get_voigt_par)

        lowest_value_ = lowest_value / intensity # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        lim0, lim1, calc_check = calc_lims(lowest_value_, doppler_width, lorentz_width, lorentz_shift, linelist["wavenumber"][i][0], vnu)

        if calc_check:

            voigt = rautian.calc_rautian(vnu[lim0:lim1], linelist["wavenumber"][i][0], doppler_width, lorentz_width, 0.0, lorentz_shift, fo_lm) * intensity

            sigma[lim0:lim1] += voigt

    return sigma
#---------------------------------------------------------------------------------------------------

# Calculate Rautian absorption cross section using Weideman algorithm from Py4CAtS -----------------

def calc_rautian(pressure, temperature, mass, mole_fraction, linelist, vnu, lowest_value, tips_ratio, intensity_factor):

    sigma = np.zeros(len(vnu), dtype = np.double)

    dnu = (vnu[-1] - vnu[0]) / (len(vnu) - 1)
    for i in range(0, len(linelist["wavenumber"]), 1):
        intensity, (lorentz_width, lorentz_shift, doppler_width, narrowing), fo_lm = get_params(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i, get_rautian_par)

        lowest_value_ = lowest_value / intensity # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        lim0, lim1, calc_check = calc_lims(lowest_value_, doppler_width, lorentz_width, lorentz_shift, linelist["wavenumber"][i][0], vnu)
        if calc_check:

            rau = rautian.calc_rautian(vnu[lim0:lim1], linelist["wavenumber"][i][0], doppler_width, lorentz_width, narrowing, lorentz_shift, fo_lm) * intensity

            sigma[lim0:lim1] += rau

    return sigma
#-------------------------------------------------------------------------------------------------

# Calculate qSDR absorption cross section using H. Tran et al's algorithm, "Efficient computation of some speed-dependent isolated line profiles", JQSRT 129 (2013) 199-203 and the corresponding erratum, JQSRT 134

def calc_qsdprofile(profile : str, pressure, temperature, mass, mole_fraction, linelist, vnu, lowest_value, tips_ratio, intensity_factor):

    sigma = np.zeros(len(vnu), dtype = np.double)



    dnu = (vnu[-1] - vnu[0]) / (len(vnu) - 1)
    for i in range(0, len(linelist["wavenumber"]), 1):
        intensity, (lorentz_width, lorentz_shift, doppler_width), fo_lm = get_params(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i, get_voigt_par)

        lowest_value_ = lowest_value / intensity # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        LS_qSDV_R = 0.0
        LS_qSDV_I = 0.0

        lim0, lim1, calc_check = calc_lims(lowest_value_, doppler_width, lorentz_width, lorentz_shift, linelist["wavenumber"][i][0], vnu)

        if calc_check:
            sd_width, sd_shift = get_qsd_lorentz_par_only(pressure, temperature, mass, mole_fraction, linelist, i)
            if profile == "qsd_voigt":
                line_ri = [qSDV.qsdv(linelist["wavenumber"][i][0], doppler_width, lorentz_width, sd_width, lorentz_shift, sd_shift, vnu[j], LS_qSDV_R, LS_qSDV_I) for j in range(lim0, lim1, 1)]
            if profile in ["qsd_rautian", "qsd_nelkin-ghatak"]:
                narrowing = get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i)
                line_ri = [qSDHC.qsdhc(linelist["wavenumber"][i][0], doppler_width, lorentz_width, sd_width, lorentz_shift, sd_shift, narrowing, vnu[j], LS_qSDV_R, LS_qSDV_I) for j in range(lim0, lim1, 1)]
            elif profile == "modified_hartmann-tran":
                narrowing = get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i)
                narrowing_i = get_narrowing_i(pressure, temperature, mass, mole_fraction, linelist, i)
                narrowing = narrowing + (1j * narrowing_i)
                line_ri = [mHT.qsdhc(linelist["wavenumber"][i][0], doppler_width, lorentz_width, sd_width, lorentz_shift, sd_shift, narrowing, vnu[j], LS_qSDV_R, LS_qSDV_I) for j in range(lim0, lim1, 1)]
            elif profile == "hartmann-tran":
                narrowing = get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i)
                sd_width, sd_shift, doppler_width, narrowing, correlation_HT = get_hartmann_tran_par(pressure, temperature, mass, mole_fraction, linelist, i)
                line_ri = [pCqSDHC.pcqsdhc(linelist["wavenumber"][i][0], doppler_width, lorentz_width, sd_width, lorentz_shift, sd_shift, narrowing, correlation_HT, vnu[j], LS_qSDV_R, LS_qSDV_I) for j in range(lim0, lim1, 1)]

            sigma[lim0:lim1] += np.array([(j[0] - (fo_lm * j[1])) * intensity for j in line_ri])

    return sigma
#---------------------------------------------------------------------------------------------------

# Calculate hard collision absorption cross section using relaxation matrix ------------------------

def calc_nosd_mat(pressure, temperature, mass, mole_fraction, linelist, offdiag, tips, vnu, lowest_value, profile, tips_ratio, intensity_factor):

    c2_constant = 1.438776877 # cm K

    # Diagonal matrices containing the line positions and the narrowing coefficients
    nu0 = np.diag([(i[0]) for i in linelist["wavenumber"]])
    n_nu0 = len(linelist["wavenumber"])

    # Vector containing the relative populations of the inferior levels
    popu = np.array([linelist["statistical weight"][i] * np.exp(-c2_constant * linelist["energy"][i] / temperature) / tips for i in range(0, n_nu0, 1)])
    rho = np.diag(popu)

    # Vector X containing the square root of the line intensities divided by relative population, wavenumber, and (1-exp(-E/kT))
    X = np.array([math.sqrt((get_intensity(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i)) / (popu[i] * linelist["wavenumber"][i][0] * (1 - np.exp(-c2_constant * linelist["wavenumber"][i][0] / temperature)))) for i in range(0, len(linelist["intensity"]), 1)], dtype=np.double)

    W = generate_W_nosd(mole_fraction, pressure, n_nu0, mass, temperature, linelist, offdiag, popu)

    # Diagonalize matrix
    h = nu0 - (W * 1j)
    eival, eimat = np.linalg.eig(h)
    eimat_ = np.linalg.inv(eimat)

    sigma1 = np.matmul(X.transpose(), eimat)
    sigma2 = np.matmul(eimat_, rho)
    sigma2 = np.matmul(sigma2, X)
    # X being a 1D array, sigma1 and sigma2 are also 1D array

    sigma = np.zeros(len(vnu), dtype = np.double)


    for i in range(0, len(linelist["wavenumber"]), 1):
        intensity, (lorentz_width, lorentz_shift, doppler_width), fo_lm = get_params(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i, get_voigt_par)
        lowest_value_ = lowest_value / intensity # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        if profile == "rautian":
            narrowing = get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i)
        else:
            narrowing = 0.0
        
        lim0, lim1, calc_check = calc_lims(lowest_value_, doppler_width, lorentz_width, lorentz_shift, linelist["wavenumber"][i][0], vnu)
        if (lim1 - lim0) % 2 == 0:
            if len(vnu) > lim1 + 1:
                lim1 += 1
            else:
                lim1 -= 1

        if calc_check:
            nu0mid = (vnu[lim0] + vnu[lim1-1])/2

            rau = rautian.calc_rautian(vnu[lim0:lim1], nu0mid, doppler_width, 0.0, narrowing, 0.0, 0.0)

            sigma_temp = -(((sigma1[i] * sigma2[i])/(vnu[lim0:lim1] - eival[i])).imag)
            sigma_temp *= vnu[lim0:lim1] * (1 - np.exp(-c2_constant * vnu[lim0:lim1] / temperature))
            if profile != "lorentz":
                sigma_temp = np.concatenate(([sigma_temp[0]]*int(len(rau)/2), sigma_temp), axis = None)
                sigma_temp = np.concatenate((sigma_temp, [sigma_temp[-1]]*int(len(rau)/2)), axis = None)
                conv = signal.fftconvolve(sigma_temp, rau, mode = "valid") * ((vnu[-1] - vnu[0]) / (len(vnu) - 1))
                sigma[lim0:lim1] += signal.fftconvolve(sigma_temp, rau, mode = "valid") * ((vnu[-1] - vnu[0]) / (len(vnu) - 1))
            else:
                sigma[lim0:lim1] += sigma_temp

    sigma /= np.pi

    return sigma
#---------------------------------------------------------------------------------------------------

# Calculate speed-dependent absorption cross section using relaxation matrix -----------------------

def calc_sd_mat(profile : str, pressure, temperature, mass, mass_p, mole_fraction, linelist, offdiag, tips, vnu, lowest_value, method, tips_ratio, intensity_factor):

    c2_constant = 1.438776877 # cm K


    # Maxwell Boltzmann distibution: Velocities and distribution (Eq. 1.18 of Lance's PhD thesis)
    n_v, u_max = 21, 4.0                      # u = v/vp where vp is the most probable speed
    boltzmann_constant, avogadro_constant = 1.380649E-16, 6.02214076E23
    vp = math.sqrt(2.0*boltzmann_constant * temperature / (mass/avogadro_constant))
    v = np.linspace(0, u_max*vp, n_v, dtype = np.double)

    # Sampling of "cos(\theta)" (range = [-1,+1])
    n_mu = 21
    mu = np.linspace(-1.0, 1.0, n_mu, dtype = np.double)
    
    # Diagonal matrices containing the line positions and the narrowing coefficients
    nu0 = np.diag([(i[0]) for i in linelist["wavenumber"]])
    n_nu0 = len(linelist["wavenumber"])

    # Vector and matrix containing the relative populations of the inferior levels
    popu = np.array([linelist["statistical weight"][i] * np.exp(-c2_constant * linelist["energy"][i] / temperature) / tips for i in range(0, n_nu0, 1)])
    rho = np.diag(popu)

    # Vector X containing the square root of the line intensities divided by relative population, wavenumber, and (1-exp(-E/kT))
    X = np.array([math.sqrt((get_intensity(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i)) / (popu[i] * linelist["wavenumber"][i][0] * (1 - np.exp(-c2_constant * linelist["wavenumber"][i][0] / temperature)))) for i in range(0, len(linelist["intensity"]), 1)], dtype=np.double)


    # The absorption coefficient is computed according to Eq. 21, i.e. assuming a speed-dependent Rautian profile with line mixing.
    # The elements involved in it are created here below, for each wavenumber
    if method == "general":
        W = generate_W(mole_fraction, pressure, n_nu0, mass, mass_p, temperature, linelist, offdiag, popu, v/vp)
    elif method == "correlation":
        W, C = generate_W_diag(mole_fraction, pressure, n_nu0, mass, mass_p, temperature, linelist, offdiag, popu, v/vp)

    sigma = np.zeros(len(vnu), dtype = np.double)
    lims = [[0.0, 0.0] for i in range(0, len(linelist["wavenumber"]), 1)]
    dnu = ((vnu[-1] - vnu[0]) / (len(vnu) - 1))
    for i in range(0, len(linelist["wavenumber"]), 1):
        intensity, (lorentz_width, lorentz_shift, doppler_width), fo_lm = get_params(pressure, temperature, mass, mole_fraction, tips_ratio, intensity_factor, linelist, i, get_voigt_par)
        lowest_value_ = lowest_value / intensity # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        if lorentz_width < 0.0008:
            #print("Lorentz width too small -> switching to first order")
            Y = get_Y(W[0], X, nu0, pressure)
            for j in range(0, len(linelist["line-mixing_fo"]), 1):
                linelist["line-mixing_fo"][j] = (Y[j], 0.0, False)
            return calc_qsdprofile(profile, pressure, temperature, mass, mole_fraction, linelist, vnu, lowest_value, tips_ratio, intensity_factor)

        lims[i][0], lims[i][1], calc_check = calc_lims(lowest_value_, doppler_width, lorentz_width, lorentz_shift, linelist["wavenumber"][i][0], vnu)

    lim0 = min(np.array(lims)[:,0])
    lim1 = max(np.array(lims)[:,1])

    # Diagonal matrix of velocity-changing collision parameters
    beta = np.diag([0.0 for i in linelist["wavenumber"]])

    if method == "general":
        if profile in ["qsd_rautian", "qsd_nelkin-ghatak"]:
            beta = np.diag([get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i) for i in range(0, len(linelist["narrowing"]), 1)])
            sigma[lim0:lim1] += calc_g.calc_abs_diag(vnu[lim0:lim1], v, mu, vp, nu0, beta, rho, X, W, temperature)
        elif profile == "modified_hartmann-tran":
            beta = np.diag([get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i) + (1j * get_narrowing_i(pressure, temperature, mass, mole_fraction, linelist, i)) for i in range(0, len(linelist["narrowing"]), 1)])
            sigma[lim0:lim1] += calc_g.calc_abs_diag(vnu[lim0:lim1], v, mu, vp, nu0, beta, rho, X, W, temperature)
        elif profile == "hartmann-tran":
            beta = [np.diag([get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i) - (linelist["correlation_HT"] * W[j][i][i]) for i in range(0, len(linelist["narrowing"]), 1)]) * pressure for j in range(0, len(W), 1)]
            sigma[lim0:lim1] += calc_g.calc_abs_diag_ht(vnu[lim0:lim1], v, mu, vp, nu0, beta, rho, X, W, temperature)
        else:
            sigma[lim0:lim1] += calc_g.calc_abs_diag(vnu[lim0:lim1], v, mu, vp, nu0, beta, rho, X, W, temperature)
    if method == "correlation":
        if profile in ["qsd_rautian", "qsd_nelkin-ghatak"]:
            beta = np.diag([get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i) for i in range(0, len(linelist["narrowing"]), 1)])
            sigma[lim0:lim1] += calc_g.calc_abs_corr(vnu[lim0:lim1], v, mu, vp, nu0, beta, rho, X, C, W, temperature)
        elif profile == "modified_hartmann-tran":
            beta = np.diag([get_narrowing(pressure, temperature, mass, mole_fraction, linelist, i) + (1j * get_narrowing_i(pressure, temperature, mass, mole_fraction, linelist, i)) for i in range(0, len(linelist["narrowing"]), 1)])
            sigma[lim0:lim1] += calc_g.calc_abs_corr(vnu[lim0:lim1], v, mu, vp, nu0, beta, rho, X, C, W, temperature)
        else:
            sigma[lim0:lim1] += calc_g.calc_abs_corr(vnu[lim0:lim1], v, mu, vp, nu0, beta, rho, X, C, W, temperature)

    return sigma

#---------------------------------------------------------------------------------------------------

# Calibrate the spectrum ---------------------------------------------------------------------------
def calibrate_spectrum(nu, trans, xcal, baseline = [], extrapolate : bool = False):
    nnu = [i * (1 - xcal) for i in nu]
    interplot = interpolate.interp1d(nnu, trans, kind = "cubic")
    nu0 = nnu[0]
    nu1 = nnu[-1]
    nnu = [i for i in nu if nu0 <= i <= nu1]
    ntrans = interplot(nnu)
    if extrapolate:
        if nu[0] < nnu[0]:
            ntrans = np.concatenate((baseline[0 : int(len(nu) - len(nnu))], ntrans), axis = None)
        elif nu[-1] > nnu[-1]:
            ntrans = np.concatenate((ntrans, baseline[len(nnu):]), axis = None)
        return [nu, np.array(ntrans)]
    return [np.array(nnu), np.array(ntrans)]

#---------------------------------------------------------------------------------------------------

# Calculate the spectrum using ILS * (e^{alpha l}) -------------------------------------------------

def calc_spectrum(spectral_data, spectrum_index, linelists, offdiags, nu, apply_xcal : bool = True):

    spectrum_form = spectral_data["spectra"][spectrum_index]["form"]

    #Resampling

    nu = nu[0::spectral_data["spectra"][spectrum_index]["downsample_factor"]]
    lnu = 0
    for j in range(0, len(nu), 1):
        if j % spectral_data["spectra"][spectrum_index]["downsample_factor"] == 0:
            lnu = j

    if spectral_data["spectra"][spectrum_index]["downsample_factor"] != 1:
        pass    # Will redo with interp1d
        #inter = np.fft.irfft(alpha)
        #inter = np.concatenate((inter[:int(len(inter)/2)], [0]* (((lnu+1) - len(nu)) * 2), inter[int(len(inter)/2):]), axis = None)
        #alpha = np.fft.rfft(inter).real
        #if len(alpha) != len(nu) - 1:
        #    alpha = np.append(alpha, [alpha[-1]] * (len(nu) - len(alpha)), axis = None)


    dnu = ((nu[-1] - nu[0]) / (len(nu) - 1))
    if apply_xcal:
        nnu = [i * (1 - spectral_data["calculation"]["x_calibration_factor"]) for i in nu]
        nu = np.array(nnu)
        del nnu

    alpha = np.zeros(nu.shape[0], dtype = np.double)

    for j in range(0, len(linelists), 1):
        for k in range(0, len(offdiags[j]), 1):
            alpha += calc_alpha(spectral_data["calculation"]["line_profile"], spectral_data["spectra"][spectrum_index], spectral_data["linelists"][j], linelists[j][k], j, nu, offdiags[j][k], spectral_data["calculation"]["method"])

        if linelists[j] != []:
            alpha += calc_alpha(spectral_data["calculation"]["line_profile"], spectral_data["spectra"][spectrum_index], spectral_data["linelists"][j], linelists[j][-1], j, nu)


    if spectrum_form == "absorption_coefficient":
        alpha += calc_baseline(spectral_data["spectra"][spectrum_index]["baseline"], nu)
        return alpha.real

    elif spectrum_form == "absorption_cross-section":
        mole_fraction = 1.0
        cross = np.zeros(nu.shape[0], dtype = np.double)
        for j in range(0, len(linelists), 1):
            for i in spectral_data["spectra"][spectrum_index]["mole_fraction"]:
                if j + 1 in i[1]:
                    mole_fraction = i[0]
            cross += alpha / (constants.physical_constants["Loschmidt constant (273.15 K, 101.325 kPa)"][0] * 1E-6 * (273.15 / spectral_data["spectra"][spectrum_index]["temperature"]) * (mole_fraction * spectral_data["spectra"][spectrum_index]["total_pressure"] / 1013.25))

        return cross.real

    elif spectrum_form == "transmittance":
        transmittance = generate_transmittance(alpha, spectral_data["spectra"][spectrum_index]["path_length"])
        if spectral_data["spectra"][spectrum_index]["ils_type"][0] == "external":
            ils = get_ILS(spectral_data["spectra"][spectrum_index]["ils_type"][1], nu[int(nu.size / 2)], nu)
            transmittance = np.concatenate(([transmittance[0]]*int(len(ils)/2), transmittance), axis = None)
            transmittance = np.concatenate((transmittance, [transmittance[-1]]*int(len(ils)/2)), axis = None)
            transmittance = signal.fftconvolve(transmittance, ils, mode = "valid") * dnu

        transmittance *= calc_baseline(spectral_data["spectra"][spectrum_index]["baseline"], nu)

        return transmittance.real


#-------------------------------------------------------------------------------------------------
# Create W matrix

# Speed-independent relaxation matrix
def generate_W_nosd(molfrac, press, n_nu0, m, temperature, linelist, offdiag, popu):

    # Loop over (N(N+1))/2 to create relaxation matrix elements (row(i) and column(i)
    # then rown(i+1) and column(i+1) minus the elements of the previous rows and columns)
    #       j1  j2
    #   i1  11  12
    #   i2  21  22
    W = np.zeros((n_nu0, n_nu0), dtype = np.cdouble)
    for i in range(n_nu0):
        # self-broadening(v) + i self-shift(v) + foreign-broadening(v) + i foreign-shift(v)
        lorentz_width, lorentz_shift = get_lorentz_par(press, temperature, m, molfrac, linelist, i)
        W[i][i] = lorentz_width + (lorentz_shift * 1j)
        for j in range(i+1, n_nu0):
            # check if upper triangle element (ij) exists in offdiag
            if "{} {}".format(linelist["name"][i], linelist["name"][j]) in offdiag["names"]:
                W[i,j] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][i], linelist["name"][j]))][0] * press * (-1)                    # upper triangle
                W[j,i] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][i], linelist["name"][j]))][0] * press * (-1) * popu[j]/popu[i]  # lower triangle: Wji = Wij rhoj/rhoi
            # check if lower triangle element (ji) exists in offdiag
            elif "{} {}".format(linelist["name"][j], linelist["name"][i]) in offdiag["names"] :
                W[i,j] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][j], linelist["name"][i]))][0] * press * (-1) * popu[i]/popu[j]  # upper triangle
                W[j,i] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][j], linelist["name"][i]))][0] * press * (-1)                    # lower triangle

    return W

# Speed-dependent relaxation matrix
def generate_W(molfrac, press, n_nu0, m, mp, temperature, linelist, offdiag, popu, u, sdkind = "quadratic") :

    # Loop over (N(N+1))/2 to create relaxation matrix elements (row(i) and column(i)
    # then rown(i+1) and column(i+1) minus the elements of the previous rows and columns)
    #       j1  j2
    #   i1  11  12
    #   i2  21  22
    W = np.zeros((n_nu0, n_nu0, len(u)), dtype = np.cdouble)
    for i in range(n_nu0) :
        sd_width, sd_shift = get_qsd_lorentz_par(press, temperature, m, molfrac, linelist, i)
        W[i][i] = sd_width + (sd_shift * 1j)
        for j in range(i+1, n_nu0) :
            # check if upper triangle element (ij) exists in offdiag
            if "{} {}".format(linelist["name"][i], linelist["name"][j]) in offdiag["names"]:
                W[i,j] = get_line_mixing(press, temperature, m, molfrac, offdiag, offdiag["names"].index("{} {}".format(linelist["name"][i], linelist["name"][j]))) * (-1)                    # upper triangle
                W[j,i] = get_line_mixing(press, temperature, m, molfrac, offdiag, offdiag["names"].index("{} {}".format(linelist["name"][i], linelist["name"][j]))) * (-1) * popu[j]/popu[i]  # lower triangle: Wji = Wij rhoj/rhoi
            # check if lower triangle element (ji) exists in offdiag
            elif "{} {}".format(linelist["name"][j], linelist["name"][i]) in offdiag["names"] :
                W[i,j] = get_line_mixing(press, temperature, m, molfrac, offdiag, offdiag["names"].index("{} {}".format(linelist["name"][j], linelist["name"][i]))) * (-1) * popu[i]/popu[j]  # upper triangle
                W[j,i] = get_line_mixing(press, temperature, m, molfrac, offdiag, offdiag["names"].index("{} {}".format(linelist["name"][j], linelist["name"][i]))) * (-1)                    # lower triangle

    return W.transpose(2, 0, 1)


# Diagonal relaxation matrix + correlation matrix
def generate_W_diag(molfrac, press, n_nu0, m, mp, temperature, linelist, offdiag, popu, u, sdkind = "quadratic") :

    # Loop over (N(N+1))/2 to create relaxation matrix elements (row(i) and column(i)
    # then rown(i+1) and column(i+1) minus the elements of the previous rows and columns)
    #       j1  j2
    #   i1  11  12
    #   i2  21  22
    W = np.zeros((n_nu0, n_nu0, len(u)), dtype = np.cdouble)
    C = np.zeros((n_nu0, n_nu0), dtype = np.cdouble)
    for i in range(n_nu0) :
        sd_width, sd_shift = get_qsd_lorentz_par(press, temperature, m, molfrac, linelist, i)
        W[i][i] = sd_width + (sd_shift * 1j)

        for j in range(i+1, n_nu0) :
            # check if upper triangle element (ij) exists in offdiag
            if "{} {}".format(linelist["name"][i], linelist["name"][j]) in offdiag["names"]:
                C[i,j] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][i], linelist["name"][j]))][0] * press * (-1)                    # upper triangle
                C[j,i] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][i], linelist["name"][j]))][0] * press * (-1) * popu[j]/popu[i]  # lower triangle: Cji = Cij rhoj/rhoi
            # check if lower triangle element (ji) exists in offdiag
            elif "{} {}".format(linelist["name"][j], linelist["name"][i]) in offdiag["names"] :
                C[i,j] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][j], linelist["name"][i]))][0] * press * (-1) * popu[i]/popu[j]  # upper triangle
                C[j,i] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][j], linelist["name"][i]))][0] * press * (-1)                    # lower triangle

    return W.transpose(2, 0, 1), C


# Return Y parameters
def get_Y(W, X, nu0, press):
    Y = np.zeros(len(X), dtype = np.double)
    for i in range(0, len(X), 1):
        Y[i] = 2 * np.sum([(X[j]/X[i]) * (W[j][i] / (nu0[i][i] - nu0[j][j])) for j in range(0, len(W[i]), 1) if j != i]).real / press
    return Y



#-------------------------------------------------------------------------------------------------
# Create ILS

def get_ILS(ils_file, nu0, nu):
    domain = ""
    par = dict()
    datablock = False
    with open(ils_file) as f:
        blocktype = ""
        for x, line in enumerate(f):
            if line.startswith("$GENERAL"):
                blocktype = "general"
            elif line.startswith("$PARAMETERS"):
                blocktype = "parameters"
            elif line.startswith("$DATA"):
                datablock = True
                break
            elif blocktype == "general":
                if line.startswith("domain ="):
                    domain = ctrl_pars.get_parameter(x, line)
            elif blocktype == "parameters":
                if line.startswith("ft_exponent_sign"):
                    par["ft_exponent_sign"] = ctrl_pars.get_parameter(x, line)
                elif line.startswith("polynomial_order"):
                    par["polynomial_order"] = int(ctrl_pars.get_parameter(x, line))
                elif line.startswith("iris_diameter"):
                    par["iris_diameter"] = float(ctrl_pars.get_parameter(x, line))
                elif line.startswith("focal_length"):
                    par["focal_length"] = float(ctrl_pars.get_parameter(x, line))
                elif line.startswith("ils_size"):
                    par["ils_size"] = float(ctrl_pars.get_parameter(x, line))
                elif line.startswith("igm_interpolation"):
                    par["igm_interpolation"] = ctrl_pars.get_parameter(x, line)
                elif line.startswith("mopd"):
                    par["mopd"] = float(ctrl_pars.get_parameter(x, line))
                elif line.startswith("apodization"):
                    par["apod"] = ctrl_pars.get_parameter(x, line)
    if domain.lower() == "time" or domain.lower() == "opd" or not datablock:
        return calc_ILS_time(ils_file, nu0, nu, par)
    elif (domain.lower() == "frequency" or domain.lower() == "wavenumber") and datablock:
        return calc_ILS_frequency(ils_file, nu0, nu, par)
    print("ERROR: unrecognized ILS domain\n\nAccepted domains are: \"time\", \"opd\", \"frequency\", or \"wavenumber\"")
    sys.exit()
    return

# ILS from modulation
def calc_ILS_time(ils_file, nu0, nu, par):

    polynu0 = []
    modeffs = []
    phases = []
    data = {"wavenumber": [], "opd": [], "modeff": modeffs, "phase": phases}

    with open(ils_file) as f:
        count = -1
        blocktype = "none"
        for x, line in enumerate(f):
            if line.startswith("$DATA"):
                blocktype = "data"
                count += 1
                data["opd"].append([])
                modeffs.append([])
                phases.append([])
            elif blocktype == "data":
                if line.startswith("wavenumber ="):
                    data["wavenumber"].append(float(ctrl_pars.get_parameter(x, line)))
                elif not line.startswith("%") and len(line.strip()) != 0:
                    data["opd"][count].append(float(line[0:10]))
                    modeffs[count].append(float(line.split()[1]))
                    phases[count].append(float(line.split()[2]))


    if len(modeffs) > 0:
        if len(data["wavenumber"]) <= par["polynomial_order"]:
            par["polynomial_order"] = len(data["wavenumber"]) - 1
        modeff = []
        phase = []
        A = np.ones((len(data["wavenumber"]), par["polynomial_order"] + 1), dtype = np.double)
        B = np.ones((1, par["polynomial_order"] + 1), dtype = np.double)
        for i in range(0, A.shape[0], 1):
            for j in range(0, A.shape[1], 1):
                A[i][j] = data["wavenumber"][i]**j
        for i in range(0, B.shape[1], 1):
            B[0][i] = nu0**i

        opd = data["opd"][0]
        for i in range(0, len(opd), 1):
            yr = np.array([modeffs[j][i] for j in range(0, len(data["wavenumber"]), 1)])
            yi = np.array([phases[j][i] for j in range(0, len(data["wavenumber"]), 1)])
            polyparrs = np.linalg.lstsq(A, yr, rcond = 1E-15)[0]
            polyparis = np.linalg.lstsq(A, yi, rcond = 1E-15)[0]
            modeff.append(np.matmul(B, polyparrs)[0])
            phase.append(np.matmul(B, polyparis)[0])
    else:
        opd = [par["mopd"]]

    nu_inter = (nu[-1] - nu[0])/(len(nu) - 1)   # Interval between each nu[i]
    n_opd = int((2 * int(par["ils_size"] / nu_inter)) + 1)  # Number of points for the ILS
    nu_max = nu_inter * (((n_opd + 1)/2) - 1)   # Nu_max of the ILS

    opd_max = 1/(2*nu_inter)
    opd_inter = 1/(2*nu_max)
    n_mopd = int((opd[-1]/opd_inter) + 1)    # Number of points up to MOPD
    mopd = (n_mopd - 1) * opd_inter # Nearest point to MOPD

    if len(modeffs) > 0:
        interplotr = interpolate.interp1d(opd, modeff, kind = "cubic")
        interploti = interpolate.interp1d(opd, phase, kind = "cubic")
        opd = np.linspace(0, mopd, n_mopd)  # OPD using right interval
        modeff = interplotr(opd)
        phase = interploti(opd)
    else:
        modeff = [1.0] * n_mopd
        phase = [0.0] * n_mopd

    opd = np.linspace(-opd_max, opd_max, n_opd)
    modulation = [0 + 0j] * n_opd
    fils = np.sinc(((par["iris_diameter"]**2) * nu0 * opd) / (8 * (par["focal_length"]**2)))
    if "apod" in par:
        if par["apod"] == "Norton-Beer Weak":
            fils *= norton_beer([0.384093, -0.087577, 0.703484], opd, mopd)
        else:
            print("WARNING: Unrecognized apodization \"{}\"".format(apod))

    modulation[int(n_opd/2)] = 1
    for i in range(1, n_mopd, 1):
        modulation[int(n_opd/2) + i] = modeff[i] * (1 + 1j * np.tan(phase[i]))
        modulation[int(n_opd/2) - i] = modeff[i] * (1 + 1j * np.tan(-phase[i]))
    modulation = np.array(modulation)

    if par["ft_exponent_sign"] == "positive":
        ils = np.fft.fftshift(np.fft.ifft(np.fft.ifftshift(modulation * fils))) * opd_inter
    else:
        ils = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(modulation * fils))) * opd_inter
    nu = np.linspace(-nu_max, nu_max, n_opd)
    # wavelength cos(ax) = 2pi/a
    wavelength = 2/(2*mopd)
    # if cos bell starts at nu[-1] - wavelength: (cos((2pi * nu /wavelength) + (pi * (nu[-1] - wavelength/2)) + 1) / 2
    for i in range(0, len(nu), 1):
        if nu[i] >= nu_max - (wavelength / 2):
            ils[i] *= (np.cos((2 * np.pi * nu[i] / wavelength) -  ((nu_max - (wavelength / 2)) * 2 * np.pi / wavelength)) + 1)/2
        elif nu[i] <= -nu_max + (wavelength / 2):
            ils[i] *= (np.cos((2 * np.pi * nu[i] / wavelength) +  ((nu_max - (wavelength / 2)) * 2 * np.pi / wavelength)) + 1)/2

    area = integrate.simpson(ils.real, dx = nu_inter)
    ils /= area


    return ils


# Pre-calculated ILS
def calc_ILS_frequency(ils_file, nu0, nu, par):

    polynu0 = []
    modeffs = []
    phases = []
    nuilss = []
    ilss = []
    data = {"wavenumber": [], "nu": nuilss, "ils": ilss}

    with open(ils_file) as f:
        count = -1
        blocktype = "none"
        for x, line in enumerate(f):
            if line.startswith("$DATA"):
                blocktype = "data"
                count += 1
                nuilss.append([])
                ilss.append([])
            elif blocktype == "data":
                if line.startswith("wavenumber ="):
                    data["wavenumber"].append(float(ctrl_pars.get_parameter(x, line)))
                elif not line.startswith("%") and len(line.strip()) != 0:
                    nuilss[count].append(float(line.split()[0]))
                    ilss[count].append(float(line.split()[1]))


    nudiff = [abs(nu0 - i) for i in data["wavenumber"]]
    ilsindex = np.argmin(nudiff)

    nu_inter = (nu[-1] - nu[0])/(len(nu) - 1)   # Interval between each nu[i]
    if "ils_size" in par and par["ils_size"] < min(abs(nuilss[ilsindex][0]), abs(nuilss[ilsindex][-1])):
        n_opd = int((2 * int(par["ils_size"] / nu_inter)) + 1)  # Number of points for the ILS
        nu_max = nu_inter * (((n_opd + 1)/2) - 1)   # Nu_max of the ILS
    else:
        par["ils_size"] = min(abs(nuilss[ilsindex][0]), abs(nuilss[ilsindex][-1]))
        n_opd = int((2 * int(par["ils_size"] / nu_inter)) + 1)  # Number of points for the ILS
        nu_max = nu_inter * (((n_opd + 1)/2) - 1)   # Nu_max of the ILS

    nu = np.linspace(-nu_max, nu_max, n_opd)

    interplot = interpolate.interp1d(nuilss[ilsindex], ilss[ilsindex], kind = "cubic")
    ils = interplot(nu)


    return ils

#-------------------------------------------------------------------------------------------------
# Convert absorption coefficient to transmittance

def generate_transmittance(alpha, l):
    tr = np.exp(-alpha * l)
    return tr

#-------------------------------------------------------------------------------------------------
# Calculate baseline

def calc_baseline(baseline_par, nu):

    baseline  = np.zeros(len(nu), dtype = np.double)
    for i in range(0, len(baseline_par), 1):
        baseline += baseline_par[i] * ((nu - ((nu[-1] + nu[0]) / 2)) ** i)

    return baseline


