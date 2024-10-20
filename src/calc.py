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
import scipy.constants as constants
import scipy.interpolate as interpolate
import scipy.integrate as integrate
import scipy.signal as signal
import sd
import math
import calc_g
import qSDV
import qSDHC
import ctrl_pars
import doppler
import rautian
import tips

#-------------------------------------------------------------------------------------------------
# Calculate limits up to which the line is calculated in the wavenumber axis

def calc_lims(lowest_value, doppler_width, lorentz_width, shift, nu0, vnu):

    dnu = (vnu[-1] - vnu[0]) / (len(vnu) - 1)

    numaxd = dnu    # Ensure that at least 3 points are calculated
    numaxl = dnu
    if lowest_value < (math.sqrt(np.log(2) / math.pi) / doppler_width):
        numaxd = doppler_width * math.sqrt( - np.log(lowest_value / (math.sqrt(np.log(2) / math.pi) / doppler_width)) / math.log(2))
    if lowest_value < (1 / (math.pi * lorentz_width)):
        numaxl = math.sqrt((lorentz_width / (math.pi * lowest_value)) - (lorentz_width ** 2))
    numax = max(numaxd, numaxl)

    return max(int((nu0 + shift - numax - vnu[0]) / dnu), 0), min(int(math.ceil((nu0 + shift + numax - vnu[0]) / dnu)) + 1, len(vnu))

#-------------------------------------------------------------------------------------------------
# Calculate absorption cross section

def calc_alpha(profile, spectrum_data, linelist_data, linelist, linelist_index, tips, vnu, xcal, offdiag = dict(), method = ""):
    mass = linelist_data["mass"]
    mass_p = linelist_data["mass_perturbing"]
    pressure = spectrum_data["total_pressure"] / 1013.25

    mole_fraction = 1
    for i in spectrum_data["mole_fraction"]:
        if linelist_index + 1 in i[1]:
            mole_fraction = i[0]
    lowest_value =  - np.log(1 - spectrum_data["lowest_value"]) / ((constants.physical_constants["Loschmidt constant (273.15 K, 100 kPa)"][0] * 1E-6 * (273.15 / spectrum_data["temperature"]) * (mole_fraction * pressure)) * spectrum_data["path_length"])   # - ln(1 - Transmittance) divided by density * path length

    alpha = np.zeros(len(vnu), dtype = np.double)
    if profile == "qsd_voigt":
        # There is not enough sampled points for v and cos(theta) at low pressure so we neglect line mixing for those
        if offdiag == dict():
            alpha = calc_qsdvoigt(pressure, spectrum_data["temperature"], mass, mole_fraction, xcal, linelist, vnu, lowest_value)
    elif profile == "qsd_rautian":
        # There is not enough sampled points for v and cos(theta) at low pressure so we neglect line mixing for those
        if offdiag == dict() or spectrum_data["total_pressure"] < 10:
            alpha = calc_qsdrautian(pressure, spectrum_data["temperature"], mass, mole_fraction, xcal, linelist, vnu, lowest_value)
        else:
            alpha = calc_sd_mat(pressure, spectrum_data["temperature"], mass, mass_p, mole_fraction, xcal, linelist, offdiag, tips, vnu, lowest_value, method)
    elif profile == "rautian":
        if offdiag == dict():
            alpha = calc_rautian(pressure, spectrum_data["temperature"], mass, mole_fraction, xcal, linelist, vnu, lowest_value)
        else:
            alpha = calc_nosd_mat(pressure, spectrum_data["temperature"], mass, mole_fraction, xcal, linelist, offdiag, tips, vnu, lowest_value, profile)
    elif profile == "voigt":
        if offdiag == dict():
            alpha = calc_voigt(pressure, spectrum_data["temperature"], mass, mole_fraction, xcal, linelist, vnu, lowest_value)
        else:
            alpha = calc_nosd_mat(pressure, spectrum_data["temperature"], mass, mole_fraction, xcal, linelist, offdiag, tips, vnu, lowest_value, profile)

    alpha *= constants.physical_constants["Loschmidt constant (273.15 K, 100 kPa)"][0] * 1E-6 * (273.15 / spectrum_data["temperature"]) * (mole_fraction * pressure)

    return alpha

#-------------------------------------------------------------------------------------------------
# Calculate Voigt absorption cross section using Weideman algorithm from Py4CAtS

def calc_voigt(pressure, temperature, mass, mole_fraction, xcal, linelist, vnu, lowest_value):

    sigma = np.zeros(len(vnu), dtype = np.double)

    for i in range(0, len(linelist["wavenumber"]), 1):

        lowest_value_ = lowest_value / linelist["intensity"][i][0] # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        doppler_width = doppler.hwhm(linelist["wavenumber"][i][0], temperature, mass)
        lorentz_width = ((linelist["self-broadening"][i][0] * mole_fraction) + (linelist["foreign-broadening"][i][0] * (1 - mole_fraction))) * pressure
        shift = ((linelist["self-shift"][i][0] * mole_fraction) + (linelist["foreign-shift"][i][0] * (1 - mole_fraction))) * pressure
        
        lim0, lim1 = calc_lims(lowest_value_, doppler_width, lorentz_width, shift, linelist["wavenumber"][i][0], vnu)

        voigt = rautian.calc_rautian(vnu[lim0:lim1], linelist["wavenumber"][i][0] * (1 + xcal), doppler_width, lorentz_width, 0.0, shift, linelist["line-mixing"][i][0] * pressure) * linelist["intensity"][i][0]

        sigma[lim0:lim1] += voigt

    return sigma

#-------------------------------------------------------------------------------------------------
# Calculate Rautian absorption cross section using Weideman algorithm from Py4CAtS

def calc_rautian(pressure, temperature, mass, mole_fraction, xcal, linelist, vnu, lowest_value):

    sigma = np.zeros(len(vnu), dtype = np.double)

    for i in range(0, len(linelist["wavenumber"]), 1):

        lowest_value_ = lowest_value / linelist["intensity"][i][0] # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        doppler_width = doppler.hwhm(linelist["wavenumber"][i][0], temperature, mass)
        lorentz_width = ((linelist["self-broadening"][i][0] * mole_fraction) + (linelist["foreign-broadening"][i][0] * (1 - mole_fraction))) * pressure
        narrowing = linelist["narrowing"][i][0] * pressure
        shift = ((linelist["self-shift"][i][0] * mole_fraction) + (linelist["foreign-shift"][i][0] * (1 - mole_fraction))) * pressure
        
        lim0, lim1 = calc_lims(lowest_value_, doppler_width, lorentz_width, shift, linelist["wavenumber"][i][0], vnu)

        rau = rautian.calc_rautian(vnu[lim0:lim1], linelist["wavenumber"][i][0] * (1 + xcal), doppler_width, lorentz_width, narrowing, shift, linelist["line-mixing"][i][0] * pressure) * linelist["intensity"][i][0]

        sigma[lim0:lim1] += rau

    return sigma

#-------------------------------------------------------------------------------------------------
# Calculate qSDV absorption cross section using H. Tran et al's algorithm, "Efficient computation of some speed-dependent isolated line profiles", JQSRT 129 (2013) 199-203 and the corresponding erratum, JQSRT 134

def calc_qsdvoigt(pressure, temperature, mass, mole_fraction, xcal, linelist, vnu, lowest_value):

    sigma = np.zeros(len(vnu), dtype = np.double)

    for i in range(0, len(linelist["wavenumber"]), 1):
        sigma_temp = [0.0] * len(vnu)
        lowest_value_ = lowest_value / linelist["intensity"][i][0] # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        doppler_width = doppler.hwhm(linelist["wavenumber"][i][0], temperature, mass)
        lorentz_width = ((linelist["self-broadening"][i][0] * mole_fraction) + (linelist["foreign-broadening"][i][0] * (1 - mole_fraction))) * pressure
        sd_width = linelist["SD-broadening"][i][0] * lorentz_width
        shift = ((linelist["self-shift"][i][0] * mole_fraction) + (linelist["foreign-shift"][i][0] * (1 - mole_fraction))) * pressure
        sd_shift = linelist["SD-shift"][i][0] * shift
        LS_qSDV_R = 0.0
        LS_qSDV_I = 0.0

        lim0, lim1 = calc_lims(lowest_value_, doppler_width, lorentz_width, shift, linelist["wavenumber"][i][0], vnu)

        for j in range(lim0, lim1, 1):
            line_r, line_i = (qSDV.qsdv(linelist["wavenumber"][i][0] * (1 + xcal), doppler_width, lorentz_width, sd_width, shift, sd_shift, vnu[j], LS_qSDV_R, LS_qSDV_I))
            sigma_temp[j] = line_r - (pressure * linelist["line-mixing"][i][0] * line_i)
        sigma_temp = np.array(sigma_temp) * linelist["intensity"][i][0]
        sigma += sigma_temp

    return sigma

#-------------------------------------------------------------------------------------------------
# Calculate qSDR absorption cross section using H. Tran et al's algorithm, "Efficient computation of some speed-dependent isolated line profiles", JQSRT 129 (2013) 199-203 and the corresponding erratum, JQSRT 134

def calc_qsdrautian(pressure, temperature, mass, mole_fraction, xcal, linelist, vnu, lowest_value):

    sigma = np.zeros(len(vnu), dtype = np.double)


    for i in range(0, len(linelist["wavenumber"]), 1):
        sigma_temp = [0.0] * len(vnu)
        lowest_value_ = lowest_value / linelist["intensity"][i][0] # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        doppler_width = doppler.hwhm(linelist["wavenumber"][i][0], temperature, mass)
        lorentz_width = ((linelist["self-broadening"][i][0] * mole_fraction) + (linelist["foreign-broadening"][i][0] * (1 - mole_fraction))) * pressure
        narrowing = linelist["narrowing"][i][0] * pressure
        sd_width = linelist["SD-broadening"][i][0] * lorentz_width
        shift = ((linelist["self-shift"][i][0] * mole_fraction) + (linelist["foreign-shift"][i][0] * (1 - mole_fraction))) * pressure
        sd_shift = linelist["SD-shift"][i][0] * shift
        LS_qSDV_R = 0.0
        LS_qSDV_I = 0.0

        lim0, lim1 = calc_lims(lowest_value_, doppler_width, lorentz_width, shift, linelist["wavenumber"][i][0], vnu)

        for j in range(lim0, lim1, 1):
            line_r, line_i = (qSDHC.qsdhc(linelist["wavenumber"][i][0] * (1 + xcal), doppler_width, lorentz_width, sd_width, shift, sd_shift, narrowing, vnu[j], LS_qSDV_R, LS_qSDV_I))
            sigma_temp[j] = line_r - (pressure * linelist["line-mixing"][i][0] * line_i)
        sigma_temp = np.array(sigma_temp) * linelist["intensity"][i][0]
        sigma += sigma_temp

    return sigma

#-------------------------------------------------------------------------------------------------
# Calculate hard collision absorption cross section using relaxation matrix

def calc_nosd_mat(pressure, temperature, mass, mole_fraction, xcal, linelist, offdiag, tips, vnu, lowest_value, profile):

    lorentz_width = ((linelist["self-broadening"][0][0] * mole_fraction) + (linelist["foreign-broadening"][0][0] * (1 - mole_fraction))) * pressure
    dnu = ((vnu[-1] - vnu[0]) / (len(vnu) - 1))
    interp = 1
    #if dnu > lorentz_width / 10 : # if dnu is bigger than 0.1 lorentzian width
    #    interp = (dnu / (lorentz_width/10))
    #    print(interp)
    nvnu = np.linspace(vnu[0], vnu[-1], int(len(vnu) / interp))

    c2_constant = 1.438776877 # cm K

    # Diagonal matrices containing the line positions and the narrowing coefficients
    nu0 = np.diag([(i[0] * (1 + xcal)) for i in linelist["wavenumber"]])
    n_nu0 = len(linelist["wavenumber"])
    beta = np.diag([i[0] for i in linelist["narrowing"]]) * pressure

    # Vector containing the relative populations of the inferior levels
    popu = np.array([linelist["statistical weight"][i] * np.exp(-c2_constant * linelist["energy"][i] / temperature) / tips for i in range(0, n_nu0, 1)])
    rho = np.diag(popu)

    # Vector X containing the square root of the line intensities divided by relative population, wavenumber, and (1-exp(-E/kT))
    X = np.array([math.sqrt(linelist["intensity"][i][0] / (popu[i] * linelist["wavenumber"][i][0] * (1 - np.exp(-c2_constant * linelist["wavenumber"][i][0])))) for i in range(0, len(linelist["intensity"]), 1)], dtype=np.double)

    W = generate_W_nosd(mole_fraction, pressure, n_nu0, mass, linelist, offdiag, popu)

    # Diagonalize matrix
    h = nu0 - (W * 1j)
    eival, eimat = np.linalg.eig(h)
    eimat_ = np.linalg.inv(eimat)

    sigma1 = np.matmul(X.transpose(), eimat)
    sigma2 = np.matmul(eimat_, rho)
    sigma2 = np.matmul(sigma2, X)
    # X being a 1D array, sigma1 and sigma2 are also 1D array

    sigma = np.zeros(len(nvnu), dtype = np.double)


    for i in range(0, len(linelist["wavenumber"]), 1):
        lowest_value_ = lowest_value / linelist["intensity"][i][0] # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        doppler_width = doppler.hwhm(linelist["wavenumber"][i][0], temperature, mass)
        lorentz_width = ((linelist["self-broadening"][i][0] * mole_fraction) + (linelist["foreign-broadening"][i][0] * (1 - mole_fraction))) * pressure
        shift = ((linelist["self-shift"][i][0] * mole_fraction) + (linelist["foreign-shift"][i][0] * (1 - mole_fraction))) * pressure
        if profile == "rautian":
            narrowing = linelist["narrowing"][i][0] * pressure
        else:
            narrowing = 0.0
        
        lim0, lim1 = calc_lims(lowest_value_, doppler_width, lorentz_width, shift, linelist["wavenumber"][i][0], nvnu)
        nu0mid = (nvnu[lim0] + nvnu[lim1-1])/2

        rau = rautian.calc_rautian(nvnu[lim0:lim1], nu0mid, doppler_width, 0.0, narrowing, 0.0)
        if len(rau) % 2 == 0:
            rau = np.delete(rau, len(rau)-1)

        sigma_temp = -(((sigma1[i] * sigma2[i])/(nvnu[lim0:lim1] - eival[i])).imag)
        sigma_temp *= nvnu[lim0:lim1] * (1 - np.exp(-c2_constant * nvnu[lim0:lim1])) 
        if profile != "lorentz":
            sigma_temp = np.concatenate(([0.0]*int(len(rau)/2), sigma_temp), axis = None)
            sigma_temp = np.concatenate((sigma_temp, [0.0]*int(len(rau)/2)), axis = None)
            sigma[lim0:lim1] += signal.fftconvolve(sigma_temp, rau, mode = "valid") * ((nvnu[-1] - nvnu[0]) / (len(nvnu) - 1))

    interplot = interpolate.interp1d(nvnu, sigma, kind = "cubic")
    sigma = interplot(vnu)
    sigma /= np.pi

    return sigma

#-------------------------------------------------------------------------------------------------
# Calculate uncorrelated speed-dependent (+ narrowing) absorption cross section using relaxation matrix

def calc_sd_mat(pressure, temperature, mass, mass_p, mole_fraction, xcal, linelist, offdiag, tips, vnu, lowest_value, method):

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
    nu0 = np.diag([(i[0] * (1 + xcal)) for i in linelist["wavenumber"]])
    n_nu0 = len(linelist["wavenumber"])
    beta = np.diag([i[0] for i in linelist["narrowing"]]) * pressure

    # Vector and matrix containing the relative populations of the inferior levels
    popu = np.array([linelist["statistical weight"][i] * np.exp(-c2_constant * linelist["energy"][i] / temperature) / tips for i in range(0, n_nu0, 1)])
    rho = np.diag(popu)

    # Vector X containing the square root of the line intensities divided by relative population, wavenumber, and (1-exp(-E/kT))
    X = np.array([math.sqrt(linelist["intensity"][i][0] / (popu[i] * linelist["wavenumber"][i][0] * (1 - np.exp(-c2_constant * linelist["wavenumber"][i][0])))) for i in range(0, len(linelist["intensity"]), 1)], dtype=np.double)


    # The absorption coefficient is computed according to Eq. 21, i.e. assuming a speed-dependent Rautian profile with line mixing.
    # The elements involved in it are created here below, for each wavenumber
    if method == "general":
        W = generate_W(mole_fraction, pressure, n_nu0, mass, mass_p, linelist, offdiag, popu, v/vp)
    elif method == "correlation":
        W, C = generate_W_diag(mole_fraction, pressure, n_nu0, mass, mass_p, linelist, offdiag, popu, v/vp)

    sigma = np.zeros(len(vnu), dtype = np.double)
    lims = [[0.0, 0.0] for i in range(0, len(linelist["wavenumber"]), 1)]
    for i in range(0, len(linelist["wavenumber"]), 1):
        lowest_value_ = lowest_value / linelist["intensity"][i][0] # Absorption cross section divided by the intensity to get the lowest value needed for the profile

        doppler_width = doppler.hwhm(linelist["wavenumber"][i][0], temperature, mass)
        lorentz_width = ((linelist["self-broadening"][i][0] * mole_fraction) + (linelist["foreign-broadening"][i][0] * (1 - mole_fraction))) * pressure
        shift = ((linelist["self-shift"][i][0] * mole_fraction) + (linelist["foreign-shift"][i][0] * (1 - mole_fraction))) * pressure

        lims[i][0], lims[i][1] = calc_lims(lowest_value_, doppler_width, lorentz_width, shift, linelist["wavenumber"][i][0], vnu)

    lim0 = min(np.array(lims)[:,0])
    lim1 = max(np.array(lims)[:,1])
    if method == "general":
        sigma[lim0:lim1] += calc_g.calc_abs_diag(vnu[lim0:lim1], v, mu, vp, nu0, beta, rho, X, W)
    if method == "correlation":
        sigma[lim0:lim1] += calc_g.calc_abs_corr(vnu[lim0:lim1], v, mu, vp, nu0, beta, rho, X, C, W)

    return sigma


#-------------------------------------------------------------------------------------------------
# Calculate the spectrum using ILS * (e^{alpha l})

def calc_spectrum(spectral_data, spectrum_index, linelists, offdiags, nu):

    #Resampling
    nu = nu[0::spectral_data["spectra"][spectrum_index]["downsample_factor"]]
    lnu = 0
    for j in range(0, len(nu), 1):
        if j % spectral_data["spectra"][spectrum_index]["downsample_factor"] == 0:
            lnu = j
    alpha = np.zeros(nu.shape[0], dtype = np.double)

    for j in range(0, len(linelists), 1):
        tips_ = tips.get_tips(spectral_data["linelists"][j]["tips"], spectral_data["spectra"][spectrum_index]["temperature"])
        for k in range(0, len(offdiags[j]), 1):
            alpha += calc_alpha(spectral_data["calculation"]["line_profile"], spectral_data["spectra"][spectrum_index], spectral_data["linelists"][j], linelists[j][k], j, tips_, nu, spectral_data["calculation"]["x_calibration_factor"], offdiags[j][k], spectral_data["calculation"]["method"])


        if linelists[j] != []:
            alpha += calc_alpha(spectral_data["calculation"]["line_profile"], spectral_data["spectra"][spectrum_index], spectral_data["linelists"][j], linelists[j][-1], j, tips_, nu, spectral_data["calculation"]["x_calibration_factor"])


    if spectral_data["spectra"][spectrum_index]["downsample_factor"] != 1:
        inter = np.fft.irfft(alpha)
        inter = np.concatenate((inter[:int(len(inter)/2)], [0]* (((lnu+1) - len(nu)) * 2), inter[int(len(inter)/2):]), axis = None)
        alpha = np.fft.rfft(inter).real
        if len(alpha) != len(nu) - 1:
            alpha = np.append(alpha, [alpha[-1]] * (len(nu) - len(alpha)), axis = None)


    transmittance = generate_transmittance(alpha, spectral_data["spectra"][spectrum_index]["path_length"])
    if spectral_data["spectra"][spectrum_index]["ils_type"][0] == "external":
        ils = get_ILS(spectral_data["spectra"][spectrum_index]["ils_type"][1], nu[int(nu.size / 2)], nu)
        transmittance = np.concatenate(([transmittance[0]]*int(len(ils)/2), transmittance), axis = None)
        transmittance = np.concatenate((transmittance, [transmittance[-1]]*int(len(ils)/2)), axis = None)
        transmittance = signal.fftconvolve(transmittance, ils, mode = "valid") * ((nu[-1] - nu[0]) / (len(nu) - 1))
    transmittance *= calc_baseline(spectral_data["spectra"][spectrum_index]["baseline"], nu)

    return transmittance.real


#-------------------------------------------------------------------------------------------------
# Create W matrix

# Speed-independent relaxation matrix
def generate_W_nosd(molfrac, press, n_nu0, m, linelist, offdiag, popu):

    # Loop over (N(N+1))/2 to create relaxation matrix elements (row(i) and column(i)
    # then rown(i+1) and column(i+1) minus the elements of the previous rows and columns)
    #       j1  j2
    #   i1  11  12
    #   i2  21  22
    W = np.zeros((n_nu0, n_nu0), dtype = np.cdouble)
    for i in range(n_nu0):
        # self-broadening(v) + i self-shift(v) + foreign-broadening(v) + i foreign-shift(v)
        b_self = linelist["self-broadening"][i][0]
        d_self = linelist["self-shift"][i][0]
        b_other = linelist["foreign-broadening"][i][0]
        d_other = linelist["foreign-shift"][i][0]
        W[i][i] = ((b_self + (d_self * 1j))*molfrac + (b_other + (d_other * 1j))*(1 - molfrac)) * press
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
def generate_W(molfrac, press, n_nu0, m, mp, linelist, offdiag, popu, u, sdkind = "quadratic") :

    # Loop over (N(N+1))/2 to create relaxation matrix elements (row(i) and column(i)
    # then rown(i+1) and column(i+1) minus the elements of the previous rows and columns)
    #       j1  j2
    #   i1  11  12
    #   i2  21  22
    W = np.zeros((n_nu0, n_nu0, len(u)), dtype = np.cdouble)
    for i in range(n_nu0) :
        # self-broadening(v) + i self-shift(v) + foreign-broadening(v) + i foreign-shift(v)
        b_self = sd.sdpar(linelist["self-broadening"][i][0], u, m, m, sdkind, dpp0 = linelist["SD-broadening"][i][0])
        d_self = sd.sdpar(linelist["self-shift"][i][0], u, m, m, sdkind, dpp0 = linelist["SD-shift"][i][0])
        b_other = sd.sdpar(linelist["foreign-broadening"][i][0], u, m, mp, sdkind, dpp0 = linelist["SD-broadening"][i][0])
        d_other = sd.sdpar(linelist["foreign-shift"][i][0], u, m, mp, sdkind, dpp0 = linelist["SD-shift"][i][0])
        W[i][i] = ((b_self + d_self * 1j)*molfrac + (b_other + d_other * 1j)*(1 - molfrac)) * press
        for j in range(i+1, n_nu0) :
            # check if upper triangle element (ij) exists in offdiag
            if "{} {}".format(linelist["name"][i], linelist["name"][j]) in offdiag["names"]:
                W[i,j] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][i], linelist["name"][j]))][0] * press * (-1)                    # upper triangle
                W[j,i] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][i], linelist["name"][j]))][0] * press * (-1) * popu[j]/popu[i]  # lower triangle: Wji = Wij rhoj/rhoi
            # check if lower triangle element (ji) exists in offdiag
            elif "{} {}".format(linelist["name"][j], linelist["name"][i]) in offdiag["names"] :
                W[i,j] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][j], linelist["name"][i]))][0] * press * (-1) * popu[i]/popu[j]  # upper triangle
                W[j,i] = offdiag["line-mixing"][offdiag["names"].index("{} {}".format(linelist["name"][j], linelist["name"][i]))][0] * press * (-1)                    # lower triangle

    return W.transpose(2, 0, 1)


# Diagonal relaxation matrix + correlation matrix
def generate_W_diag(molfrac, press, n_nu0, m, mp, linelist, offdiag, popu, u, sdkind = "quadratic") :

    # Loop over (N(N+1))/2 to create relaxation matrix elements (row(i) and column(i)
    # then rown(i+1) and column(i+1) minus the elements of the previous rows and columns)
    #       j1  j2
    #   i1  11  12
    #   i2  21  22
    W = np.zeros((n_nu0, n_nu0, len(u)), dtype = np.cdouble)
    C = np.zeros((n_nu0, n_nu0), dtype = np.cdouble)
    for i in range(n_nu0) :
        # self-broadening(v) + i self-shift(v) + foreign-broadening(v) + i foreign-shift(v)
        b_self = sd.sdpar(linelist["self-broadening"][i][0], u, m, m, sdkind, dpp0 = linelist["SD-broadening"][i][0])
        d_self = sd.sdpar(linelist["self-shift"][i][0], u, m, m, sdkind, dpp0 = linelist["SD-shift"][i][0])
        b_other = sd.sdpar(linelist["foreign-broadening"][i][0], u, m, mp, sdkind, dpp0 = linelist["SD-broadening"][i][0])
        d_other = sd.sdpar(linelist["foreign-shift"][i][0], u, m, mp, sdkind, dpp0 = linelist["SD-shift"][i][0])
        W[i][i] = ((b_self + d_self * 1j)*molfrac + (b_other + d_other * 1j)*(1 - molfrac)) * press

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

    modulation[int(n_opd/2)] = 1
    for i in range(1, n_mopd, 1):
        modulation[int(n_opd/2) + i] = modeff[i] * (1 + 1j * np.tan(phase[i]))
        modulation[int(n_opd/2) - i] = modeff[i] * (1 + 1j * np.tan(-phase[i]))
    modulation = np.array(modulation)


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

    area = 0
    for i in range(0, len(nu)-1, 1):
        area += (ils[i] + ils[i+1]).real * nu_inter / 2
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


