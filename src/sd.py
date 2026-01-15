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

#-------------------------------------------------------------------------------------------------

def hypergeo(a, b, z):
    Mv = np.zeros(z.shape[0], dtype = np.double)
    for h in range(0, z.shape[0], 1):
        M = 0.0
        n = 1 + int(5*max(np.abs(a), np.abs(b), np.abs(z[h])))
        for i in range(0, n+1, 1):
            temp = 1.0
            for j in range(0, i, 1):
                temp = temp * ((a + j)/(b + j))
            M = M + (temp * ((z[h]**i)/np.math.factorial(i)))
        Mv[h] = M
    return Mv

#-------------------------------------------------------------------------------------------------

def cumulative_hypergeo(n, a, b, z):
    arr = [0] * (n+1)
    M = 0
    for i in range(0, n+1, 1):
        temp = 1.0
        for j in range(0, i, 1):
            temp = temp * ((a + j)/(b + j))
        M = M + (temp * ((z**i)/np.math.factorial(i)))
        arr[i] = M
    return arr

#-------------------------------------------------------------------------------------------------

def sdpar(p0, u, m, mp, kind = "quadratic", q = 6, dpp0 = 0.1):
    if kind == "quadratic":
        par = p0 * (1 + (dpp0 * ((u**2) - 1.5)))
        return par
    elif kind == "hypergeometric":
        alpha = (q - 3)/(q - 1)
        par = p0 * ((1 + mp/m)**(-alpha/2)) * hypergeo(-alpha/2, 1.5, -(mp/m)*(u**2))
        return par
    print("unrecognized SD kind")
    return 0.0

#-------------------------------------------------------------------------------------------------

#n = 40
#m = 17
#mp = 40
#temp = 296
#x = np.linspace(0, n, n+1)
#vp = np.sqrt(2 * constants.Boltzmann * temp / (m / (constants.Avogadro * 1000)))
#v = np.linspace(0, 3*vp, 21)
##y = cumulative_hypergeo(n, -0.3, 1.5, -(mp/m)*((2*vp/vp)**2))
#q = 4.8
#a = -(q-3)/((2*q)-2)
#b = 3/2
##yv = ((1 + mp/m)**a) * hypergeo(a, b, -(mp/m)*((v/vp)**2))
##yvq = 1 + (0.123 * (((v/vp)**2) - b))
#
#yv = sdpar(0.07, v, vp, m, mp, "hypergeometric", 4.8)
#yvq = sdpar(0.07, v, vp, m, mp, "quadratic", dpp0 = 0.123)
#plt.figure(1)
#plt.plot(v, yv)
#plt.plot(v, yvq)
#plt.show()

#-------------------------------------------------------------------------------------------------
