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

def sdpar(p0, n, T, u, m, mp, kind = "quadratic", q = 6, dpp0 = 0.1):
    if kind == "quadratic":
        par = (p0 * ((296.0 / T)**n)) * (1 + (dpp0 * ((u**2) - 1.5)))
        return par
    elif kind == "hypergeometric":
        alpha = (q - 3)/(q - 1)
        par = (p0 * ((296.0 / T)**n)) * ((1 + mp/m)**(-alpha/2)) * hypergeo(-alpha/2, 1.5, -(mp/m)*(u**2))
        return par
    print("unrecognized SD kind")
    return 0.0


#-------------------------------------------------------------------------------------------------
