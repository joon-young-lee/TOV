import numpy as np
import matplotlib.pyplot as plt
import units_cgs as cgs
from scipy.interpolate import CubicSpline

def eos(file_name):
    data = np.loadtxt(file_name)

    e_arr = data[:, 0]
    p_arr = data[:, 1]
    e_arr = np.flip(e_arr * cgs.MeV_fm_to_cgs)
    p_arr = np.flip(p_arr * cgs.MeV_fm_to_cgs)

    eos = CubicSpline(p_arr, e_arr, bc_type='natural', extrapolate=False)
    
    
    return eos, p_arr

def pure_eos(p):
    e_0 = cgs.m_n ** 4 * cgs.c ** 5/(cgs.pi ** 2 * cgs.hbar3)
    P_x = lambda x: e_0/24 * ((2 * x ** 3 - 3 * x)*np.sqrt(x**2 + 1) \
                              + 3 * np.log(x + np.sqrt(x**2 + 1.0)))
    E_x = lambda x: e_0/8 * ((2 * x ** 3 + x)*np.sqrt(x**2 + 1) - \
                             np.log(x + np.sqrt(x**2 + 1.0)))
    h = 1.0e-6
    deriv_P = lambda x: e_0/12 * (4 * x ** 4 - 3 * x ** 2) / np.sqrt(1+x**2)
    
    x_old = 1.0e-1
    x_new = x_old - (P_x(x_old) - p) / deriv_P(x_old)

    j = 1
    maxit = 100
    tol = 1.e-5
    
    while j <= maxit and np.abs(x_new - x_old) > tol:
        
        x_old = x_new
        x_new = x_old - (P_x(x_old) - p) / deriv_P(x_old)
        
        
        j += 1
    #print(x_new)
    return E_x(x_new)

e_0 = cgs.m_n ** 4 * cgs.c ** 5/(cgs.pi ** 2 * cgs.hbar3)
P_x = lambda x: e_0/24 * ((2 * x ** 3 - 3 * x)*np.sqrt(x**2 + 1) \
                              + 3 * np.log(x + np.sqrt(x**2 + 1.0)))



data = np.loadtxt('Archive/eft_pnm32_000078_ldm_eos.dat')
e = data[:, 0] * cgs.MeV_fm_to_cgs
p = data[:, 1] * cgs.MeV_fm_to_cgs
print(np.max(p))
print(np.min(p))

