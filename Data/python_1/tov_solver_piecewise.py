import numpy as np
from scipy.interpolate import CubicSpline
import units_cgs as cgs

# All units in cgs

# Using Runge-Kutta 4th order
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

    deriv_P = lambda x: e_0/12 * (4 * x ** 4 - 3 * x ** 2) / np.sqrt(1+x**2)
    
    x_old = 1.0
    x_new = x_old - (P_x(x_old) - p) / deriv_P(x_old)

    j = 1
    maxit = 100
    tol = 1.e-4
    
    while j <= maxit and np.abs(x_new - x_old) > tol:
        
        x_old = x_new
        x_new = x_old - (P_x(x_old) - p) / deriv_P(x_old)
        
        
        j += 1
    #print(x_new)
    return E_x(x_new)

def TOV(EoS, p_arr, p0, del_r, maxit): # inputs: gamma, max iteration,  ∆r, initial pressure
    
    

    dp_dr = lambda r, p, m: -cgs.G * EoS(p) * m * \
        (1 + p/EoS(p)) * \
        (1 + 4 * cgs.pi * p * r ** 3/(m * cgs.c2)) / \
        (1 - 2 * cgs.G * m / (cgs.c2 * r)) / \
        (cgs.c * r) ** 2

    dm_dr = lambda r, p, m: 4 * cgs.pi * EoS(p)  * r ** 2 / cgs.c2
        
    pure_dp_dr = lambda r, p, m: -cgs.G * pure_eos(p) * m * \
        (1 + p/pure_eos(p)) * \
        (1 + 4 * cgs.pi * p * r ** 3/(m * cgs.c2)) / \
        (1 - 2 * cgs.G * m / (cgs.c2 * r)) / \
        (cgs.c * r) ** 2

    pure_dm_dr = lambda r, p, m: 4 * cgs.pi * pure_eos(p)  * r ** 2 / cgs.c2

    r = 1.0e-10
    p = p0
    m = 1.0e-10
    # Runge-Kutta 4th order
    i = 0
    while i <= maxit and p > 0:
        
        if np.min(p_arr) < p < np.max(p_arr):
            k1_p = dp_dr(r, p, m)
            k1_m = dm_dr(r, p, m)
        
            k2_p = dp_dr(r + del_r/2, p + k1_p/2, m + k1_m/2)
            k2_m = dm_dr(r + del_r/2, p + k1_p/2, m + k1_m/2)
        
            k3_p = dp_dr(r + del_r/2, p + k2_p/2, m + k2_m/2)
            k3_m = dm_dr(r + del_r/2, p + k2_p/2, m + k2_m/2)
        
            k4_p = dp_dr(r + del_r, p+ k3_p, m + k3_m)
            k4_m = dm_dr(r + del_r, p + k3_p, m + k3_m)
        
            p = p + (k1_p + 2*k2_p + 2*k3_p + k4_p) / 6 * del_r
            m = m + (k1_m + 2*k2_m + 2*k3_m + k4_m) / 6 * del_r
            r = r + del_r
            i += 1

        else:
            
            k1_p = pure_dp_dr(r, p, m)
            k1_m = pure_dm_dr(r, p, m)
        
            k2_p = pure_dp_dr(r + del_r/2, p + k1_p/2, m + k1_m/2)
            k2_m = pure_dm_dr(r + del_r/2, p + k1_p/2, m + k1_m/2)
        
            k3_p = pure_dp_dr(r + del_r/2, p + k2_p/2, m + k2_m/2)
            k3_m = pure_dm_dr(r + del_r/2, p + k2_p/2, m + k2_m/2)
        
            k4_p = pure_dp_dr(r + del_r, p+ k3_p, m + k3_m)
            k4_m = pure_dm_dr(r + del_r, p + k3_p, m + k3_m)
        
            p = p + (k1_p + 2*k2_p + 2*k3_p + k4_p) / 6 * del_r
            m = m + (k1_m + 2*k2_m + 2*k3_m + k4_m) / 6 * del_r
            r = r + del_r
            i += 1
        
        
    m /= cgs.M0
    r /= 1e5
    return m, r

# 1.e32 ~ 1.e41, divide each interval by 1e3, 2e3 ...

del_r = 100.0
maxit = 1e6  


for j in range(27, 101):
    EoS, p_arr = eos(f'Archive/eft_pnm32_{j:0>6}_ldm_eos.dat')
    M = np.zeros(1000)
    R = np.zeros(1000)
    P = np.zeros(1000)  
    step = 1000
    p_i = 1.e34
    p_f = 1.602179e+36
    for l in range(1, step + 1):
        P[l-1] = p_i + (l-1) * (p_f-p_i)/step
        M[l-1], R[l-1] = TOV(EoS, p_arr, P[l-1], del_r, maxit)
        P[l-1] *= cgs.erg_cm_to_MeV_fm
        print(f'Archive/eft_pnm32_{j:0>6}_ldm_eos.dat: {l/10}%')
    
    print(f'Initial Pressure {p_i} ~ {p_f}')
    data = np.column_stack((P, M, R))
    np.savetxt(f'Result/eft_pnm32_{j:0>6}_ldm_result.dat', data, fmt='%.4E', delimiter=' ', 
                header='', comments='')
    