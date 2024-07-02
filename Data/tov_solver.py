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


def TOV(EoS, p_arr, p0, del_r, maxit): # inputs: gamma, max iteration,  ∆r, initial pressure
    
    

    dp_dr = lambda r, p, m: -cgs.G * EoS(p) * m * \
        (1 + p/EoS(p)) * \
        (1 + 4 * cgs.pi * p * r ** 3/(m * cgs.c2)) / \
        (1 - 2 * cgs.G * m / (cgs.c2 * r)) / \
        (cgs.c * r) ** 2

    dm_dr = lambda r, p, m: 4 * cgs.pi * EoS(p)  * r ** 2 / cgs.c2

    r = 1.0e-10
    p = p0
    m = 1.0e-10
    # Runge-Kutta 4th order
    i = 0
    while i <= maxit and p > np.min(p_arr):
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

        
        
    m /= cgs.M0
    r /= 1e5
    return m, r

# 1.e32 ~ 1.e41, divide each interval by 1e3, 2e3 ...

del_r = 100.0
maxit = 1e6  


for j in range(1, 101):
    EoS, p_arr = eos(f'Archive/eft_pnm32_{j:0>6}_ldm_eos.dat')
    M = np.zeros(1000)
    R = np.zeros(1000)
    P = np.zeros(1000)  
    step = 1000
    p_i = 1.e34
    p_f = np.max(p_arr)
    for l in range(1, step + 1):
        P[l-1] = p_i + (l-1) * (p_f-p_i)/step
        M[l-1], R[l-1] = TOV(EoS, p_arr, P[l-1], del_r, maxit)
        P[l-1] *= cgs.erg_cm_to_MeV_fm
        print('-------------------------------------------')
        print(f'Archive/eft_pnm32_{j:0>6}_ldm_eos.dat: {l/10}%')
        print(M[l-1], 'M_0')
        print(R[l-1], 'km')
    
    print(f'Initial Pressure {p_i} ~ {p_f}')
    data = np.column_stack((P, M, R))
    np.savetxt(f'Result_1/eft_pnm32_{j:0>6}_ldm_result.dat', data, fmt='%.6E', delimiter=' ', 
                header='', comments='')
    