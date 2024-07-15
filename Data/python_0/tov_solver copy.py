import numpy as np
from scipy.interpolate import CubicSpline
import units_cgs as cgs

# All units in cgs

# Using Runge-Kutta 4th order
def eos(file_name):
    data = np.loadtxt(file_name)

    ## e_arr = data[:, 0]
    ## p_arr = data[:, 1]
    e_arr = cgs.MeV_fm_to_cgs * np.flip(data[:, 0])
    p_arr = cgs.MeV_fm_to_cgs * np.flip(data[:, 1])

    eos = CubicSpline(p_arr, e_arr, bc_type='natural', extrapolate=False)
    
    
    return eos, p_arr


def TOV(EoS, p0, del_r, maxit): # inputs: gamma, max iteration,  ∆r, initial pressure
    
    

    dp_dr = lambda r, p, m: -cgs.G * EoS(p) * m * \
        (1 + p/EoS(p)) * \
        (1 + 4 * cgs.pi * p * r ** 3/(m * cgs.c2)) / \
        (1 - 2 * cgs.G * m / (cgs.c2 * r)) / \
        (cgs.c * r) ** 2

    dm_dr = lambda r, p, m: 4 * cgs.pi * EoS(p) * r ** 2 / cgs.c2

    r = 1.0e-10
    p = p0
    m = 1.0e-10
    # Runge-Kutta 4th order
    i = 0
    cut = int(4000)
    for i in range(0, cut): # final pressure 3.3587e-18 MeV/fm^3

        k1_p = dp_dr(r, p, m)
        k1_m = dm_dr(r, p, m)
        
        r2 = r + del_r/2
        p2 = p + del_r * k1_p / 2
        m2 = m + del_r * k1_m / 2
        k2_p = dp_dr(r2, p2, m2)
        k2_m = dm_dr(r2, p2, m2)
        
        r3 = r2
        p3 = p + del_r * k2_p / 2
        m3 = m + del_r * k2_m / 2
        k3_p = dp_dr(r3, p3, m3)
        k3_m = dm_dr(r3, p3, m3)
        
        r4 = r2 + del_r/2
        p4 = p + del_r * k3_p
        m4 = m + del_r * k3_m
        k4_p = dp_dr(r4, p4, m4)
        k4_m = dm_dr(r4, p4, m4)
        
        p += (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_r
        m += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_r
        r += del_r
        

    i = cut
    while i <= maxit and p >= 5.38125462909e15: # final pressure 3.3587e-18 MeV/fm^3, end of iteration

        k1_p = dp_dr(r, p, m)
        k1_m = dm_dr(r, p, m)
        
        r2 = r + del_r/2
        p2 = p + del_r * k1_p / 2
        m2 = m + del_r * k1_m / 2
        k2_p = dp_dr(r2, p2, m2)
        k2_m = dm_dr(r2, p2, m2)
        
        r3 = r2
        p3 = p + del_r * k2_p / 2
        m3 = m + del_r * k2_m / 2
        k3_p = dp_dr(r3, p3, m3)
        k3_m = dm_dr(r3, p3, m3)
        
        r4 = r2 + del_r/2
        p4 = p + del_r * k3_p
        m4 = m + del_r * k3_m
        k4_p = dp_dr(r4, p4, m4)
        k4_m = dm_dr(r4, p4, m4)
        
        p += (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_r
        m += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_r
        r += del_r
        i += 1
        
        
        # if i % 100 == 0:
        #     print(r)
        #     print(m)
        #     print(p)


     
    #m /= cgs.M0
    r /= 1e5
    return m, r

# 1.e32 ~ 1.e41, divide each interval by 1e3, 2e3 ...

del_r = 100.0
maxit = 1e1 


for j in range(1, 101):
    EoS, p_arr = eos(f'Archive/eft_pnm32_{j:0>6}_ldm_eos.dat')
    step = 10
    M = np.zeros(step)
    R = np.zeros(step)
    P = np.zeros(step)  
    
    p_i = 1.602179e+33 # 1.0 MeV/fm^3
    p_f = 1.602179e+36 # 1000 MeV/fm^3
    print('-------------------------------------------')
    print(f'Archive/eft_pnm32_{j:0>6}_ldm_eos.dat')
    print('-------------------------------------------')
    for l in range(1, step + 1):
        P[l-1] = p_i + (l-1) * (p_f-p_i)/step
        M[l-1], R[l-1] = TOV(EoS, P[l-1], del_r, maxit)
        P[l-1] *= cgs.erg_cm_to_MeV_fm # MeV/fm^3
        print(f'{l/10}%: Initial Pressure: {P[l-1]:.4f} MeV/fm^3 Mass: {M[l - 1]}M_0 {R[l - 1]}km')

    print('-------------------------------------------')
    print(f'Initial Pressure {p_i} ~ {p_f}')
    data = np.column_stack((P, M, R))
    # np.savetxt(f'Result_2/eft_pnm32_{j:0>6}_ldm_result.dat', data, fmt='%.6E', delimiter=' ', 
                # header='', comments='')
    