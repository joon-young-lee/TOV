import numpy as np
from scipy.interpolate import CubicSpline
import units_cgs as cgs

# All units in km ( c = G = 1 )

# Using Runge-Kutta 4th order
def eos(file_name):
    data = np.loadtxt(file_name)

    e_arr = np.flip(data[:, 0] * cgs.MeV_fm_to_km)
    p_arr = np.flip(data[:, 1] * cgs.MeV_fm_to_km)

    eos = CubicSpline(p_arr, e_arr, bc_type='natural', extrapolate=False)
    
    
    return eos


def TOV(EoS, p0, del_h, maxit): # inputs: gamma, max iteration,  ∆r, initial pressure
    
    dp_dh = lambda p: -(EoS(p) + p) 
    dr2_dh = lambda r, p, m: -2 * r ** 2 * (r - 2 * m)/(m + 4 * np.pi * p * r ** 3)
    dm_dh = lambda r, p, m: -4 * np.pi * EoS(p) * r ** 3 * (r - 2 * m)/(m + 4 * np.pi * p * r ** 3)
    

    p = p0
    r2 = -3/(2 * np.pi * (EoS(p0) + 3 * p0)) * del_h
    m = 0.0
    # Runge-Kutta 4th order
    
    i = 1
    while i < maxit and p > 1.e-11: # final pressure 3.3587e-18 MeV/fm^3

        dr2 = del_h * dr2_dh(r, p, m)
        

    
    return m, r2

# 1.e32 ~ 1.e41, divide each interval by 1e3, 2e3 ...
file_name = '../Archive/eft_pnm32_000001_ldm_eos.dat'
data = np.loadtxt(file_name)

e_arr = np.flip(data[:, 0] * cgs.MeV_fm_to_km)
p_arr = np.flip(data[:, 1] * cgs.MeV_fm_to_km)
# print(e_arr, 'E_D')
# print(p_arr, 'pressure')



del_h = -1.e-3
maxit = int(1e3)  

p0 = 0.001720545994953356 # 1.3E+3MeV/fm^3
EoS = eos(file_name)
print(TOV(EoS, p0, del_h, maxit))

print(max(p_arr), min(p_arr))



for j in range(1, 2):
    EoS = eos(f'../Archive/eft_pnm32_{j:0>6}_ldm_eos.dat')
    step = 1000
    M = np.zeros(step)
    R = np.zeros(step)
    P = np.zeros(step)  
    
    p_i = 1.720545994953356e-6 # 1.3 MeV/fm^3
    p_f = 1.720545994953356e-3 # 1.3E+3MeV/fm^3
    for l in range(1, step + 1):
        P[l-1] = p_i + (l-1) * (p_f-p_i)/step
        M[l-1], R[l-1] = TOV(EoS, P[l-1], del_h, maxit)
        print('-------------------------------------------')
        print(f'Archive/eft_pnm32_{j:0>6}_ldm_eos.dat: {l/10}%')
        print(f'Initial pressure: {P[l-1] / cgs.MeV_fm_to_km} MeV/fm^3')
        print(M[l-1], 'M_0')
        print(R[l-1], 'km')
        print('-------------------------------------------')

    print(f'Initial Pressure {p_i / cgs.MeV_fm_to_km}MeV/fm^3 ~ {p_f / cgs.MeV_fm_to_km}/MeV/fm^3')
    data = np.column_stack((P, M, R))
    np.savetxt(f'../Result_2/eft_pnm32_{j:0>6}_ldm_result.dat', data, fmt='%.6E', delimiter=' ', 
                header='Pressure, Mass, Radius', comments='')
    