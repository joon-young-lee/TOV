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
    
    dp_dh = lambda p: EoS(p) + p 
    dr2_dh = lambda r2, p, m: -2 * r2 * (np.sqrt(r2) - 2 * m)/(m + 4 * np.pi * p * r2 ** (3.0/2.0))
    dm_dh = lambda r2, p, m: -4 * np.pi * EoS(p) * r2 ** (3.0/2.0) * \
        (np.sqrt(r2) - 2 * m)/(m + 4 * np.pi * p * r2 ** (3.0/2.0))
    

    p = p0
    r2 = -3 / (2 * np.pi * (EoS(p0) + 3 * p0)) * del_h
    m = 0.0
    # Runge-Kutta 4th order
    # print(r2)
    cut = int(300)
    for j in range(cut):
            
        
        k1_p = dp_dh(p)
        k1_m = dm_dh(r2, p, m)
        k1_r2 = dr2_dh(r2, p, m)
        
        
        p_2 = p + del_h * k1_p / 2
        m_2 = m + del_h * k1_m / 2
        r_2 = r2 + del_h * k1_r2 / 2
        k2_p = dp_dh(p_2)
        k2_m = dm_dh(r_2, p_2, m_2)
        k2_r2 = dr2_dh(r_2, p_2, m_2)
        
        
        p_3 = p + del_h * k2_p / 2
        m_3 = m + del_h * k2_m / 2
        r_3 = r2 + del_h * k2_r2 / 2
        k3_p = dp_dh(p_3)
        k3_m = dm_dh(r_3, p_3, m_3)
        k3_r2 = dr2_dh(r_3, p_3, m_3)
        
        
        p_4 = p + del_h * k3_p
        m_4 = m + del_h * k3_m
        r_4 = r2 + del_h * k3_r2
        k4_p = dp_dh(p_4)
        k4_m = dm_dh(r_4, p_4, m_4)
        k4_r2 = dr2_dh(r_4, p_4, m_4)
        

        p += (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_h
        # if i  % 100 == 0:
        #     print(p, 'pressure')
        m += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_h
        # print(m, 'mass')
        r2 += (k1_r2 + 2 * k2_r2 + 2 * k3_r2 + k4_r2) / 6 * del_h
    
    i = cut
    
    while i < maxit and p >= 4.445229102499874e-13: # final pressure 3.3587e-5 MeV/fm^3
        
        
        k1_p = dp_dh(p)
        k1_m = dm_dh(r2, p, m)
        k1_r2 = dr2_dh(r2, p, m)
        
        
        p_2 = p + del_h * k1_p / 2
        m_2 = m + del_h * k1_m / 2
        r_2 = r2 + del_h * k1_r2 / 2
        k2_p = dp_dh(p_2)
        k2_m = dm_dh(r_2, p_2, m_2)
        k2_r2 = dr2_dh(r_2, p_2, m_2)
        
        
        p_3 = p + del_h * k2_p / 2
        m_3 = m + del_h * k2_m / 2
        r_3 = r2 + del_h * k2_r2 / 2
        k3_p = dp_dh(p_3)
        k3_m = dm_dh(r_3, p_3, m_3)
        k3_r2 = dr2_dh(r_3, p_3, m_3)
        
        
        p_4 = p + del_h * k3_p
        m_4 = m + del_h * k3_m
        r_4 = r2 + del_h * k3_r2
        k4_p = dp_dh(p_4)
        k4_m = dm_dh(r_4, p_4, m_4)
        k4_r2 = dr2_dh(r_4, p_4, m_4)
        

        p += (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_h
        # if i  % 100 == 0:
        #     print(p, 'pressure')
        m += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_h
        # print(m, 'mass')
        r2 += (k1_r2 + 2 * k2_r2 + 2 * k3_r2 + k4_r2) / 6 * del_h
        
        # print(r2, 'Radius squared')
        i += 1
    
    print(i)
    return m, np.sqrt(r2)

# 1.e32 ~ 1.e41, divide each interval by 1e3, 2e3 ...
# file_name = '../Archive/eft_pnm32_000001_ldm_eos.dat'
# data = np.loadtxt(file_name)

# e_arr = np.flip(data[:, 0] * cgs.MeV_fm_to_km)
# p_arr = np.flip(data[:, 1] * cgs.MeV_fm_to_km)
# print(e_arr, 'E_D')
# print(p_arr, 'pressure')



del_h = -1.e-4
maxit = int(1e6)  

# p0 = 0.001720545994953356 # 1.3E+3MeV/fm^3
# EoS = eos(file_name)
# print(TOV(EoS, p0, del_h, maxit))

# print(max(p_arr), min(p_arr))



for j in range(1, 101):
    EoS = eos(f'../Archive/eft_pnm32_{j:0>6}_ldm_eos.dat')
    step = 1000
    M = np.zeros(step)
    R = np.zeros(step)
    P = np.zeros(step)  
    
    p_i = 1.3234969191948892e-6 # 1.0 MeV/fm^3
    p_f = 1.3234969191948892e-4 # 100.0 MeV/fm^3
    for l in range(1, step + 1):
        P[l-1] = p_i + (l-1) * (p_f-p_i)/step
        M[l-1], R[l-1] = TOV(EoS, P[l-1], del_h, maxit)
        print('-------------------------------------------')
        print(f'Archive/eft_pnm32_{j:0>6}_ldm_eos.dat: {l/10}%')
        print(f'Initial pressure: {P[l-1] / cgs.MeV_fm_to_km} MeV/fm^3')
        print(M[l-1] * cgs.km_to_M0, 'M_0')
        print(R[l-1], 'km')
        print('-------------------------------------------')

    M *= cgs.km_to_M0
    P *= 1/cgs.MeV_fm_to_km
    print(f'Initial Pressure {p_i / cgs.MeV_fm_to_km}MeV/fm^3 ~ {p_f / cgs.MeV_fm_to_km}MeV/fm^3')
    data = np.column_stack((P, M, R))
    np.savetxt(f'../Result_2/eft_pnm32_{j:0>6}_ldm_result.dat', data, fmt='%.6E', delimiter=' ', 
                header='', comments='')
    