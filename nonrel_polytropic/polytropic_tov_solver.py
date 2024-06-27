import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd
import units_cgs as cgs

# All units in cgs

# Using Runge-Kutta 4th order

def TOV(gamma, maxit, del_r, p0): # inputs: gamma, max iteration,  ∆r, initial pressure
    
    
    K = 6.428e-26
    # Polytropic EOS    
    # p = K * e ** gamma
    # e = (p/K) ** (1/gamma)
    # substitue for e in are equations

    # m = mass, e = energy density, r = distance, p = pressure
    # rho = e/c**2, not mass density!

    dp_dr = lambda r, p, m: -cgs.G * (p/K) ** (1/gamma) * m * \
        (1 + p/((p/K) ** (1/gamma))) * \
        (1 + 4 * cgs.pi * p * r ** 3/(m * cgs.c2)) / \
        (1 - 2 * cgs.G * m / (cgs.c2 * r)) / \
        (cgs.c * r) ** 2

    dm_dr = lambda r, p, m: 4 * cgs.pi * (p/K) ** (1/gamma)  * r ** 2 / cgs.c2
        # np.sqrt(1 - 2 * cgs.G * m/r)
    
    r_arr = [1.0e-10]
    p_arr = [p0]
    m_arr = [1.0e-10]
    # Runge-Kutta 4th order
    i = 0
    while i <= maxit and p_arr[i] > 0:
        

        k1_p = dp_dr(r_arr[i], p_arr[i], m_arr[i])
        k1_m = dm_dr(r_arr[i], p_arr[i], m_arr[i])
        
        k2_p = dp_dr(r_arr[i] + del_r/2, p_arr[i] + k1_p/2, m_arr[i] + k1_m/2)
        k2_m = dm_dr(r_arr[i] + del_r/2, p_arr[i] + k1_p/2, m_arr[i] + k1_p/2)
        
        k3_p = dp_dr(r_arr[i] + del_r/2, p_arr[i] + k2_p/2, m_arr[i] + k2_m/2)
        k3_m = dm_dr(r_arr[i] + del_r/2, p_arr[i] + k2_p/2, m_arr[i] + k2_p/2)
        
        k4_p = dp_dr(r_arr[i] + del_r, p_arr[i] + k3_p, m_arr[i] + k3_m)
        k4_m = dm_dr(r_arr[i] + del_r, p_arr[i] + k3_p, m_arr[i] + k3_p)
        
        p_new = p_arr[i] + (k1_p + 2*k2_p + 2*k3_p + k4_p) / 6 * del_r
        m_new = m_arr[i] + (k1_m + 2*k2_m + 2*k3_m + k4_m) / 6 * del_r
        r_new = r_arr[i] + del_r

        p_arr.append(p_new)
        m_arr.append(m_new)
        r_arr.append(r_new)
        
        i += 1
    print(f'Total iterations: {i}')
    return r_arr, p_arr, m_arr, 
    
    
r, p, m = TOV(gamma = 5/3, maxit = int(1.0e7), del_r = 1., p0 = 1.0e32)
del_r = 1.
maxit = int(1.0e7)
# plt.plot(r, p, color = 'b')
# plt.show()
# plt.plot(r, m, color = 'r')
# plt.show()

print(f'Iteration to {maxit * del_r/1.0e5} km')