import numpy as np


for j in range(1, 101):
    file_name_2 = f'Result_2/eft_pnm32_{j:0>6}_ldm_result.dat'
    file_name_3 = f'Result_3/eft_pnm32_{j:0>6}_ldm_result.dat'
    data2 = np.loadtxt(file_name_2)
    data3 = np.loadtxt(file_name_3)
    P1 = data2[:, 0]
    M1 = data2[:, 1]
    R1 = data2[:, 2]
    P2 = data3[:, 0]
    M2 = data3[:, 1]
    R2 = data3[:, 2]
    M = np.concatenate((M1, M2))
    R = np.concatenate((R1, R2))
    P = np.concatenate((P1, P2))
    file_name = f'Result_metric/eft_pnm32_{j:0>6}_ldm_result.dat'
    data = np.column_stack((P, M, R))
    np.savetxt(file_name, data, fmt='%.6E', delimiter=' ', header='', comments='')
