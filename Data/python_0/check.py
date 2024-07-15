import numpy as np
import units_cgs as cgs
k = 2.71396278250104e-24
for i in range(1, 101):
    data = np.loadtxt(f'../Archive/eft_pnm32_{i:0>6}_ldm_eos.dat')
    press = data[:, 1] * cgs.MeV_fm_to_km

    print(np.min(press))
    if k < np.min(press):
        k = np.min(press)

print(k)

