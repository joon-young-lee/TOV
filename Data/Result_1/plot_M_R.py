import numpy as np
import matplotlib.pyplot as plt


file_name = 'eft_pnm32_000078_ldm_result.dat'
data1 = np.loadtxt(file_name)



M = data1[:, 1]
R = data1[:, 2]


plt.figure(figsize=(10, 8))
plt.title('Radius vs. Mass')
plt.plot(R, M)
plt.xlabel('Radius in km')
plt.ylabel(r'Mass in $M_0$')
plt.xlim(8.0, 30.0)
plt.ylim(0, 2.5)

plt.grid()
plt.show()