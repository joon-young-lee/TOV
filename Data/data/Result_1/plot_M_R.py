import numpy as np
import matplotlib.pyplot as plt
# Result_1
# p_i = 1.602179e+33 # 1.0 MeV/fm^3
# p_f = 1.602179e+36 # 1000 MeV/fm^3
# step = 1000

file_name = 'eft_pnm32_000027_ldm_result.dat'
data1 = np.loadtxt(file_name)



M = data1[:, 1]
R = data1[:, 2]


plt.figure(figsize=(10, 8))
plt.title('Radius vs. Mass')
plt.plot(R, M, 'o')
plt.xlabel('Radius in km')
plt.ylabel(r'Mass in $M_0$')
plt.xlim(8.0, 30.0)
plt.ylim(0, 2.5)

plt.grid()
#plt.savefig(f'{file_name}_plot.png')
plt.show()