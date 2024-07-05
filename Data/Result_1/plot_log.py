import numpy as np
import matplotlib.pyplot as plt


file_name = 'eft_pnm32_000078_ldm_result.dat'
data1 = np.loadtxt(file_name)


p0 = data1[:, 0]
M = data1[:, 1]
R = data1[:, 2]



# Create the First Plot
fig, ax1 = plt.subplots(figsize=(10, 6))
xticks = [10**i for i in range(0, 4)]

plt.title('Radius, Mass vs. Initial Pressure')
ax1.set_ylabel('Radius in km') 
ax1.plot(p0, R, '--', color='k', label=r'Radius')
ax1.tick_params(axis='y')
ax1.set_xscale('log')
ax1.legend(loc='upper left')  # Adjust legend location
ax1.set_xlabel(r'Initial Pressure $MeV/fm^3$')
ax1.set_xticks(xticks)
ax1.set_xlim([1.e0, 1.e4])

ax2 = ax1.twinx()
ax2.set_xscale('log')
ax2.set_xlabel(r'Initial Pressure ($dyne/cm^2$)')
ax2.set_ylabel(r'Mass in $M_\odot$')
ax2.plot(p0, M, color='k', label=r'Mass')

ax2.tick_params(axis='y')
ax2.legend(loc='upper right')


plt.show()