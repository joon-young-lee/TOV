import numpy as np
import matplotlib.pyplot as plt


# Read the data from the file
data1 = np.loadtxt('./R_M_gp.txt')
data2 = np.loadtxt('./R_M.txt')
# Separate the columns into variables
p0 = data1[:, 0]
M = data1[:, 1]
R = data1[:, 2]
p0_ = data2[:, 0]
M_ = data2[:, 1]
R_ = data2[:, 2]


# Create plots
# plt.figure(figsize=(12, 6))

# # Plot initial pressure vs mass
# plt.subplot(1, 2, 1)
# plt.plot(p0, M,label=r'Mass in $M_\odot$')
# plt.xlabel(r'Initial Pressure($p_0$)')
# plt.ylabel(r'Total Mass($M_\odot$)')
# plt.xscale('log')
# plt.title('Initial pressure vs Mass')
# plt.legend()

# # Plot mass vs radius
# plt.subplot(1, 2, 2)
# plt.plot(p0, R,label='Radius in km', color='red')
# plt.xlabel(r'Initial Pressure($p_0$)')
# plt.ylabel('Total Radius(km)')
# plt.xscale('log')
# plt.title('Initial pressure vs Radius')
# plt.legend()

# # Show plots
# plt.tight_layout()
# plt.savefig('p0_M_R.png')
# plt.show()


# Create the First Plot
fig, ax1 = plt.subplots(figsize=(10, 6))
xticks = [10**i for i in range(29, 43)]

plt.title('Radius and Mass vs. Initial Pressure')
ax1.set_ylabel('Radius in km') 
ax1.plot(p0, R, '--', color='b', 
         label=r'Radius for $dV=\dfrac{4\pi r^2dr}{\sqrt{1-\frac{2Gm}{c^2r}}}$')
ax1.plot(p0_, R_, '--', color='r', 
         label=r'Radius for $dV=4\pi r^2dr$')
ax1.tick_params(axis='y')
ax1.set_xscale('log')
ax1.legend(loc='upper left')  # Adjust legend location
ax1.set_xlabel('Initial Pressure')
ax1.set_xticks(xticks)
ax1.set_xlim([1.e30, 1.e42])

ax2 = ax1.twinx()
ax2.set_xscale('log')
ax2.set_xlabel(r'Initial Pressure ($dyne/cm^2$)')
ax2.set_ylabel(r'Mass in $M_\odot$')
ax2.plot(p0, M, color='b', 
         label=r'Mass for $dV=\dfrac{4\pi r^2dr}{\sqrt{1-\frac{2Gm}{c^2r}}}$')
ax2.plot(p0_, M_, color='r', 
         label=r'Mass $dV=4\pi r^2dr$')
ax2.tick_params(axis='y')
ax2.legend(loc='upper right')  # Adjust legend location
ax2.set_ylim([0, 0.8])

plt.savefig('p0_M_R_2.png')
plt.show()