import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, ScalarFormatter

# Read the data from the file
data = np.loadtxt('./R_M_data.txt')

# Separate the columns into variables
p0 = data[:, 0]
M = data[:, 1]
R = data[:, 2]



# Create plots
plt.figure(figsize=(12, 6))

# Plot initial pressure vs mass
plt.subplot(1, 2, 1)
plt.plot(p0, M,label=r'Mass in $M_\odot$')
plt.xlabel(r'Initial Pressure($p_0$)')
plt.ylabel(r'Total Mass($M_\odot$)')
plt.xscale('log')
plt.title('Initial pressure vs Mass')
plt.legend()

# Plot mass vs radius
plt.subplot(1, 2, 2)
plt.plot(p0, R,label='Radius in km', color='red')
plt.xlabel(r'Initial Pressure($p_0$)')
plt.ylabel('Total Radius(km)')
plt.xscale('log')
plt.title('Initial pressure vs Radius')
plt.legend()

# Show plots
plt.tight_layout()
plt.savefig('p0_M_R.png')
plt.show()


# Create the First Plot
fig, ax1 = plt.subplots(figsize=(10, 6))
xticks = [10**i for i in range(30, 43)]

plt.title('Radius and Mass vs. Initial Pressure')
ax1.set_ylabel('Radius in km') 
ax1.plot(p0, R, '*', color='r', markersize=1, label='Radius in km')
ax1.tick_params(axis='y')
ax1.set_xscale('log')
ax1.legend(loc='upper left')  # Adjust legend location
ax1.set_xticks(xticks)
ax1.xaxis.set_major_formatter(ScalarFormatter())
ax1.xaxis.set_major_locator(FixedLocator(xticks))

ax2 = ax1.twinx()
ax2.set_xscale('log')
ax2.set_xlabel(r'Initial Pressure ($dyne/cm^2$)')
ax2.set_ylabel(r'Mass in $M_\odot$')
ax2.plot(p0, M, '*', color='b', markersize=1, label='Mass')
ax2.tick_params(axis='y')
ax2.legend(loc='upper right')  # Adjust legend location

# Set xtick labels
ax1.set_xticklabels([f'$10^{{{i}}}$' for i in range(30, 43)])

plt.savefig('p0_M_R_.png')
plt.show()