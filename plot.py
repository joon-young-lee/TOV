import numpy as np
import matplotlib.pyplot as plt
import units_cgs as cgs

# Read the data from the file
data = np.loadtxt('tov_data.txt')

# Separate the columns into variables
r = data[:, 0]/1e5
p = data[:, 1]
m = data[:, 2]/cgs.M0

# Create plots
plt.figure(figsize=(12, 6))

# Plot pressure vs radius
plt.subplot(1, 2, 1)
plt.plot(r, p, label=rf'initial presssure: {p[0]} $dyne/cm^2$')
plt.xlabel('r(km)')
plt.ylabel(r'Pressure (p) $dyne/cm^2$')
plt.title('Pressure vs Radius')
plt.legend()

# Plot mass vs radius
plt.subplot(1, 2, 2)
plt.plot(r, m, label=rf'Mass {m[-1]:.5}$M_\odot$' + f'\n Radius {r[-1]}km', color='red')
plt.xlabel('r(km)')
plt.ylabel(r'Mass (m) in $M_{\odot}$')
plt.title('Mass vs Radius')
plt.legend()

# Show plots
plt.tight_layout()
plt.savefig('10^32.png')
plt.show()
