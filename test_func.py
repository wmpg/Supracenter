import numpy as np
import matplotlib.pyplot as plt

nominal = np.array([1.60e10, 4.1e9, 9.8e8, 2.01e9, 1.70e10])
max_arr = np.array([1.62e10, 4.2e9, 1.0e9, 2.10e9, 1.77e10])
min_arr = np.array([1.58e10, 4.0e9, 9.5e8, 1.96e9, 1.64e10])

h_nom = np.array([21.2, 22.2, 25.7, 30.1, 32.4])
h_max = np.array([22.2, 23.8, 26.5, 32.2, 34.3])
h_min = np.array([20.1, 20.6, 24.3, 28.9, 31.3])

plt.semilogy()
plt.xlabel('Height of Fragmentation [km]')
plt.ylabel('Energy Released [J]')
plt.scatter(h_nom, nominal, c='g', label='Energy Estimates (This Work)')
plt.errorbar(h_nom, nominal, c='g', xerr=[h_nom-h_min, h_max-h_nom], yerr=[nominal-min_arr, max_arr-nominal], fmt='none')
plt.axvline(x=20.5, c='k')
plt.axvline(x=21.5, c='k')
plt.axvline(x=25.5, c='k')
plt.axvline(x=30.4, c='k', label='Light Curve Peaks (Spurny)')
plt.legend()
plt.show()
