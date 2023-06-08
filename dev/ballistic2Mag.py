import numpy as np


tau = 0.08
R_0 = 60
v = 21000
P = 100
D = 100 #km
L = 1500 #m

R = tau*v*R_0**2*P/4/np.pi


print("Energy: {:.2E} J/m".format(R_0**2*P))
print("Energy: {:.2E} J/s".format(R_0**2*P*v))


print("Radiance: {:.2E} Watts/str".format(R/tau))

print("Length of Region: {:.2E} m".format(L))

E = R*L
print("Total Energy: {:.2E}".format(E))

#convert to cm
D_cm = D*1000*100

M = -2.5*np.log10(E/D_cm**2/2.48e-12)

print(M)
