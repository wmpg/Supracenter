import numpy as np
import matplotlib.pyplot as plt



TAU_FILE = "F:\\Documents\\Meteor_Research\\Event\\Romania\\Romania_tau.csv"


with open(TAU_FILE, "r+") as f:
	lines = f.readlines()


h = []
t = []
for line in lines:
	line.strip("\n")
	a = line.split(",")
	t.append(float(a[0]))
	h.append(float(a[1].strip("\n")))

t = np.array(t)
h = np.array(h)

heights = np.array([42.125, 47.775, 55.600])
taus = np.array([10.1, 1.8, 4.2])

tau_mins = np.array([9.77, 1.77, 4.13])
tau_maxs = np.array([10.1, 1.85, 4.27])

heights_opt = np.array([42.7, 47.0, 54.4])
taus_opt = np.array([12.58, 2.2, 10.2])

tau_opt_mins = np.array([12.18, 2.16, 10.04])
tau_opt_maxs = np.array([13.00, 2.26, 10.37])

h_mins = np.array([41.950, 47.650, 55.500])
h_maxs = np.array([42.300, 47.900, 55.700])

h_opt_mins = np.array([42.525, 46.875, 54.300])
h_opt_maxs = np.array([42.875, 47.125, 54.500])


h_ball = np.array([38.535])


R = 18
P = 2.879480535097855*100
v = 27760


E_l = R**2*P*np.pi
E_t = E_l*v


I_curve = 2.56E+10
I_mod = 1500*10**(-17.371/-2.5)

tau_ball = np.array([(I_mod/E_t*100)])
print(E_t, I_mod, I_mod/E_t)
plt.style.use('ggplot')
# plt.rcParams['text.usetex'] = True


plt.plot(t, h, label="Borovicka Luminous Efficiency Model")
plt.errorbar(taus, heights, yerr=np.array([abs(h_mins-heights), abs(h_maxs-heights)]), xerr=np.array([abs(tau_mins-taus), abs(tau_maxs-taus)]), fmt='o', label="Acoustic Luminous Efficiency - Acoustic Height", color="g")
plt.errorbar(taus_opt, heights_opt, yerr=np.array([abs(h_opt_mins-heights_opt), abs(h_opt_maxs-heights_opt)]), xerr=np.array([abs(tau_opt_mins-taus_opt), abs(tau_opt_maxs-taus_opt)]), fmt='o', label="Acoustic Luminous Efficiency - Optical Height", color="k")
# plt.scatter(tau_ball, h_ball, label="Ballistic")
plt.ylabel("Height [km]")
plt.xlabel("Luminous Efficiency [%]")
plt.ylim([30, 60])
plt.legend()
plt.show()