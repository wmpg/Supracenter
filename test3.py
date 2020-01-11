import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import multiprocessing

from supra.Utils.pso import pso

def func_error(Z, *args):
	p_ratio = args
	return np.abs(func(Z[0]) - p_ratio[0])

def gunc_error(W, *args):
	p_ratio, f_d, d, W_0, k, b, I, J_m, P_0, c, c_m, h, h_0, cf = args
	return np.abs(gunc(W[0], f_d, d, W_0, k, b, I, J_m, P_0, c, c_m, h, h_0, cf) - p_ratio)

def inverse_func(p_ratio):
	a, b = pso(func_error, [1], [1000000000], args=([p_ratio]), processes=multiprocessing.cpu_count())

	return a[0]

def inverse_gunc(p_ratio, f_d, d, W_0, k, b, I, J_m, P_0, c, c_m, h, h_0, c_f):
	a, b = pso(gunc_error, [8], [15], args=([p_ratio, f_d, d, W_0, k, b, I, J_m, P_0, c, c_m, h, h_0, cf]), processes=multiprocessing.cpu_count())

	return 10**a[0]

def func(Z):
	# KG85
	return 3.2e6*Z**(-3)*(1 + (Z/87)**2)**(0.5)*(1+Z/800)

def chem_func(Z):
	return 808*(1 + (Z/4.5)**2)/(1 + (Z/0.048)**2)**0.5/(1 + (Z/0.32)**2)**0.5/(1 + (Z/1.35)**2)**0.5

def gunc(W, f_d, d, W_0, k, b, I, J_m, P_0, c, c_m, h, h_0, cf):
	W = 10**W

	v = 1/2/J_m*(W_0*P_0/W/P_a)**(1/3)*(c/c_m)
	#return func(f_d*d/(W/W_0)**(1/3))*P_a*cf*integration_full(k, v, b, I)
	return func(f_d*d/(W/W_0)**(1/3))*P_a*cf*integration(k, v, b, P_0, h, d, h_0=h_0)

def sin_theta(h, l):
	return h/(np.sqrt(h**2 + l**2)) 

def integration(k, v, b, P_0, h, d, h_0=0):
	#attenuation
	return np.exp(-k*v**2/b/P_0/h*d*(np.exp(b*h) - np.exp(b*h_0)))

def integration_full(k, v, b, I, P_a):
	return np.exp(-k*v**2/b/P_a*I)

def phi(f_d, d, W, W_0):
	#scaled distance
	return f_d*d/(W/W_0)**(1/3)

P_a = 101325

z = np.linspace(np.log10(100), np.log10(20000000), num=100)
Z = 10**z

p_ratio = chem_func(Z)*P_a

W_0 = 4.2e12

plt.loglog()
# plt.plot(KG85[:, 0], KG85[:, 1])
# plt.plot(KGSCALED[:, 0], KGSCALED[:, 1])
plt.plot(Z, p_ratio, label='KG85 1kT')
plt.plot(Z*1000**(1/3), p_ratio, label='KG85 1T')
plt.plot(Z*1e6**(1/3), p_ratio, label='KG85 1kg')

# plt.plot(Z/(W_0/W)**(1/3), p_ratio, label='KG85 {:.2E} J'.format((W)))
plt.xlabel('Scaled Distance (m/kg^1/3)')
plt.ylabel('Pressure Ratio (del_p/P_A)')
# plt.legend()
# plt.show()


# plt.loglog()
# # plt.plot(KG85[:, 0], KG85[:, 1])
# # plt.plot(KGSCALED[:, 0], KGSCALED[:, 1])
# plt.plot(p_ratio, Z, label='KG85 1kT')
# plt.ylabel('Actual Distance (m/kg^1/3)')
# plt.xlabel('Pressure Ratio (del_p/P_A)')
# plt.axvline(x=0.28/88800)
# plt.legend()
# plt.show()

P_0 = 868.229
f_d = 0.573
R = 89186
h = 32400
h_0 = 1098
c = 295

k = 2e-4
b = 1.19e-4
J_m = 0.375
c_m = 347

v = 1
I = 113#0.0815
#I = 0.004

geometric_attenuation_angle = 3.07

cf = np.sqrt(geometric_attenuation_angle/5)

p_ratio = chem_func(Z)*P_a*integration(k, v, b, P_0, h, Z, h_0=h_0)

plt.loglog()
# plt.plot(KG85[:, 0], KG85[:, 1])
# plt.plot(KGSCALED[:, 0], KGSCALED[:, 1])
plt.plot(Z, p_ratio, label='KG85 1kT Attenuation Correction')
plt.plot(Z*1000**(1/3), p_ratio, label='KG85 1T Attenuation Correction')
plt.plot(Z*1e6**(1/3), p_ratio, label='KG85 1kg Attenuation Correction')

p_ratio = chem_func(Z)*P_a*cf
plt.plot(Z, p_ratio, label='KG85 1kT Geometric Correction')
plt.plot(Z*1000**(1/3), p_ratio, label='KG85 1T Geometric Correction')
plt.plot(Z*1e6**(1/3), p_ratio, label='KG85 1kg Geometric Correction')

p_ratio = chem_func(Z)*P_a*cf*integration(k, v, b, P_0, h, Z, h_0=h_0)
plt.plot(Z, p_ratio, label='KG85 1kT Combined Correction')
plt.plot(Z*1000**(1/3), p_ratio, label='KG85 1T Geometric Correction')
plt.plot(Z*1e6**(1/3), p_ratio, label='KG85 1kg Geometric Correction')

#plt.plot(Z/(W_0/W)**(1/3), p_ratio, label='KG85 {:.2E} J att'.format((W)))
plt.xlabel(r'Scaled Distance, Z (m kg^{-1/3})')
plt.ylabel(r'Overpressure, $\Delta$p (Pa)')
# plt.title("Effects of the Attenuation on Overpressure for Scaled Distances")
plt.axhline(y=0.28)
plt.gca().set_ylim(bottom=1e-3, top=1e4)
plt.legend()
plt.show()


w = np.linspace(np.log10(1e4), np.log10(1e12))
W = 10**w

d = R
P_s = 88800
# v = 1/2/J_m*(W_0*P_0/W/P_a)**(1/3)*(c/c_m)

# p_nuc = func(phi(f_d, d, W, W_0))*P_s
# p_chem = chem_func(phi(f_d, d, W, W_0))*P_s

plt.loglog()
plt.plot(W, phi(f_d, d, W, W_0))
plt.plot(W, phi(f_d, d, W, W_0))
plt.show()


# p_nuc = func(Z)
# p_chem = chem_func(Z)

# plt.loglog()
# plt.plot(Z, p_nuc)
# plt.plot(Z, p_chem)
# plt.show()
exit()

p_ratio = chem_func(phi(f_d, d, W, W_0))*P_a
plt.loglog()
plt.plot(p_ratio, W, label='Unattenuated R = {:2} km'.format(R/1000))

p_ratio = chem_func(phi(f_d, d, W, W_0))*P_a*cf
plt.plot(p_ratio, W, label='Geometric Factor Only R = {:2} km'.format(R/1000))

p_ratio = chem_func(phi(f_d, d, W, W_0))*P_a*(integration_full(k, v, b, I, P_s))
plt.plot(p_ratio, W, label='Attenuated Factor Only R = {:2} km'.format(R/1000))

p_ratio = chem_func(phi(f_d, d, W, W_0))*P_a*(integration_full(k, v, b, I, P_s))*cf
plt.plot(p_ratio, W, label='Combined Factor R = {:2} km'.format(R/1000))

plt.axvline(x=0.284, c='k')
plt.gca().set_xlim(left=1e-3, right=1e2)
plt.ylabel('Yield, W (J)')
plt.xlabel('Overpressure, del_p (Pa)')
plt.title("Effects of the Attenuation on Yield for a Given Overpressure")
plt.legend()
plt.show()

plt.loglog()
# 5 Frags

color = ['k', 'c', 'r', 'g', 'm']
h = np.array([21200, 22200, 25700, 30100, 32400])
P_0 = np.array([4687, 4142, 2353, 1218, 868])
c = np.array([293.48, 293.27, 290.59, 290.18, 295])
R = np.array([82800, 83342, 85214, 87818, 89186])
I = np.array([113, 113, 113, 113, 113])
p = np.array([0.094975, 0.057085, 0.042187, 0.102567, 0.284065])
cf = np.array([1/3.3333, 1/3.1623, 1/2.3570, 1/1.4086, 1/1.2762])
sol = np.array([5.97e8, 5.89e8, 1.58e9, 1.048e10, 4.26e10])
for i in range(5):
	v = 1/2/J_m*(W_0*P_0[i]/W/P_a)**(1/3)*(c[i]/c_m)
	p_ratio = chem_func(phi(f_d, R[i], W, W_0))*P_s*(integration_full(k, v, b, I[i], P_a))*cf[i]
	plt.plot(p_ratio, W, label='Fragmentation {:d}, Height {:.2f} km'.format(i+1, h[i]/1000), c=color[i])
	plt.axvline(x=p[i], c=color[i])
	plt.text(p[i], 1e12/(i+1), "{:.2f} Pa".format(p[i]))
	plt.scatter(p[i], sol[i], c=color[i])
	plt.text(p[i], sol[i], "{:.2E} J".format(sol[i]))

plt.text(1e1, 1e8, "Total Yield: {:.2E} J".format(np.sum(sol)))
plt.text(1e1, 1.5e8, "1/2 m v^2 Yield: {:.2E} J".format(5.88E10))
plt.gca().set_xlim(left=1e-3, right=1e2)
plt.ylabel('Yield, W (J)')
plt.xlabel('Overpressure, del_p (Pa)')
plt.title("Yield for a Given Overpressure (Attenuation and Geometric Factors)")
plt.legend()
plt.show()

plt.loglog()
color = ['k', 'c', 'r', 'g', 'm']
h = np.array([21200, 22200, 25700, 30100, 32400])
P_0 = np.array([4687, 4142, 2353, 1218, 868])
c = np.array([293.48, 293.27, 290.59, 290.18, 295])
R = np.array([82800, 83342, 85214, 87818, 89186])
I = np.array([0.003096, 0.003813, 0.0100, 0.03748, 0.0815])
p = np.array([0.094975, 0.057085, 0.042187, 0.102567, 0.284065])
cf = np.array([1/3.3333, 1/3.1623, 1/2.3570, 1/1.4086, 1/1.2762])
sol_1 = np.array([3.29e8, 3.68e8, 1.16e9, 9.23e9, 3.81e10])
for i in range(5):
	v = 1/2/J_m*(W_0*P_0[i]/W/P_a)**(1/3)*(c[i]/c_m)
	p_ratio = chem_func(phi(f_d, R[i], W, W_0))*P_a*(integration_full(k, v, b, I[i], P_s))
	plt.plot(p_ratio, W, label='Fragmentation {:d}, Height {:.2f} km'.format(i+1, h[i]/1000), c=color[i])
	plt.axvline(x=p[i], c=color[i])
	plt.text(p[i], 1e12/(i+1), "{:.2f} Pa".format(p[i]))
	plt.scatter(p[i], sol[i], c=color[i])
	plt.text(p[i], sol[i], "{:.2E} J".format(sol_1[i]))

plt.text(1e1, 1e8, "Total Yield: {:.2E} J".format(np.sum(sol_1)))
plt.text(1e1, 1.5e8, "1/2 m v^2 Yield: {:.2E} J".format(5.88E10))
plt.gca().set_xlim(left=1e-3, right=1e2)
plt.ylabel('Yield, W (J)')
plt.xlabel('Overpressure, del_p (Pa)')
plt.title("Yield for a Given Overpressure (Attenuation Factor Only)")
plt.legend()
plt.show()
# ##### HEIGHT AND RANGE ESTIMATES ##################
# p_ratio = np.linspace(np.log10(1e-3), np.log10(1e2))
# p = 10**p_ratio
# plt.loglog()

# for i in range(4):
# 	d = 25000*(i + 1)
# 	W = []
# 	for pp in p:
# 		print(pp, f_d, d, W_0, k, b, I, J_m, P_0, c, c_m, h, h_0, cf)
# 		W.append(inverse_gunc(pp, f_d, d, W_0, k, b, I, J_m, P_0, c, c_m, h, h_0, cf))
	

# 	plt.plot(p_ratio, np.array(W), label='{:2} km Range'.format(d/1000))


# plt.ylabel('Yield, W (J)')
# plt.xlabel('Overpressure, p, (Pa)')
# plt.title("Effects of the Attenuation on Yield given Overpressure")
# plt.legend()
# plt.show()

# plt.loglog()
# d = R
# for i in range(5):
# 	h = 2500*(i) + 22000
# 	W = []
# 	for pp in p:
# 		W.append(inverse_gunc(pp, f_d, d, W_0, k, b, I, J_m, P_0, c, c_m, h, h_0, cf))
	

# 	plt.plot(p_ratio, np.array(W), label='{:2} km Height'.format(h/1000))


# plt.ylabel('Yield, W (J)')
# plt.xlabel('Overpressure, p, (Pa)')
# plt.title("Effects of the Attenuation on Yield given Overpressure")
# plt.legend()
# plt.show()
# #########################
