import numpy as np
import matplotlib.pyplot as plt

# data = np.array([
# [1,	21200,	20100,	22200,	663.0,	0.094975,	4687.744,	5560.057,	4017.002,	0.072534,	0.086634,	0.061750,	82800,	82274,	83342,	5.31E+09,	5.25E+09,	5.38E+09,	5.52E+08,	5.42E+08,	5.63E+08,	5.66E+09,	5.55E+09,	5.77E+09,	7.69E+09,	7.64E+09,	7.75E+09],
# [2,	22200,	20800,	23800,	398.5,	0.057085,	4017.002,	4987.542,	3140.440,	0.061750,	0.077376,	0.047770,	83342,	82625,	84166,	3.24E+09,	3.18E+09,	3.30E+09,	1.58E+08,	1.54E+08,	1.62E+08,	3.47E+09,	3.38E+09,	3.57E+09,	4.66E+09,	4.61E+09,	4.71E+09],
# [3,	25700,	24300,	26500,	294.5,	0.042187,	2353.861,	2909.731,	2084.172,	0.035359,	0.044115,	0.031139,	85214,	84453,	85604,	2.50E+09,	2.45E+09,	2.52E+09,	7.92E+07,	7.71E+07,	8.03E+07,	2.74E+09,	2.67E+09,	2.78E+09,	3.52E+09,	3.49E+09,	3.54E+09],
# [4,	30100,	28900,	32200,	716.0,	0.102567,	1218.085,	1455.144,	893.7981,	0.017763,	0.021393,	0.012826,	87818,	87101,	89072,	6.45E+09,	6.35E+09,	6.62E+09,	7.99E+08,	7.79E+08,	8.33E+08,	7.29E+09,	7.12E+09,	7.61E+09,	8.85E+09,	8.77E+09,	8.99E+09],
# [5,	32400,	31300,	34300,	1983.0,	0.284065,	868.2219,	1019.529,	659.5404,	0.012424,	0.014745,	0.009180,	89186,	88545,	90374,	1.84E+10,	1.82E+10,	1.86E+10,	1.07E+10,	1.04E+10,	1.11E+10,	2.12E+10,	2.07E+10,	2.20E+10,	2.49E+10,	2.47E+10,	2.53E+10]
# ])

# h_min = data[:, 2]/1000
# h_max = data[:, 3]/1000
# h = data[:, 1]/1000
# h_avg = (h_max + h_min)/2
# herr = h_max - h_avg

# Ell_min = data[:, 16]
# Ell_max = data[:, 17]
# Ell = data[:, 15]
# Ell_avg = (Ell_max + Ell_min)/2
# Ell_err = Ell_max - (Ell_max + Ell_min)/2

# Er_min = data[:, 19]
# Er_max = data[:, 20]
# Er = data[:, 18]
# Er_avg = (Er_max + Er_min)/2
# Er_err = Er_max - (Er_max + Er_min)/2

# Esr_min = data[:, 22]
# Esr_max = data[:, 23]
# Esr = data[:, 21]
# Esr_avg = (Esr_max + Esr_min)/2
# Esr_err = Esr_max - (Esr_max + Esr_min)/2

# Esp_min = data[:, 25]
# Esp_max = data[:, 26]
# Esp = data[:, 24]
# Esp_avg = (Esp_max + Esp_min)/2
# Esp_err = Esp_max - (Esp_max + Esp_min)/2

# # plt.scatter(h, Ell, c='r', label='Lutzky & Lehto 1968')
# # #plt.plot(h, Ell, c='r')
# # plt.errorbar(h, Ell, xerr=[h_max-h, h-h_min], yerr=[Ell_max-Ell, Ell-Ell_min], fmt='none', c='r')

# plt.scatter(h, Er, c='r', label='Fragmentation Energy (Reed 1972)')
# #plt.plot(h, Er, c='g')
# plt.errorbar(h, Er, xerr=[h_max-h, h-h_min], yerr=[Er_max-Er, Er-Er_min], fmt='none', c='r')

# # plt.scatter(h, Esr, c='b', label='Sach Scaling (R)')
# # #plt.plot(h, Esr, c='g')
# # plt.errorbar(h, Esr, xerr=[h_max-h, h-h_min], yerr=[Esr_max-Esr, Esr-Esr_min], fmt='none', c='b')

# # plt.scatter(h, Esp, c='m', label='Sach Scaling (P)')
# # #plt.plot(h, Esp, c='m')
# # plt.errorbar(h, Esp, xerr=[h_max-h, h-h_min], yerr=[Esp_max-Esp, Esp-Esp_min], fmt='none', c='m')

plt.xlabel('Height of Fragmentation [km]')
plt.ylabel('Energy Released [J]')
# plt.title('Energy of Fragmentations \n Deposited at Various Altitudes')
h = np.array([21200, 22200, 25700, 30100, 32400])/1000
plt.axvline(x=20.5, c='k')
plt.axvline(x=21.5, c='k')
plt.axvline(x=25.5, c='k')
plt.axvline(x=30.4, c='k', label='Light Curve Peaks (Spurny)')
# plt.legend(loc='upper left', fancybox=True, framealpha=1, edgecolor='k')
# plt.show()
data = np.array([
[8.08E+09,	8.02E+09,	8.14E+09,	5.52E+08,	5.42E+08,	5.63E+08,	2.70E+09,	2.70E+09,	2.71E+09,		5.95E+09,	5.91E+09,	5.99E+09,	2.57E+08,	2.52E+08,	2.62E+08,	2.52E+09,	2.51E+09,	2.52E+09],
[4.89E+09,	4.85E+09,	4.94E+09,	1.58E+08,	1.54E+08,	1.62E+08,	3.16E+09,	3.16E+09,	3.17E+09,		3.60E+09,	3.57E+09,	3.64E+09,	7.34E+07,	7.15E+07,	7.55E+07,	2.25E+09,	2.25E+09,	2.26E+09],
[3.70E+09,	3.67E+09,	3.72E+09,	7.92E+07,	7.71E+07,	8.03E+07,	3.75E+09,	3.74E+09,	3.76E+09,		2.73E+09,	2.70E+09,	2.74E+09,	3.68E+07,	3.58E+07,	3.73E+07,	2.13E+09,	2.12E+09,	2.13E+09],
[9.30E+09,	9.22E+09,	9.44E+09,	7.99E+08,	7.79E+08,	8.33E+08,	5.66E+09,	5.65E+09,	5.68E+09,		6.84E+09,	6.78E+09,	6.95E+09,	3.71E+08,	3.62E+08,	3.87E+08,	2.60E+09,	2.59E+09,	2.61E+09],
[2.62E+10,	2.60E+10,	2.66E+10,	1.07E+10,	1.04E+10,	1.11E+10,	8.92E+09,	8.90E+09,	8.96E+09,		1.93E+10,	1.91E+10,	1.96E+10,	4.97E+09,	4.86E+09,	5.17E+09,	3.38E+09,	3.38E+09,	3.40E+09]])

sach = data[:, 0]
reed = data[:, 3]
reed_fix = data[:, 6]
sach_cf = data[:, 9]
reed_cf = data[:, 12]
reed_fix_cf = data[:, 15]

print('TOTALS:')
print("Sach: {:.2E} J".format(np.sum(sach)))
print("Reed: {:.2E} J".format(np.sum(reed)))
print("Reed Att: {:.2E} J".format(np.sum(reed_fix)))
print("Sach Latunde: {:.2E} J".format(np.sum(sach_cf)))
print("Reed Latunde: {:.2E} J".format(np.sum(reed_cf)))
print("Reed Att Latunde: {:.2E} J".format(np.sum(reed_fix_cf)))
print("Kinetic Energy: 5.88E+10 J")


plt.scatter(h, sach, label='Sach')
plt.scatter(h, reed, label='Reed')
plt.scatter(h, reed_fix, label='Reed Attenuation')
plt.scatter(h, sach_cf, label='Sach-Latunde')
plt.scatter(h, reed_cf, label='Reed-Latunde')
plt.scatter(h, reed_fix_cf, label='Reed Attenuation-Latunde')
plt.legend(loc='upper left', fancybox=True, framealpha=1, edgecolor='k')
plt.show()


