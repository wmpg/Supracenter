import numpy as np
import matplotlib.pyplot as plt

FILE_NAME = "F:\\Documents\\Meteor_Research\\Event\\Alaska\\metsim_results.csv"

time = []
b_ht = []
b_len = []
b_vel = []
l_ht = []
l_len = []
l_vel = []
l_dyn_press = []
m_ht = []
m_len = []
m_vel = []
m_dyn_press = []
tau_total = []
tau_main = []
tau_grain = []
abs_mag_total = []
abs_mag_main = []
abs_mag_grain = []
lum_total = []
lum_main = []
lum_grain = []
elec_line_dens = []
mass_total = []
mass_main = []

data = []

with open(FILE_NAME, "r+") as f:

    for ll, line in enumerate(f):
            
        # Header
        if ll < 2:
            continue

        line = line.split(',')

        new_line = []
        for element in line:
            element = float(element)

            new_line.append(element)

        data.append(new_line)

data = np.array(data)

with open("F:\\Desktop\\Alaska_tau_curve.csv", "w+") as f:
    for ii in range(len(data[:, 12])):

        f.write("{:}, {:}\n".format(data[ii, 12]*100, data[ii, 8]))

plt.plot(data[:, 12]*100, data[:, 8])
plt.ylim([20, 60])
plt.show()
