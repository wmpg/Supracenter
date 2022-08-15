import numpy as np
import matplotlib.pyplot as plt

# FILE_NAME = "F:\\Desktop\\trajectory_samples.csv"

FILE_NAME = "F:\\Desktop\\NZ_fragmentation.txt"

data = []

with open(FILE_NAME, "r+") as f:

    for ll, line in enumerate(f):
            
        # Header
        if ll < 0:
            continue

        line = line.split(',')

        new_line = []
        for element in line:
            try:
                element = float(element)
            except:
                element = np.nan

            new_line.append(element)

        data.append(new_line)


data = np.array(data)


#temp_traj.t, temp_traj.v, temp_traj.zenith.deg, temp_traj.azimuth.deg, temp_traj.pos_f.lat, temp_traj.pos_f.lon, total_error, failed_stats

# plt.subplots(231)
plt.scatter(data[:, 1], data[:, 0])
plt.show()


