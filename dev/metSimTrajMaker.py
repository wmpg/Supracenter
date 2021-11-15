import numpy as np

from wmpl.Utils.TrajConversions import date2JD
from wmpl.Trajectory.Trajectory import Trajectory

file_name = "C:\\Users\\lmcfd\\Documents\\Fireballs\\BC\\metsimtraj.csv"

jdt_ref = date2JD(2017, 9, 5, 5, 11, 27, millisecond=0)

station = [50.0635, -117.0855, 36927.4988]
station_id = "GLM-16"

# 1 = Right Ascension for meas1, declination for meas2, NOTE: epoch of date, NOT J2000!
# 2 = Azimuth +east of due north for meas1, Elevation angle above the horizon for meas2
# 3 = Azimuth +west of due south for meas1, Zenith angle for meas2
# 4 = Azimuth +north of due east for meas1, Zenith angle for meas2
meastype = 2

# Init new trajectory solving
traj_solve = Trajectory(jdt_ref, meastype=meastype, save_results=True,\
 monte_carlo=False, show_plots=False, output_dir="C:\\Users\\lmcfd\\Desktop\\denis_output\\Crawford_Bay")


data = []
with open(file_name, "r+") as f:

	for line in f:
		a = line.split(",")
		data.append([float(a[0]), float(a[1]), float(a[2])])

data = np.array(data)


traj_solve.infillTrajectory(np.radians(data[:, 1]), np.radians(data[:, 2]), data[:, 0], \
			np.radians(station[0]), np.radians(station[1]), station[2], station_id=station_id)

traj_solve.run()



