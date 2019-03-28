import numpy as np
import matplotlib.pyplot as plt

def readAllTimes(file_name, dists_name=[], typ=4):

	#typ = 4 [perturbation, station, 0 - ballistic/ 1 - fragmentation, frag number (0 for ballistic)]
	#typ = 6 [perturbation, station, 0 - ballistic/ 1 - fragmentation, zeniths, velocity] 
	zenith_list = [5, 25, 45, 65, 85]
	velocity_list = [11, 16, 21, 26, 31]
	lat_list = [49.5, 50.5, 49.0, 49.5, 49.5]
	lon_list = [10.0, 10.0, 10.0, 9.5, 10.5]
	titles = ['lat: 49.5 lon: 10.0', 'lat: 50.5 lon: 10.0', 'lat: 49.0 lon: 10.0', 'lat: 49.5 lon: 9.5', 'lat: 49.5 lon: 10.5']

	allTimes = np.load(file_name)
	allDists = np.load(dists_name)

	
	if typ == 6:
		for stn in range(5):

			fig = plt.figure()
			ax = fig.add_subplot(111)
			print('//////////////////////////')
			for ze_i, ze in enumerate(zenith_list):
				print("###############")
				print("Ze: {:} Lat: {:} Lon: {:}".format(ze, lat_list[stn], lon_list[stn]))
				for j in range(len(velocity_list)):
					print("V = {:} km/s".format(velocity_list[j]))
					print(allTimes[0, stn, 0, ze_i, j])
					print(allDists[0, stn, 0, ze_i, j])
				ax.plot(velocity_list, allTimes[0, stn, 0, ze_i, :5], label='Zenith: {:}'.format(ze))
				ax.set_xlabel('Velocity')
				ax.set_ylabel('Arrival Time')
			
			ax.set_title(titles[stn])
			plt.legend()
			plt.show()


if __name__ == '__main__':
	file_name = '/home/luke/Desktop/Seismic_data/Theoretical_Fireball/all_pick_times.npy'
	dists_name = '/home/luke/Desktop/Seismic_data/Theoretical_Fireball/all_pick_dists.npy'
	typ = 6
	readAllTimes(file_name, dists_name=dists_name, typ=typ)