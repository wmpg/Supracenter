import numpy as np
import matplotlib.pyplot as plt

from supra.Fireballs.Program import configRead, configParse, position, station

def readAllTimes(file_name, dists_name=[], typ=4, stn_list=[]):

	#typ = 4 [perturbation, station, 0 - ballistic/ 1 - fragmentation, frag number (0 for ballistic)]
	#typ = 6 [perturbation, station, 0 - ballistic/ 1 - fragmentation, zeniths, velocity] 
	zenith_list = [5, 45, 85]
	velocity_list = [11, 21, 31]
	# lat_list = [49.5, 50.5, 49.0, 49.5, 49.5]
	# lon_list = [10.0, 10.0, 10.0, 9.5, 10.5]
	# titles = ['lat: 49.5 lon: 10.0', 'lat: 50.5 lon: 10.0', 'lat: 49.0 lon: 10.0', 'lat: 49.5 lon: 9.5', 'lat: 49.5 lon: 10.5']

	allTimes = np.load(file_name)
	allDists = np.load(dists_name)

	
	if typ == 6:
		for ii, stn in enumerate(range(stn_list)):

			fig = plt.figure()
			ax = fig.add_subplot(111)
			print('//////////////////////////')
			for ze_i, ze in enumerate(zenith_list):
				print("###############")
				print("Station: {:}".format(stn))
				for j in range(len(velocity_list)):
					print("V = {:} km/s".format(velocity_list[j]))
					print(allTimes[0, ii, 0, ze_i, j])
					print(allDists[0, ii, 0, ze_i, j])
				ax.plot(velocity_list, allTimes[0, ii, 0, ze_i, :], label='Zenith: {:}'.format(ze))
				ax.set_xlabel('Velocity')
				ax.set_ylabel('Arrival Time')
			
			ax.set_title(stn.code)
			plt.legend()
			plt.show()


if __name__ == '__main__':
	file_name = '/home/luke/Desktop/Seismic_data/Theoretical_Fireball/all_pick_times.npy'
	dists_name = '/home/luke/Desktop/Seismic_data/Theoretical_Fireball/all_pick_dists.npy'
	typ = 6

	arg_parser.add_argument('input_file', type=str, help='Path to Supracenter input file.')

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    #################

    setup = configRead(cml_args.input_file)
    configParse(setup, 'picks')

    ##########################################################################################################
    if not os.path.exists(setup.working_directory):
        os.makedirs(setup.working_directory)

        #Build seismic data path
    dir_path = os.path.join(setup.working_directory, setup.fireball_name)

    # Load the station and waveform files list
    data_file_path = os.path.join(dir_path, DATA_FILE)
    if os.path.isfile(data_file_path):
        
        stn_list = readStationAndWaveformsListFile(data_file_path, rm_stat=setup.rm_stat)

    else:
        print('Station and waveform data file not found! Download the waveform files first!')
        sys.exit()
	readAllTimes(file_name, dists_name=dists_name, typ=typ, stn_list=stn_list)