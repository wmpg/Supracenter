import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

from supra.Fireballs.Program import configRead, configParse, position, station
from supra.Fireballs.GetIRISData import readStationAndWaveformsListFile

def readAllTimes(file_name, dists_name=[], typ=4, stn_list=[]):

    #typ = 4 [perturbation, station, 0 - ballistic/ 1 - fragmentation, frag number (0 for ballistic)]
    #typ = 6 [perturbation, station, 0 - ballistic/ 1 - fragmentation, zeniths, velocity] 
    zenith_list = [5, 45, 85]
    velocity_list = [11, 15, 19, 23, 27, 31]
    # lat_list = [49.5, 50.5, 49.0, 49.5, 49.5]
    # lon_list = [10.0, 10.0, 10.0, 9.5, 10.5]
    # titles = ['lat: 49.5 lon: 10.0', 'lat: 50.5 lon: 10.0', 'lat: 49.0 lon: 10.0', 'lat: 49.5 lon: 9.5', 'lat: 49.5 lon: 10.5']
    c = ['r', 'y', 'b', 'g', 'k']
    allTimes = np.load(file_name)
    #print(allTimes.shape)
    if typ == 6:
        plt.style.use('dark_background')
        for ii, stn in enumerate(stn_list):

            fig = plt.figure()
            ax = fig.add_subplot(111)
            print('//////////////////////////')
            for ze_i, ze in enumerate(zenith_list):
                print("###############")
                print("Station: {:}".format(stn.code))
                for j in range(len(velocity_list)):
                    print("V = {:} km/s".format(velocity_list[j]))
                    print(allTimes[0, ii, 0, ze_i, j])
                    #print(allDists[0, ii, 0, ze_i, j])
                plt.plot(velocity_list, allTimes[0, ii, 0, ze_i, :] - (allTimes[0, ii, 0, ze_i, 2] + allTimes[0, ii, 0, ze_i, 3])/2, label='Zenith: {:}'.format(ze))
                plt.plot([11, 31], [0, 0], 'w--')
            ax.set_xlabel('Average Velocity (km/s)')
            ax.set_ylabel('Relative Arrival Time (s)')
            ax.set_ylim([-2, 2])
            ax.set_title("{:}-{:}".format(stn.network.strip(), stn.code.strip()))
            plt.legend()
        plt.show()


if __name__ == '__main__':
    file_name = '/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball THEO/all_pick_times.npy'
    dists_name = '/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball THEO/all_pick_dists.npy'
    typ = 6
    arg_parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

    arg_parser.add_argument('input_file', type=str, help='Path to Supracenter input file.')

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    #################

    setup = configRead(cml_args.input_file)
    configParse(setup, 'picks')

    picks_name = os.path.join(setup.working_directory, setup.fireball_name, 'data.txt')

    picks_name = setup.station_picks_file
    
    stn_list = []

    with open(picks_name) as f:
        for ii, line in enumerate(f):
            if ii > 0:
                line = line.split(',')
                stn = station(line[1], line[2], position(float(line[3]), float(line[4]), float(line[5])), '...', '...', '...')
                stn_list.append(stn)

    readAllTimes(file_name, dists_name=dists_name, typ=typ, stn_list=stn_list)