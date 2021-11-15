""" Tool for marking times of arrival of seismic/air waves in IRIS data. """

from __future__ import print_function, division, absolute_import

import os
import sys
import datetime
import argparse
import pickle

import obspy
import numpy as np
import scipy.signal
import matplotlib
import time
import copy

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button
from matplotlib.widgets import Slider

# # Force backend
# plt.switch_backend("TkAgg")

DATA_FILE = 'data.txt'
OUTPUT_CSV = 'data_picks.csv'


from supra.Fireballs.GetIRISData import readStationAndWaveformsListFile, butterworthBandpassFilter, \
    plotAllWaveforms, convolutionDifferenceFilter
from supra.Fireballs.SeismicTrajectory import timeOfArrival, waveReleasePointWinds, \
    parseWeather, Constants
from supra.Utils.Classes import Position, Station


global arrTimes 
global sounding


def trajCalc(setup):
    """ Creates trajectory between point A and the ground (B) based off of the initial position and the angle of travel

    Arguments:
    setup: [Object] ini file parameters

    Returns: 
    A [list] lat/lon/elev of the tail of the trajectory
    B [list] lat/lon/elev of the head of the trajectory
    """


    B = np.array([0, 0, 0])

    # convert angles to radians
    ze = np.radians(setup.zangle)
    az = np.radians(setup.azim)

    # Create trajectory vector
    traj = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

    # backwards propegate the trajectory until it reaches 100000 m up
    n = 85920/traj[2]

    # B is the intersection between the trajectory vector and the ground
    A = n*traj

    # Convert back to geo coordinates
    B = np.array(loc2Geo(setup.lat_centre, setup.lon_centre, 0, B))
    A = np.array(loc2Geo(setup.lat_centre, setup.lon_centre, 0, A))

    # print("Created Trajectory between A and B:")
    # print("     A = {:10.4f}N {:10.4f}E {:10.2f}m".format(A[0], A[1], A[2]))
    # print("     B = {:10.4f}N {:10.4f}E {:10.2f}m".format(B[0], B[1], B[2]))

    A[2] /= 1000
    B[2] /= 1000

    setup.lat_i = A[0]
    setup.lon_i = A[1]
    setup.elev_i = A[2]

    return A, B


def calcAllTimes(stn_list, setup, sounding):
        """ 
        Method to calculate all arrival times and place them into an array, so that they are not recalculated every
        time a new waveform is opened
                
        Arguments:
        data_list: [ndarray] array containing all station names and locations
        setup: [object] ini file parameters
        sounding: [ndarray] atmospheric profile

        Returns:
        allTimes: [ndarray] a ndarray that stores the arrival times for all stations over every perturbation
        [perturbation, station, 0 - ballistic/ 1 - fragmentation, frag number (0 for ballistic)]
        """

        zenith_list = [5, 45, 85]
        velocity_list = [11, 15, 19, 23, 27, 31]
        ze_array = [0]*len(zenith_list)
        ze_array_dist = [0]*len(zenith_list)
        v_array = [0]*len(velocity_list)
        v_array_dist = [0]*len(velocity_list)

        #All perturbation happens here
        allTimes = [0]*setup.perturb_times
        allDists = [0]*setup.perturb_times

        # Ballistic Prediction
        ref_pos = Position(setup.lat_centre, setup.lon_centre, 0)
        #lat0, lon0, elev0 = data_list[0][2], data_list[0][3], data_list[0][4]
        # if setup.fragmentation_point == '':
        #     setup.fragmentation_point = []

        no_of_frags = len(setup.fragmentation_point)

        # array of frags and ballistic arrivals have to be the same size. So, minimum can be 1
        if no_of_frags == 0:
            no_of_frags = 1

        # Initialize variables
        b_time = 0
        b_dist = 0

        consts = Constants()

        # convert angles to radians
        az = np.radians(setup.azim)
        ze = np.radians(setup.zangle)

        # vector of the trajectory of the fireball
        traj_vect = np.array([np.cos(az)*np.sin(ze), np.sin(az)*np.cos(ze), -np.cos(ze)])

        # For temporal perturbations, fetch the soudning data for the hour before and after the event
        if setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(setup, consts, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(setup, consts, time=-1)

        else:
            sounding_u = []
            sounding_l = []

        if setup.perturb_method == 'ensemble':
            ensemble_file = setup.perturbation_spread_file

        d_time = (setup.perturb_times*len(stn_list)*no_of_frags*len(zenith_list)*len(velocity_list))
        count = 0

        #number of perturbations
        for ptb_n in range(setup.perturb_times-1):

            if ptb_n > 0:
                
                if setup.debug:
                    print("STATUS: Perturbation {:}".format(ptb_n))

                # generate a perturbed sounding profile
                sounding_p = perturb(setup, sounding, setup.perturb_method, \
                    sounding_u=sounding_u, sounding_l=sounding_l, \
                    spread_file=setup.perturbation_spread_file, lat=setup.lat_centre, lon=setup.lon_centre, ensemble_file=ensemble_file, ensemble_no=ptb_n)
            else:

                # if not using perturbations on this current step, then return the original sounding profile
                sounding_p = sounding

            # Initialize station times array
            stnTimes = [0]*len(stn_list)
            stnDists = [0]*len(stn_list)

            #number of stations
            for n, stn in enumerate(stn_list):
                for ze_idx, ZE in enumerate(zenith_list):
                    for v_idx, V in enumerate(velocity_list): 
                        # For ballistic arrivals
                        if setup.show_ballistic_waveform:
                            bTimes = [0]*no_of_frags
                            bDist = [0]*no_of_frags
                            for i in range(no_of_frags):

                                #need filler values to make this a numpy array with fragmentation
                                if i == 0:

                                    setup.zangle = ZE
                                    setup.v = V
                                    A, B = trajCalc(setup)
                                    setup.traj_i = Position(A[0], A[1], A[2])
                                    setup.traj_f = Position(B[0], B[1], B[2])

                                    stn.position.pos_loc(ref_pos)
                                    setup.traj_f.pos_loc(ref_pos)
                                    # Station location in local coordinates
                                    # xx, yy, zz  = geo2Loc(lat0, lon0, elev0, \
                                    #                             stat_lat, stat_lon, stat_elev)
                                    #xx, yy, zz = rotateVector(np.array([xx, yy, zz]), np.array([0, 0, 1]), -np.pi/2)

                                    # Time to travel from trajectory to station
                                    b_time = timeOfArrival([stn.position.x, stn.position.y, stn.position.z], setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*V, \
                                                                np.radians(setup.azim), np.radians(ZE), setup, sounding=sounding_p, travel=False, fast=True, ref_loc=[ref_pos.lat, ref_pos.lon, ref_pos.elev], theo=True)# + setup.t 
                                    S = waveReleasePointWinds([stn.position.x, stn.position.y, stn.position.z], setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*V, \
                                                               np.radians(setup.azim), np.radians(ZE), setup, sounding_p, [ref_pos.lat, ref_pos.lon, ref_pos.elev])
                                   
                                    # S = waveReleasePoint([stn.position.x, stn.position.y, stn.position.z], setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*V, \
                                    #                             np.radians(setup.azim), np.radians(ZE), 310)
                                    # Distance from the source location to the station
                                    #b_dist = ((stn.position.x)**2 + (stn.position.y)**2)**0.5
                                    #bDist[i] = ((S[0] - stn.position.x)**2 + (S[1] - stn.position.y)**2 + (S[2] - stn.position.z)**2)**0.5
                                    bDist[i] = S[2]
                                    bTimes[i] = b_time
                                else:
                                    bDist[i] = np.nan
                                    bTimes[i] = np.nan
                        else:
                            bTimes = [np.nan]*no_of_frags
                            bDist = [np.nan]*no_of_frags
                        v_array_dist[v_idx] = np.array(bDist)
                        v_array[v_idx] = np.array(bTimes)
                    ze_array_dist[ze_idx] = np.array(v_array_dist)
                    ze_array[ze_idx] = np.array(v_array)
                stnDists[n] = ([np.array(ze_array_dist)])
                stnTimes[n] = ([np.array(ze_array)])
            allDists[ptb_n] = np.array(stnDists)
            allTimes[ptb_n] = np.array(stnTimes)

        allTimes = np.array(allTimes)
        allDists = np.array(allDists)

        # Save as .npy file to be reused in SeismicTrajectory and Supracenter
        np.save(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times'), allTimes)
        np.save(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_dists'), allDists)
        print("All Times File saved as {:}".format(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times.npy')))
        print("All Dists File saved as {:}".format(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_dists.npy')))

        return allTimes

if __name__ == "__main__":  

    with open('/home/luke/Desktop/pkl/Theoretical.pkl', 'rb') as f:
        setup = pickle.load(f)

    ##########################################################################################################
    if not os.path.exists(setup.working_directory):
        os.makedirs(setup.working_directory)

    picks_name = setup.station_picks_file
    
    stn_list = []

    with open(picks_name) as f:
        for ii, line in enumerate(f):
            if ii > 0:
                line = line.split(',')
                stn = Station(line[1], line[2], Position(float(line[3]), float(line[4]), float(line[5])), '...', '...', '...')
                stn_list.append(stn)

    # Remove duplicate lines recieved from different sources
    # stn_list = set(tuple(element) for element in stn_list)
    # stn_list = [list(t) for t in set(tuple(element) for element in stn_list)]


    if setup.stations is not None:
        stn_list = stn_list + setup.stations

    # Init the constants
    consts = Constants()
    setup.search_area = [0, 0, 0, 0]
    setup.search_area[0] = setup.lat_centre - setup.deg_radius 
    setup.search_area[1] = setup.lat_centre + setup.deg_radius
    setup.search_area[2] = setup.lon_centre - setup.deg_radius
    setup.search_area[3] = setup.lon_centre + setup.deg_radius

    sounding = parseWeather(setup)

    if setup.weather_type != 'none':
        S = Position(setup.lat_centre, setup.lon_centre, 45000)
        D = Position(setup.lat_centre, setup.lon_centre,     0)
        # try:
        #     graphWeather(setup.weather_name, setup.weather_type, [S.lat, S.lon, S.elev], [D.lat, D.lon, D.elev], 'spread_r', 100, spread_file=setup.spread_file, setup=setup)
        # except:
        #     pass
    arrTimes = np.array([])

    if len(stn_list) == 0:
        print('ERROR: No stations!')
        exit()


    # Init the wave\from picker
    arrTimes = calcAllTimes(stn_list, setup, sounding)


    # plt.tight_layout()

    # #plt.waitforbuttonpress(timeout=-1)

    # plt.show()

if __name__ == '__main__':

    pass