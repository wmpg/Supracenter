""" Tool for marking times of arrival of seismic/air waves in IRIS data. """

from __future__ import print_function, division, absolute_import

import os
import sys
import datetime
import argparse

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

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Fireballs.GetIRISData import readStationAndWaveformsListFile, butterworthBandpassFilter, \
    plotAllWaveforms, convolutionDifferenceFilter
from supra.Fireballs.SeismicTrajectory import latLon2Local, local2LatLon, timeOfArrival, waveReleasePoint, waveReleasePointWinds, \
    parseWeather, Constants
from supra.Fireballs.Program import configRead, configParse, position, station
from supra.Supracenter.angleConv import geo2Loc, loc2Geo
from supra.Supracenter.cyscan import cyscan
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.SPPT import perturb
from supra.Supracenter.weatherGraph import graphWeather
from wmpl.Utils.Earth import greatCircleDistance
from wmpl.Utils.PlotMap import GroundMap
from wmpl.Utils.Math import rotateVector
from wmpl.Utils.TrajConversions import datetime2JD

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
    n = 100000/traj[2]

    # B is the intersection between the trajectory vector and the ground
    A = n*traj

    # Convert back to geo coordinates
    B = np.array(loc2Geo(setup.lat_centre, setup.lon_centre, 0, B))
    A = np.array(loc2Geo(setup.lat_centre, setup.lon_centre, 0, A))

    print("Created Trajectory between A and B:")
    print("     A = {:10.4f}N {:10.4f}E {:10.2f}m".format(A[0], A[1], A[2]))
    print("     B = {:10.4f}N {:10.4f}E {:10.2f}m".format(B[0], B[1], B[2]))

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
        ref_pos = position(setup.lat_centre, setup.lon_centre, 0)
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

        d_time = 2*(setup.perturb_times*len(stn_list)*no_of_frags)
        count = 0

        #number of perturbations
        for ptb_n in range(setup.perturb_times):

            if ptb_n > 0:
                
                if setup.debug:
                    print("STATUS: Perturbation {:}".format(ptb_n))

                # generate a perturbed sounding profile
                sounding_p = perturb(setup, sounding, setup.perturb_method, \
                    sounding_u=sounding_u, sounding_l=sounding_l, \
                    spread_file=setup.perturbation_spread_file, lat=setup.lat_centre, lon=setup.lon_centre)
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
                                count += 1
                                sys.stdout.write("\rCalculating all times: {:5.2f} % ".format(count/d_time * 100))
                                sys.stdout.flush()
                                time.sleep(0.001)
                                #need filler values to make this a numpy array with fragmentation
                                if i == 0:

                                    setup.zangle = ZE
                                    setup.v = V
                                    A, B = trajCalc(setup)
                                    setup.traj_i = position(A[0], A[1], A[2])
                                    setup.traj_f = position(B[0], B[1], B[2])

                                    stn.position.pos_loc(ref_pos)
                                    setup.traj_f.pos_loc(ref_pos)
                                    # Station location in local coordinates
                                    # xx, yy, zz  = geo2Loc(lat0, lon0, elev0, \
                                    #                             stat_lat, stat_lon, stat_elev)

                                    #xx, yy, zz = rotateVector(np.array([xx, yy, zz]), np.array([0, 0, 1]), -np.pi/2)

                                    # Time to travel from trajectory to station
                                    b_time = timeOfArrival([stn.position.x, stn.position.y, stn.position.z], setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*V, \
                                                                np.radians(setup.azim), np.radians(ZE), setup, sounding=sounding_p, travel=False, fast=False, ref_loc=[ref_pos.lat, ref_pos.lon, ref_pos.elev])# + setup.t 
                                    S = waveReleasePointWinds([stn.position.x, stn.position.y, stn.position.z], setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*V, \
                                                                np.radians(setup.azim), np.radians(ZE), setup, sounding_p, [ref_pos.lat, ref_pos.lon, ref_pos.elev])
                                    # Distance from the source location to the station
                                    #b_dist = ((stn.position.x)**2 + (stn.position.y)**2)**0.5
                                    bDist[i] = ((S[0] - stn.position.x)**2 + (S[1] - stn.position.y)**2 + (S[2] - stn.position.z)**2)**0.5
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

        ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description="""
            ~~MakeIRISPicks~~ 
    Tool for marking times of arrival of 
    seismic/infrasound data in IRIS data. 

    Denis Vida
    """,
        formatter_class=argparse.RawTextHelpFormatter)

    arg_parser.add_argument('input_file', type=str, help='Path to Supracenter input file.')

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    #################

    setup = configRead(cml_args.input_file)
    configParse(setup, 'picks')

    ##########################################################################################################
    if not os.path.exists(setup.working_directory):
        os.makedirs(setup.working_directory)

    picks_name = setup.station_picks_file
    
    stn_list = []

    with open(picks_name) as f:
        for ii, line in enumerate(f):
            if ii > 0:
                line = line.split(',')
                stn = station(line[1], line[2], position(float(line[3]), float(line[4]), float(line[5])), '...', '...', '...')
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

    sounding = parseWeather(setup, consts)

    if setup.weather_type != 'none':
        S = position(setup.lat_centre, setup.lon_centre, 45000)
        D = position(setup.lat_centre, setup.lon_centre,     0)
        try:
            graphWeather(setup.weather_name, setup.weather_type, [S.lat, S.lon, S.elev], [D.lat, D.lon, D.elev], 'spread_r', 100, spread_file=setup.spread_file, setup=setup)
        except:
            pass
    arrTimes = np.array([])

    if len(stn_list) == 0:
        print('ERROR: No stations!')
        exit()


    # Init the wave\from picker
    arrTimes = calcAllTimes(stn_list, setup, sounding)


    # plt.tight_layout()

    # #plt.waitforbuttonpress(timeout=-1)

    # plt.show()