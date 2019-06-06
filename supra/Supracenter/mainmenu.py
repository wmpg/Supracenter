# Supracenter - for spherical sound sources
# last edited July 2018
# Denis Vida, Peter Brown, Luke McFadden
# original Supracenter by Wayne Edwards
# searching algorithm was changed from genetic search to particle swarm optimization\

import sys
import os
import datetime
import argparse

import numpy as np

from supra.Supracenter.fetchECMWF import fetchECMWF
from supra.Supracenter.fetchMERRA import fetchMERRA
from supra.Supracenter.netCDFconv import storeNetCDFECMWF, storeHDF, storeNetCDFUKMO, readCustAtm
from supra.Supracenter.angleConv import roundToNearest, geo2Loc, loc2Geo, angle2NDE
from supra.Supracenter.stationDat import convStationDat, readTimes
from supra.Supracenter.psoSearch import psoSearch
from supra.Supracenter.geneticSearch import gLocator
from supra.Supracenter.plot import scatterPlot, residPlot, outputText
from supra.Fireballs.SeismicTrajectory import Constants
from supra.Fireballs.Program import configRead, configParse, position

def main(setup):

    # Read the config file
    #setup = configRead(input_name)

    # DATA CHECKS
    if setup.weather_type not in ['none', 'custom', 'merra', 'ecmwf', 'ukmo']:
        print('VARIABLE ERROR: Incorrect weather_type, must be none, custom, merra, ecmwf, or ukmo!')
        sys.exit()

    if setup.n_theta <= 0 or setup.n_phi <= 0:
        print('VARIABLE ERROR: Incorrect n_theta or n_phi, must be integer greater than 0!')
        sys.exit()

    if setup.angle_precision <= 0:
        print('VARIABLE ERROR: Incorrect precision, must be float greater than 0!')
        sys.exit()

    if setup.angle_error_tol <= 0:
        print('VARIABLE ERROR: Incorrect tol, must be float greater than 0!')
        sys.exit()

    if setup.v_sound <= 0:
        print('VARIABLE ERROR: Incorrect speed_of_sound, must be a float greater than 0!')

    if setup.max_error <= 0:
        print('VARIABLE ERROR: Incorrect max_error, must be a float greater than 0!')

    setup.grid_size = float(setup.grid_size)
    if setup.grid_size <= 0 or setup.grid_size > 10:
        print('VARIABLE ERROR: Incorrect grid_size, must be float greater than 0, and must give enough spaces in the search area!')
        sys.exit()

    if setup.weather_type == 'merra' and setup.fireball_datetime.year < 1980 and setup.get_data == True:
        a = input('WARNING: MERRA-2 data unavailable earlier than 1980! Bypass? (y/n) ')
        if a.lower() != 'y': 
            sys.exit()

    if setup.weather_type == 'ecmwf' and setup.fireball_datetime.year < 2008 and setup.get_data == True:
        a = input('WARNING: ECMWF data unavailable earlier than 2008! Bypass? (y/n) ')
        if a.lower() != 'y': 
            sys.exit()

    if setup.weather_type == 'ukmo' and setup.fireball_datetime.year < 1991 and setup.get_data == True:
        a = input('WARNING: UKMO data unavailable earlier than 1991! Bypass? (y/n) ')
        if a.lower() != 'y': 
            sys.exit()

    if (setup.swarmsize < 20 or setup.maxiter < 10) and len(setup.fragmentation_point) == 0:
        a = input('WARNING: swarm_size or max_iter not large enough, solution may be inaccurate or break! Continue? (y/n) ')
        if a.lower() != 'y': 
            sys.exit()

    if setup.weight_distance_max < setup.weight_distance_min:
        print('VARIABLE ERROR: dmax must be greater than or equal to dmin')
        sys.exit()

    # convert string to boolean
    # setup.enable_winds = (setup.enable_winds.lower() == 'true')
    # setup.get_data = (setup.get_data.lower() == 'true')
    # setup.debug = (setup.debug.lower() == 'true')
    # setup.perturb = (setup.perturb.lower() == 'true')


    # location to output data
    setup.output_name = os.path.join(setup.working_directory, setup.fireball_name)
    
    # weather file_name
    setup.file_name = os.path.join(setup.working_directory, setup.fireball_name, setup.sounding_file)

    # station data file name
    setup.station_name = os.path.join(setup.working_directory, setup.fireball_name, setup.station_picks_file)

    # all picks file name
    setup.perturb_picks = os.path.join(setup.working_directory, setup.fireball_name, setup.arrival_times_file)

    # Create fireball folder
    if not os.path.exists(setup.output_name):
        os.makedirs(setup.output_name)

    # Check that weather paths exist
    if not os.path.exists(setup.file_name) and setup.weather_type != 'none':
        print("FILE ERROR: Weather data file does not exist!")
        sys.exit()

    if not os.path.exists(setup.station_name):
        print("FILE ERROR: Station data file does not exist!")
        sys.exit()
    
    # Convert to integer if float
    setup.n_theta = int(setup.n_theta)
    setup.n_phi = int(setup.n_phi)

    # Init the constants
    consts = Constants()

    # Convert user defined occurence time to seconds
    # time is not restricted 
    if setup.restricted_time == '':
        setup.restricted_time = None

    # time is restriced
    
    # key setting for if not using a manual search
    if setup.manual_fragmentation_search == '':
        setup.manual_fragmentation_search = []

    if setup.restricted_time == '' and setup.manual_fragmentation_search != '':
        print("ERROR: time must be restricted if doing a manual search!")
        exit()


    ### Get atmospheric data ###

    # Isotropic weather profile
    if setup.weather_type == 'none':
        dataset = np.array([[0, setup.v_sound, 0, 0], [1000, setup.v_sound, 0, 0]])

    # Custom weather data    
    elif setup.weather_type == 'custom':
        dataset = readCustAtm(setup.file_name, consts)

        # Check file name
        if len(dataset) < 2:
            print('FILE ERROR: file must have at least two data points!')
            sys.exit()

        # Check file type
        if '.txt' not in setup.file_name:
            print("FILE ERROR: custom data set must be a .txt file!")
            sys.exit()

    # MERRA-2
    elif setup.weather_type == 'merra':

        # Check file type
        if '.nc' not in setup.file_name:
            print("FILE ERROR: custom data set must be a .nc file!")
            sys.exit()

        # Run data fetching script
        # if setup.get_data == True:
        #     print('Getting data...')
        #     fetchMERRA(setup)
        dataset = storeHDF(setup.file_name, consts)

    # ECMWF
    elif setup.weather_type == 'ecmwf':

        # Check file type
        if '.nc' not in setup.file_name:
            print("FILE ERROR: custom data set must be a .nc file from ECMWF!")
            sys.exit()

        lat_centre = (setup.search_area[0] + setup.search_area[1])/2
        lon_centre = (setup.search_area[2] + setup.search_area[3])/2
        # Run data fetching script
        # if setup.get_data == True:
        #     fetchECMWF(setup, setup.file_name)
        dataset = storeNetCDFECMWF(setup.file_name, lat_centre, lon_centre, consts)

    # UKMO
    elif setup.weather_type == 'ukmo':

        # Check file type
        if '.nc' not in setup.file_name:
            print("FILE ERROR: custom data set must be a .nc file from UKMO!")
            sys.exit()

        dataset = storeNetCDFUKMO(setup.file_name, setup.search_area, consts)

    print("Status: Atmospheric Data Loaded")
    ##########################
    
    # Grab station info from file
    s_info, s_name, weights, ref_pos = convStationDat(setup.station_name, setup=setup, \
        d_min=setup.weight_distance_min, d_max=setup.weight_distance_max)
    station_no = s_info[:, 5]

    ref_pos = position(ref_pos[0], ref_pos[1], ref_pos[2])

    # Set up perturbed array
    if setup.perturb == True:
        allTimes = readTimes(setup.perturb_picks)
        print("Status: Loaded picks from perturbations")

    search_min = position(setup.search_area[0], setup.search_area[2], setup.search_height[0])
    search_max = position(setup.search_area[1], setup.search_area[3], setup.search_height[1])
    # if setup.restrict_to_trajectory:
    #     setup.search_area = min(setup.lat_i, setup.lat_f), max(setup.lat_i, setup.lat_f), min(setup.lon_i, setup.lon_f), max(setup.lon_i, setup.lon_f)

    # Convert search bounds to local coordinates
    search_min.pos_loc(ref_pos)
    search_max.pos_loc(ref_pos)

    # Search boundary
    search_area = [search_min.lat, search_max.lat, search_min.lon, search_max.lon, search_min.elev, search_max.elev]

    setup.search = search_area
    setup.ref_pos = ref_pos
    
    if setup.perturb == True:
        run_times = len(allTimes)
    else:
        run_times = 0
    
    results_arr = [0]*(run_times + 1)

    for i in range(run_times + 1):
        # pass all options into main program
        # particle swarm optimization search
        if i != 0:
            if setup.perturb:
                for j in range(len(station_no)):
                    s_info[j, 4] = allTimes[i-1, int(station_no[j]), 1, setup.observe_frag_no-1]

        #if setup.search_type == 'pso':
        results_arr[i] = psoSearch(s_info, weights, s_name, setup, dataset, consts)

        # # grid search
        # elif setup.search_type == 'gs':
        #     gLocator(search_area, s_info, weights, ref_pos, s_name, setup.reported_points, setup.start_datetime, \
        #              setup.rest_time, dataset, setup.output_name, setup.single_point, consts)

        # error in name
        # else:
        #     print('INI ERROR: Incorrect search type, use "pso" for particle swarm optimization or "gs" for grid search. (gs currently not supported)')
        #     print("ERROR: grid search (gs) mode currently not supported. Please use Particle Swarm Optimization (pso)")
        #     sys.exit()
        
    n_stations = len(s_info)
    xstn = s_info[0:n_stations, 0:3]
    tstn = s_info[0:n_stations, 4]

    # scatter plot(s)
    min_search, max_search = scatterPlot(setup, results_arr, n_stations, xstn, s_name, dataset)

    # residual plot
    residPlot(results_arr, s_name, xstn, setup.output_name, n_stations)

    # output results
    outputText(min_search, max_search, setup, results_arr, n_stations, s_name, xstn, tstn, weights)

    print("Status: Supracenter Complete!")

if __name__ == "__main__":

    ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description="""
                ~~Supracenter~~ 
    Find the point of fragmentation of a fireball 
    with atmospheric data and seismic data.

    Denis Vida, Peter Brown, Luke McFadden
    Original Supracenter by Wayne Edwards
        """,
        formatter_class=argparse.RawTextHelpFormatter)

    arg_parser.add_argument('input_file', type=str, help='Path to Supracenter input file.')

    # arg_parser.add_argument('-t', '--maxtoffset', metavar='MAX_TOFFSET', nargs='1', \
    #     help='Maximum time offset between the stations.', type=float, default=1.0)
    
    # arg_parser.add_argument('-g', '--disablegravity', \
    #     help='Disable gravity compensation.', action="store_true")

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    ############################

    setup = configRead(cml_args.input_file)
    configParse(setup, 'supracenter')
    main(setup)
