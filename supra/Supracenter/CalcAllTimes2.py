import numpy as np
import os
import copy
from functools import partial
import multiprocessing

from supra.Utils.Classes import Position
from supra.Fireballs.SeismicTrajectory import parseWeather, timeOfArrival
from supra.Supracenter.SPPT import perturb
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.cyscan2 import cyscan
#from supra.Supracenter.faultscan import cyscan as faultscan

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})
from supra.Supracenter.cyzInteg import zInterp

def findPoints(setup):    
    GRID_SPACE = 250
    MIN_HEIGHT = 10
    MAX_HEIGHT = 50000

    u = setup.trajectory.vector.xyz
    ground_point = np.array([0, 0, 0])
    # find top boundary of line given maximum elevation of trajectory
    if setup.trajectory.pos_i.elev != None:
        scale = -setup.trajectory.pos_i.elev/u[2]

    else:
        scale = -100000/u[2]

    # define line top boundary
    top_point = ground_point - scale*u

    ds = scale / (GRID_SPACE)

    points = []

    for i in range(GRID_SPACE + 1):
        points.append(top_point + i*ds*u)

    points = np.array(points)

    offset = np.argmin(np.abs(points[:, 2] - MAX_HEIGHT))
    bottom_offset = np.argmin(np.abs(points[:, 2] - MIN_HEIGHT))

    points = np.array(points[offset:(bottom_offset+1)])

    return points

def loop(setup, station_list, sounding_p, no_of_frags, points, ref_pos, theo, n):

    # For ballistic arrivals
    stn = station_list[n]

    # convert station coordinates to local coordinates based on the ref_pos
    stn.position.pos_loc(ref_pos)

    if setup.show_ballistic_waveform:
        setup.pos_f.pos_loc(ref_pos)
        bTimes = [0]*no_of_frags
        for i in range(no_of_frags):
            # count += 1
            # loadingBar('Calculating all times:', count, d_time)
            # sys.stdout.write("\rCalculating all times: {:5.2f} % ".format(count/d_time * 100))
            # sys.stdout.flush()
            #need filler values to make this a numpy array with fragmentation
            if i == 0:

                
                # Time to travel from trajectory to station
                b_time = timeOfArrival(stn.position.xyz, setup.trajectory.pos_f.x, setup.trajectory.pos_f.y, setup.trajectory.t, setup.trajectory.v, \
                                            setup.trajectory.azimuth.rad, setup.trajectory.zenith.rad, setup, points, setup.trajectory.vector.xyz, sounding=sounding_p, \
                                            travel=False, fast=False, ref_loc=ref_pos, theo=theo, div=37)# + setup.t 
                bTimes[i] = b_time
            else:
                bTimes[i] = np.nan

    else:
        bTimes = [np.nan]*no_of_frags
    # Fragmentation Prediction
    f_time = np.array([0]*no_of_frags)
    # Arrival times to the station
    fTimes = [0]*no_of_frags
    
    # Azimuths of initial rays
    fAz =    [0]*no_of_frags
    
    # Takeoff angles of initial rays
    fTf =    [0]*no_of_frags

    # Error of ray trace
    fEr =    [0]*no_of_frags

    # If manual fragmentation search is on
    if setup.show_fragmentation_waveform:

        for i, frag in enumerate(setup.fragmentation_point):
            # count += 1
            # loadingBar('Calculating all times:', count, d_time)

            # location of supracenter
            supra = frag.position
            
            # convert to local coordinates based off of the ref_pos
            supra.pos_loc(ref_pos)

            # Cut down atmospheric profile to the correct heights, and interp
            zProfile, _ = getWeather(np.array([supra.lat, supra.lon, supra.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), setup.weather_type, \
                    [ref_pos.lat, ref_pos.lon, ref_pos.elev], copy.copy(sounding_p), convert=False)
            
            #zProfile = zInterp(stn.position.z, supra.z, zProfile, div=37)


            # Travel time of the fragmentation wave
            f_time, frag_azimuth, frag_takeoff, frag_err = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True, \
                n_theta=setup.n_theta, n_phi=setup.n_phi, h_tol=setup.h_tol, v_tol=setup.v_tol)
            #f_time, frag_azimuth, frag_takeoff, _ = psoRayTrace(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile)

            fTimes[i] = f_time + frag.time
            fAz[i]    = frag_azimuth
            fTf[i]    = frag_takeoff
            fEr[i]    = frag_err

    else:

        # Repack all arrays into allTimes array
        fTimes = [np.nan]*no_of_frags

    return ([np.array(bTimes), np.array(fTimes), np.array(fAz), np.array(fTf), np.array(fEr)])

def calcAllTimes(stn_list, setup, sounding, theo=False, file_name='all_pick_times'):
    """ 
    Method to calculate all arrival times and place them into an array, so that they are not recalculated every
    time a new waveform is opened
            
    Arguments:
    data_list: [ndarray] array containing all station names and locations
    setup: [object] ini file parameters
    sounding: [ndarray] atmospheric profile

    Returns:
    allTimes: [ndarray] a ndarray that stores the arrival times for all stations over every perturbation
    [perturbation, station, 0 - ballistic/ 1 - fragmentation/ 2 - Azimuths of Initial Rays, 3 - Takeoff angles of Initial Rays, frag number (0 for ballistic)]
    """


    # define line bottom boundary
    if setup.show_ballistic_waveform:
        try:
            points = findPoints(setup)
        except AttributeError:
            points = []
    else:
        points = []

    if setup.perturb_times == 0:
        setup.perturb_times = 1

    #All perturbation happens here
    allTimes = [0]*setup.perturb_times

    # Ballistic Prediction
    ref_pos = Position(setup.lat_centre, setup.lon_centre, 0)
    no_of_frags = len(setup.fragmentation_point)

    # array of frags and ballistic arrivals have to be the same size. So, minimum can be 1
    if no_of_frags == 0:
        no_of_frags = 1

    # Initialize variables
    b_time = 0
    if not setup.show_ballistic_waveform and not setup.show_fragmentation_waveform:
        x = np.empty([setup.perturb_times, len(stn_list), 4, no_of_frags])
        A = np.full_like(x, np.nan, dtype=np.double)
        return A
    # For temporal perturbations, fetch the soudning data for the hour before and after the event
    if setup.perturb_method == 'temporal':

        # sounding data one hour later
        sounding_u = parseWeather(setup, time= 1)

        # sounding data one hour earlier
        sounding_l = parseWeather(setup, time=-1)

    else:
        sounding_u = []
        sounding_l = []

    if setup.perturb_method == 'ensemble':
        ensemble_file = setup.perturbation_spread_file
    else:
        ensemble_file = ''

    # if setup.show_fragmentation_waveform and setup.show_ballistic_waveform:
    #     d_time = 2*(setup.perturb_times*len(stn_list)*no_of_frags)
    # else:
    #     d_time = (setup.perturb_times*len(stn_list)*no_of_frags)
    # count = 0
    
    pool = multiprocessing.Pool(multiprocessing.cpu_count())
    iterable = range(len(stn_list))
    
    #number of perturbations
    for ptb_n in range(setup.perturb_times):

        if ptb_n > 0:
            
            # if setup.debug:
            #     print("STATUS: Perturbation {:}".format(ptb_n))

            # generate a perturbed sounding profile
            sounding_p = perturb(setup, sounding, setup.perturb_method, \
                sounding_u=sounding_u, sounding_l=sounding_l, \
                spread_file=setup.perturbation_spread_file, lat=setup.lat_centre, \
                lon=setup.lon_centre, ensemble_file=ensemble_file, ensemble_no=ptb_n)
        else:

            # if not using perturbations on this current step, then return the original sounding profile
            sounding_p = sounding

        # Initialize station times array
        stnTimes = [0]*len(stn_list)

        func = partial(loop, setup, stn_list, sounding_p, no_of_frags, points, ref_pos, theo) 
        stnTimes = pool.map(func, iterable)

        allTimes[ptb_n] = np.array(stnTimes)

    pool.close()
    pool.join()

    allTimes = np.array(allTimes)

    # Save as .npy file to be reused in SeismicTrajectory and Supracenter
    np.save(os.path.join(setup.working_directory, setup.fireball_name, file_name), allTimes)
    #errorMessage("All Times File saved as {:}".format(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times.npy')), 0)

    return allTimes

if __name__ == '__main__':

    from supra.Utils.Classes import Position, Config, Angle, Trajectory
    import pickle
    from supra.Fireballs.GetIRISData import readStationAndWaveformsListFile
    from wmpl.Utils.Earth import greatCircleDistance

    v_l = [11, 16, 21, 26, 31]
    ze_l = [5, 25, 45, 65, 85]

    for v in v_l:
        for ze in ze_l:
            print('STATUS: v = {:} | ze = {:}'.format(v, ze))
            DATA_FILE = 'data.txt'
            setup = Config()
            with open('/home/luke/Desktop/pkl/Theoretical2.pkl', 'rb') as f:
                setup = pickle.load(f)
            setup.trajectory = Trajectory(0, v*1000, pos_i=Position(48.05977, 13.10846, 85920.0), azimuth=Angle(354.67), zenith=Angle(ze))

            sounding = parseWeather(setup)
            # Theoretical Fireballs
            dir_path = os.path.join(setup.working_directory, setup.fireball_name)
            # Load the station and waveform files list
            data_file_path = os.path.join(dir_path, DATA_FILE)


            stn_list = readStationAndWaveformsListFile(data_file_path, rm_stat=setup.rm_stat)
            
            # Filter out all stations for which the mseed file does not exist
            filtered_stn_list = []

            names = []
            lats = []
            lons = []

            stn_list = stn_list + setup.stations

            for stn in stn_list:
                
                mseed_file = stn.file_name
                mseed_file_path = os.path.join(dir_path, mseed_file)

                if os.path.isfile(mseed_file_path) and stn.code not in setup.rm_stat:
                    filtered_stn_list.append(stn)


            stn_list = filtered_stn_list
            source_dists = []
            for stn in stn_list:

                stat_name, stat_lat, stat_lon = stn.code, stn.position.lat, stn.position.lon

                names.append(stat_name)
                lats.append(stat_lat)
                lons.append(stat_lon)

                # Calculate the distance in kilometers
                dist = greatCircleDistance(np.radians(setup.lat_centre), np.radians(setup.lon_centre), \
                    np.radians(stat_lat), np.radians(stat_lon))

                source_dists.append(dist)

            # Get sorted arguments
            dist_sorted_args = np.argsort(source_dists)

            # Sort the stations by distance
            stn_list = [stn_list[i] for i in dist_sorted_args]

            calcAllTimes(stn_list, setup, sounding, theo=True, file_name='v={:}-ze={:}'.format(v, ze))