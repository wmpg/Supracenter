import numpy as np
import os
import copy

from supra.Utils.Classes import Position
from supra.GUI.GUITools import errorMessage
from supra.Fireballs.SeismicTrajectory import parseWeather, timeOfArrival
from supra.Supracenter.SPPT import perturb
from supra.Utils.Formatting import loadingBar
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.cyscan import cyscan

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

    if setup.show_fragmentation_waveform and setup.show_ballistic_waveform:
        d_time = 2*(setup.perturb_times*len(stn_list)*no_of_frags)
    else:
        d_time = (setup.perturb_times*len(stn_list)*no_of_frags)
    count = 0

    #number of perturbations
    for ptb_n in range(setup.perturb_times):

        if ptb_n > 0:
            
            if setup.debug:
                print("STATUS: Perturbation {:}".format(ptb_n))

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

        #number of stations
        for n, stn in enumerate(stn_list):

            # For ballistic arrivals
            if setup.show_ballistic_waveform:
                bTimes = [0]*no_of_frags
                for i in range(no_of_frags):
                    count += 1
                    loadingBar('Calculating all times:', count, d_time)
                    # sys.stdout.write("\rCalculating all times: {:5.2f} % ".format(count/d_time * 100))
                    # sys.stdout.flush()
                    #need filler values to make this a numpy array with fragmentation
                    if i == 0:

                        stn.position.pos_loc(ref_pos)
                        setup.traj_f.pos_loc(ref_pos)

                        # Time to travel from trajectory to station
                        b_time = timeOfArrival(stn.position.xyz, setup.trajectory.pos_f.x, setup.trajectory.pos_f.y, setup.trajectory.t, setup.trajectory.v, \
                                                    setup.trajectory.azimuth.rad, setup.trajectory.zenith.rad, setup, sounding=sounding_p, \
                                                    travel=False, fast=False, ref_loc=ref_pos)# + setup.t 

                        bTimes[i] = b_time
                    else:
                        bTimes[i] = np.nan

            else:
                bTimes = [np.nan]*no_of_frags

            # Fragmentation Prediction
            f_time = np.array([0]*no_of_frags)

            # If manual fragmentation search is on
            if setup.show_fragmentation_waveform:
                fTimes = [0]*no_of_frags
                for i, frag in enumerate(setup.fragmentation_point):
                    count += 1
                    loadingBar('Calculating all times:', count, d_time)

                    # location of supracenter
                    supra = frag.position
                    
                    # convert to local coordinates based off of the ref_pos
                    supra.pos_loc(ref_pos)

                    # convert station coordinates to local coordinates based on the ref_pos
                    stn.position.pos_loc(ref_pos)

                    # Cut down atmospheric profile to the correct heights, and interp
                    zProfile, _ = getWeather(np.array([supra.lat, supra.lon, supra.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), setup.weather_type, \
                            [ref_pos.lat, ref_pos.lon, ref_pos.elev], copy.copy(sounding_p), convert=False)

                    # Travel time of the fragmentation wave
                    f_time, _, _ = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True, \
                        n_theta=setup.n_theta, n_phi=setup.n_phi, precision=setup.angle_precision, tol=setup.angle_error_tol)

                    fTimes[i] = f_time + frag.time

            else:

                # Repack all arrays into allTimes array
                fTimes = [np.nan]*no_of_frags

            stnTimes[n] = ([np.array(bTimes), np.array(fTimes)])

        allTimes[ptb_n] = np.array(stnTimes)

    allTimes = np.array(allTimes)

    # Save as .npy file to be reused in SeismicTrajectory and Supracenter
    np.save(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times'), allTimes)
    errorMessage("All Times File saved as {:}".format(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times.npy')), 0)

    return allTimes