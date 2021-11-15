"""Finds the optimal Supracenter using particle swarm optimization and graphs and exports the data"""

import copy
import multiprocessing
import datetime

import numpy as np
from supra.Utils.pso import pso
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Supracenter.cyscan5 import cyscan
from supra.Utils.Classes import Position
from supra.Supracenter.plot import outputWeather
from supra.GUI.Tools.GUITools import *
from supra.Utils.Formatting import *

def timeFunction(x, *args):
    ''' Helper function for PSO
    Takes in supracenter ranges, and calculates the travel time, which is used to find the error.
    Returns the residual error with each guess of a supracenter from pso()

    Arguments:
        x: [ndarray] position to search with
        args: list of passed arguments
            stns: [list] list of station positions and arrival times
            w: [list] list of station weights
            kotc: [float] user-defined occurence time
            tweaks: [Object] user-defined option
            ref_pos: [list] mean station location used for converting to/from local coordinates
            dataset: [ndarray] atmospheric profile for the entire search area
            pool: [multiprocessing pool] pool of workers for multiprocessing

    Returns:
        err: [float] error in the current position, x, searched
    '''

    # Retrieve passed arguments
    stns, w, kotc, setup, ref_pos, atmos, prefs, v, pert_num = args

    # number of stations total
    n_stations = len(stns)

    # Residuals to each station
    sotc = np.empty(n_stations)

    # Initialize variables
    # Weight of each station
    wn = w

    # total weight
    nwn = sum(w)

    # Mean occurrence time
    motc = 0

    S = Position(0, 0, 0)
    S.x, S.y, S.z = x[0], x[1], x[2]
    S.pos_geo(ref_pos)

    # number of stations total
    n_stations = len(stns)

    # Station Times
    tobs = stns[0:n_stations, 4]

    # Station Location
    xstn = stns[0:n_stations, 0:3]

    # Initialize arrays
    # Simulated travel times to each station
    time3D = np.empty(n_stations)

    # Residuals to each station
    sotc = np.empty(n_stations)
    
    for j in range(n_stations):
        # if station has weight
        if w[j] > 0:
            
            D = Position(0, 0, 0)
            D.x, D.y, D.z = xstn[j, 0], xstn[j, 1], xstn[j, 2]
            D.pos_geo(ref_pos)

            if pert_num == 0:

                # No perturbations used here
                sounding, _ = atmos.getSounding(lat=[S.lat, D.lat], lon=[S.lon, D.lon], heights=[S.elev, D.elev])
           
            else:

                # Sounding is a specfic perturbation number
                # TODO: Go back and fix how this is done, not all perts need to be generated, just one here 
                nom, sounding = atmos.getSounding(lat=[S.lat, D.lat], lon=[S.lon, D.lon], heights=[S.elev, D.elev])
                
                # sounding is none when there is an error in getting sounding
                if sounding is not None:
                
                    sounding = sounding[pert_num - 1]
                
                else:
                
                    sounding = nom

            # Use distance and atmospheric data to find path time
            time3D[j], _, _, _ = cyscan(S.xyz, D.xyz, sounding, \
                        wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,\
                        h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)
            # Residual time for each station
            sotc[j] = tobs[j] - time3D[j]

        # If station has no weight
        else:
            sotc[j] = tobs[j]
    
    motc = np.mean(sotc)
    ##########
    N_s = np.count_nonzero(~np.isnan(sotc))
    # User defined occurrence time
    if kotc != None:

        # make kOTC a list
        err = np.dot(wn, np.absolute(sotc - np.array([kotc]*n_stations)))/nwn

    # Unknown occurrence time
    else:
        #err = np.dot(wn, np.absolute(sotc - motc))/nwn
        err = 0
        for s in sotc:
            err += (1 + (s - motc)**2)**0.5 - 1

        err = err/N_s

    # if setup.debug:
    #     # print out current search location
    #     print("Supracenter: {:10.2f} m x {:10.2f} m y {:10.2f} m z  Time: {:8.2f} Error: {:25.2f}".format(x[0], x[1], x[2], motc, err))

    # variable to be minimized by the particle swarm optimization
    
    perc_fail = 100 - (n_stations-N_s)/n_stations*100

    # temporary adjustment to try and get the most stations
    if N_s >= 3:
        total_error = err# + 2*max(error_list)*(failed_stats)
    else:
        total_error = np.inf

    if np.isnan(total_error):
        total_error = np.inf

    if prefs.debug:
        print("Error {:10.4f} | Solution {:10.4f}N {:10.4f}E {:8.2f} km {:8.2f} s | Failed Stats {:3} {:}".format(total_error, S.lat, S.lon, S.elev/1000, motc, n_stations-N_s, printPercent(perc_fail, N_s)))
        # Quick adjustment to try and better include stations


    return total_error

def timeConstraints(x, *args):
    
    # Retrieve passed arguments
    stns, w, kotc, setup, ref_pos, dataset, v = args

    # number of stations total
    n_stations = len(stns)

    # Residuals to each station
    sotc = np.empty(n_stations)

    # Initialize variables
    # Weight of each station
    wn = w

    # Mean occurrence time
    motc = 0

    ### multiprocessing

    # number of stations total
    n_stations = len(stns)

    # Station Times
    tobs = stns[0:n_stations, 4]

    # Station Location
    xstn = stns[0:n_stations, 0:3]

    # Initialize arrays
    # Simulated travel times to each station
    time3D = np.empty(n_stations)

    # Residuals to each station
    sotc = np.empty(n_stations)
    
    for j in range(n_stations):
        # if station has weight
        if w[j] > 0:

            # Create interpolated atmospheric profile for use with cyscan
            sounding, points = getWeather([x[0], x[1], x[2]], xstn[j, :], setup.weather_type, ref_pos, copy.copy(dataset))


            # Use distance and atmospheric data to find path time
            time3D[j], _, _, _ = cyscan(np.array([x[0], x[1], x[2]]), np.array(xstn[j, :]), sounding, \
                                                wind=setup.enable_winds, n_theta=setup.n_theta, n_phi=setup.n_phi, h_tol=setup.h_tol, v_tol=setup.v_tol)
            # Residual time for each station
            sotc[j] = tobs[j] - time3D[j]

        # If station has no weight
        else:
            sotc[j] = tobs[j]
    
    t_avg = (setup.t_min + setup.t_max)/2


    diff = abs(motc - t_avg)
    TOL = setup.t_max - t_avg
    motc = np.dot(wn, sotc)/sum(wn)

    
    if diff < TOL:
        return diff
    else:
        return -diff

def trajConstraints(x, *args):
    # Retrieve passed arguments
    stns, w, kotc, setup, ref_pos, dataset, v = args

    TOL = 1000

    diff_vect = np.array(x) - x[2]/v[2] * v

    diff = (diff_vect[0]**2 + diff_vect[1]**2 + diff_vect[2]**2)**0.5


    # print(diff)
    if diff <= TOL:
        return diff
    else:
        return -diff

def psoSearch(stns, w, s_name, bam, prefs, ref_pos, manual=False, pert_num=0, override_supra=[]):
    """ Optimizes the paths between the detector stations and a supracenter to find the best fit for the 
        position of the supracenter, within the given search area. The supracenter is found with an initial guess,
        in a given grid, and is moved closer to points of better fit through particle swarm optimization.
        Produces a graph with the stations and residual results, and prints out the optimal supracenter location
    """

    print('Data converted. Searching...')

    setup = bam.setup
    atmos = bam.atmos

    search_min = Position(setup.lat_min, setup.lon_min, setup.elev_min)
    search_max = Position(setup.lat_max, setup.lon_max, setup.elev_max)

    if not manual:
        try:
            search_min.pos_loc(ref_pos)
            search_max.pos_loc(ref_pos)
        except (AttributeError, TypeError) as e:
            errorMessage('Search min and search max have not been defined! Aborting search!', 2, info='Please define a search area in the "Sources" tab on the left side of the screen!', detail='{:}'.format(e))
            return None


    output_name = prefs.workdir
    
    if not isinstance(override_supra, list):
        single_point = override_supra
    else:    
        single_point = setup.manual_fragmentation_search[0]

    if single_point.toList()[0] is None and manual:
        errorMessage('Manual Fragmentation Point undefined', 2, info='Unable to parse: Lat: {:} Lon: {:} Elev: {:} Time: {:}'.format(*single_point.toList()))
        return None        


    ref_time = setup.fireball_datetime

    if setup.enable_restricted_time:
        kotc = setup.restricted_time
    else:
        kotc = None

    # check if user defined occurrence time is used
    if kotc != None:
    
        kotc = (kotc - ref_time).total_seconds()

    # number of stations total
    n_stations = len(stns)

    # Station Location
    xstn = stns[0:n_stations, 0:3]

    # Initialize arrays
    # Travel time to each station
    time3D = np.zeros(n_stations)

    # Initial azimuths angles of each station
    az = np.zeros(n_stations)

    # Initial takeoff angles of each station
    tf = np.zeros(n_stations)

    # difference in theoretical and simulated travel times
    sotc = np.zeros_like(n_stations)

    # Initialize variables
    # combined weights
    nwn = sum(w)

    if prefs.ballistic_en:
        try:
            v = -setup.trajectory.vector.xyz
            setup.ref_pos = setup.trajectory.pos_f
            if prefs.debug:
                print("Constraining Trajectory")
        except:
            v = [None]
            if prefs.debug:
                print("Free Search")
    else:
        v = [None]
        if prefs.debug:
            print("Free Search")


    # If automatic search
    if not manual:
        
        # Prevent search below stations
        if search_min.elev < max(xstn[:, 2]):

            # Must be just above the stations
            search_min.elev = max(xstn[:, 2]) + 0.0001

        # Boundaries of search volume
        #  [x, y, z] local coordinates

        # arguments to be passed to timeFunction()
        args = (stns, w, kotc, setup, ref_pos, atmos, prefs, v, pert_num)

        # Particle Swarm Optimization
        # x_opt - optimal supracenter location
        # f_opt - optimal supracenter error

        # if setup.restrict_to_trajectory:
        #     #cons = [lineConstraintx, lineConstrainty, lineConstraintz]
        #     x_opt, f_opt = pso(timeFunction, lb, ub, f_ieqcons=lineConstraint, args=args, swarmsize=int(setup.swarmsize), maxiter=int(setup.maxiter), \
        #                 phip=setup.phip, phig=setup.phig, debug=False, omega=setup.omega, minfunc=setup.minfunc, minstep=setup.minstep) 
        # else:

        # # Restricted to trajectory
        # if v[0] != None:
        #     x_opt_temp, f_opt, sup, errors = pso(timeFunction, search_min.xyz, search_max.xyz, \
        #             ieqcons=[trajConstraints, timeConstraints], args=args,\
        #             swarmsize=int(prefs.pso_swarm_size), maxiter=int(prefs.pso_max_iter), \
        #             phip=prefs.pso_phi_p, phig=prefs.pso_phi_g, debug=False, omega=prefs.pso_omega, \
        #             minfunc=prefs.pso_min_error, minstep=prefs.pso_min_step, processes=1, particle_output=True) 
        
        # Unrestricted

        x_opt_temp, f_opt, sup, errors = pso(timeFunction, search_min.xyz, search_max.xyz, \
                args=args, swarmsize=int(prefs.pso_swarm_size), maxiter=int(prefs.pso_max_iter), \
                phip=prefs.pso_phi_p, phig=prefs.pso_phi_g, debug=False, omega=prefs.pso_omega,\
                minfunc=prefs.pso_min_error, minstep=prefs.pso_min_step, processes=1, particle_output=True) 


        print('Done Searching')
        
        x_opt = Position(0, 0, 0)
        x_opt.x = x_opt_temp[0]
        x_opt.y = x_opt_temp[1]
        x_opt.z = x_opt_temp[2]
        x_opt.pos_geo(ref_pos)

    # If manual search
    else:
        single_point.position.pos_loc(ref_pos)
        x_opt = single_point.position
        sup = single_point.position.xyz
        errors=0


    # Get results for current Supracenter
    time3D, az, tf, r, motc, sotc, trace = outputWeather(n_stations, x_opt, stns, setup, \
                                                ref_pos, atmos, output_name, s_name, kotc, w, prefs)

    for ii, element in enumerate(time3D):
        if np.isnan(element):
            w[ii] = 0
            sotc[ii] = 0

    # Find error for manual searches
    if manual:

        f_opt = np.dot(w, np.absolute(sotc - np.dot(w, sotc)/nwn))/nwn

    # x, y distance from Supracenter to each station
    horz_dist = np.zeros(n_stations)
    for i in range(n_stations):

        horz_dist[i] = np.sqrt((x_opt.x - xstn[i, 0])**2 + (x_opt.y - xstn[i, 1])**2)/1000


    # Calculate and Set the Occurrence Time into HH:MM:SS
    time_diff = motc + ref_time.microsecond/1e6 + ref_time.second + ref_time.minute*60 + ref_time.hour*3600

    try:
        otc = (datetime.datetime.min + datetime.timedelta(seconds=time_diff)).time()
    except (ValueError, OverflowError):
        print('STATUS: Unable to parse otc')
        otc = None

    try:
        #while max(errors) > setup.max_error/100:
        a = []
        std_error = np.std(errors)
        lim = np.mean(errors) + 0*std_error

        for i in range(len(errors)):
            if errors[i] >= lim:
                a.append(i)

        errors = np.delete(errors, (a), axis=0)
        sup = np.delete(sup, (a), axis=0)
    except:
        print("WARNING: Unable to filter errors")

    class Results:

        def __init__(self):
            pass

    results = Results()

    results.r = r
    results.w = w
    results.x_opt = x_opt
    results.f_opt = f_opt
    results.sup = sup
    results.errors = errors
    results.horz_dist = horz_dist
    results.time3D = time3D
    results.az = az
    results.tf = tf
    results.motc = motc
    results.sotc = sotc
    results.kotc = kotc
    results.otc = otc
    results.trace = trace

    return results

    # # scatter plot(s)
    # min_search, max_search = scatterPlot(single_point, n_stations, xstn, s_name, r, x_opt, \
    #                                         reported_points, search, output_name, ref_pos, sup, errors, tweaks, dataset)

    # # residual plot
    # residPlot(x_opt, s_name, xstn, r, output_name, n_stations)

    # # output results
    # outputText(min_search, max_search, single_point, ref_time, otc, kotc, x_opt, f_opt, n_stations, tweaks, s_name, xstn, \
    #                                                                     r, w, az, tf, time3D, horz_dist, output_name, tstn)