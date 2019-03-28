"""Finds the optimal Supracenter using particle swarm optimization and graphs and exports the data"""

import os
import copy
import time
import multiprocessing
import sys
import datetime

import numpy as np
from pyswarm import pso
from functools import partial
import matplotlib.pyplot as plt
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.cyscan import cyscan
from supra.Supracenter.angleConv import loc2Geo, geo2Loc, trajRestriction, point2LineDist3D, angle2NDE
from supra.Supracenter.plot import residPlot, scatterPlot, outputText, outputWeather

sup = [0, 0, 0]
errors = [0] 

def loop(x, stns, w, tweaks, ref_pos, dataset, j):
    """ Helper function that is the loop to use with multiprocessing
    returns the travel time to each station as calculated with cyscan

    Arguments:
        h: passed variables and iterable variable
            j: [int] iterable variable
            x: [list] current Supracenter position
            stns: [list] station location and time of arrivals
            w: [list] weights of each station
            tweaks: [Object] user-defined options
            ref_pos: [list] mean position of stations for use when converting to/from local coordinates
            dataset: [ndarray] atmospheric profile of the search area

    Returns:
        sotc: calculated time for the current search position signal to hit each station
    """

    # number of stations total
    n_stations = len(stns)

    # Station Times
    tobs = stns[0:n_stations, 4]

    # Station Location
    xstn = stns[0:n_stations, 0:3]

    # Initialize arrays
    # Simulated travel times to each station
    time3D = np.zeros(n_stations)

    # Residuals to each station
    sotc = np.zeros(n_stations)

    # 0 - custom weather, 1 - MERRA, 2 - ECMWF, 3 - UKMO
    weather_type = tweaks[0]
    
    # if station has weight
    if w[j] > 0:

        # Create interpolated atmospheric profile for use with cyscan
        sounding, points = getWeather([x[0],x[1],x[2]], xstn[j, :], weather_type, ref_pos, copy.copy(dataset))

        # Rotate winds to match with coordinate system
        #sounding[:, 3] = angle2NDE(sounding[:, 3])

        # Use distance and atmospheric data to find path time
        time3D[j], _, _ = cyscan(np.array([x[0],x[1],x[2]]), np.array(xstn[j, :]), sounding, \
                                            wind=tweaks[1], n_theta=tweaks[2], n_phi=tweaks[3], \
                                            precision=tweaks[4], tol=tweaks[5])

        # Residual time for each station
        sotc[j] = tobs[j] - time3D[j]

    # If station has no weight
    else:
        sotc[j] = tobs[j]

    return sotc[j]

# def lineConstraintx(x, *args):
#     _, _, _, setup, _, _, _ = args

#     # convert angles to radians
#     ze = np.radians(setup.zangle)
#     az = np.radians(setup.azim)
#     t = np.array([np.sin(az)*np.cos(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])
#     #t = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

#     A = geo2Loc(setup.ref_pos[0], setup.ref_pos[1], setup.ref_pos[2], setup.lat_i, setup.lon_i, setup.elev_i)
#     A = np.array(A)
#     A[2] *= 1000
#     s = (x[0] - A[0]) / t[0]

#     result = s*t[0] + A[0] - x[0]

#     # pso will allow the point if the constraint is >= 0.0
#     if abs(result) <= setup.traj_tol:
#         return 1
#     else:
#         return -1   

# def lineConstrainty(x, *args):
#     _, _, _, setup, _, _, _ = args

#     # convert angles to radians
#     ze = np.radians(setup.zangle)
#     az = np.radians(setup.azim)
#     t = np.array([np.sin(az)*np.cos(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])
#     #t = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

#     A = geo2Loc(setup.ref_pos[0], setup.ref_pos[1], setup.ref_pos[2], setup.lat_i, setup.lon_i, setup.elev_i)
#     A = np.array(A)
#     A[2] *= 1000
#     s = (x[1] - A[1]) / t[1]

#     result = s*t[1] + A[1] - x[1]
#     # pso will allow the point if the constraint is >= 0.0
#     if abs(result) <= setup.traj_tol:
#         return 1
#     else:
#         return -1   


# def lineConstraintz(x, *args):
#     _, _, _, setup, _, _, _ = args

#     # convert angles to radians
#     ze = np.radians(setup.zangle)
#     az = np.radians(setup.azim)

#     t = np.array([np.sin(az)*np.cos(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])
#     A = geo2Loc(setup.ref_pos[0], setup.ref_pos[1], setup.ref_pos[2], setup.lat_i, setup.lon_i, setup.elev_i)
#     A = np.array(A)
#     A[2] *= 1000
#     s = (x[2] - A[2]) / t[2]

#     result = s*t[2] + A[2] - x[2]

#     # pso will allow the point if the constraint is >= 0.0
#     if abs(result) <= setup.traj_tol:
#         return 1
#     else:
#         return -1  

def lineConstraint(x, *args):

    _, _, _, setup, ref_pos, _, _ = args

    # convert angles to radians
    ze = np.radians(setup.zangle)
    az = np.radians(setup.azim)

    t = np.array([np.cos(az)*np.sin(ze), np.sin(az)*np.sin(ze), -np.cos(ze)])
    A = geo2Loc(ref_pos[1], ref_pos[0], ref_pos[2], setup.lat_i, setup.lon_i, setup.elev_i)
    A = np.array(A)
    A[2] *= 1000
    s = (x[2] - A[2]) / t[2]

    result = s*np.sum(t) + np.sum(A) - np.sum(x)

    if abs(result) > setup.traj_tol:
        result = -1
    else:
        result = 0
        
    # pso will allow the point if the constraint is >= 0.0
    return result

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
    stns, w, kotc, setup, ref_pos, dataset, pool = args

    # number of stations total
    n_stations = len(stns)

    # Initialize arrays
    # Simulated travel times to each station
    time3D = np.zeros(n_stations)

    # Initial azimuthal angles to each station
    az = np.zeros(n_stations)

    # Initial takeoff angles to each station
    tf = np.zeros(n_stations)

    # Residuals to each station
    sotc = np.zeros(n_stations)

    # Initialize variables
    # Weight of each station
    wn = w

    # total weight
    nwn =sum(w)

    # Mean occurrence time
    motc = 0

    ### multiprocessing

    # Pass data through to multiprocess loop
    iterable = range(n_stations)

    # Store functions to be used with multiprocessing
    func = partial(loop, x, stns, w, [setup.weather_type, setup.enable_winds, setup.n_theta, setup.n_phi,\
                                        setup.angle_precision, setup.angle_error_tol], ref_pos, dataset)
    
    # Compute time of flight residuals for all stations
    sotc = pool.map(func, iterable)

    if setup.max_time != 0:
        motc = np.dot(wn, sotc)/sum(wn)
    else:
        motc = 1
    ##########

    # User defined occurrence time
    if kotc != None:

        # make kOTC a list
        err = np.dot(wn, np.absolute(sotc - np.array([kotc]*n_stations))**setup.fit_type)/nwn

    # Unknown occurrence time
    else:
        err = np.dot(wn, np.absolute(sotc - np.dot(wn, sotc)/nwn)**setup.fit_type)/nwn

    global sup 
    global errors

    # if setup.restrict_to_trajectory:
    #     a, b = trajRestriction(setup)
    #     c = np.array([x[0], x[1], x[2]])
    #     d = point2LineDist3D(a, b, c)

    # Save points and errors for plotting

    # check: Err == nan
    if not err > 0:
        err = setup.max_error

    # Method to get times to be inside the restrictions. Set errors of time outside of range to be proportional
    # to how far outside of the range it is

    if motc < setup.min_time or motc > setup.max_time:
        err = min(abs(motc - setup.min_time), abs(motc - setup.max_time))*setup.max_error + setup.max_error

    # if setup.restrict_to_trajectory:
    #     if d > setup.traj_tol:
    #         err = abs(setup.traj_tol - d)*min(abs(motc - setup.min_time), abs(motc - setup.max_time))*setup.max_error + setup.max_error

    # check large error
    if err > setup.max_error:

        errors = np.hstack((errors, setup.max_error))
        sup = np.vstack((sup, [x[0], x[1], x[2]]))

    else:
        errors = np.hstack((errors, err))
        sup = np.vstack((sup, [x[0], x[1], x[2]]))

    if setup.debug:
        # print out current search location
        print("Supracenter: {:10.2f} m x {:10.2f} m y {:10.2f} m z  Time: {:8.2f} Error: {:25.2f}".format(x[0], x[1], x[2], motc, err))

    # variable to be minimized by the particle swarm optimization
    return err

def psoSearch(stns, w, s_name, setup, dataset, consts):
    """ Optimizes the paths between the detector stations and a supracenter to find the best fit for the 
        position of the supracenter, within the given search area. The supracenter is found with an initial guess,
        in a given grid, and is moved closer to points of better fit through particle swarm optimization.
        Produces a graph with the stations and residual results, and prints out the optimal supracenter location

    Arguments:
        search: [list] [x_min, x_max, y_min, y_max, z_min, z_max] limits of coordinates to search within
        stns: [array] [Latitude, Longitude, Elevation] for each station (converted into locatl coordinates)
        w: [list] list of weights for each station (0 < w < 1, 1 - default, 0 - ignore station)
        ref_pos: [list] [latitude, longitude, elevation] reference coordinate to be used for the local coordinate system
                    by default, is the average latitude and longitude between the stations, and an elevation of 0.
        s_name: [lst] names of stations for plotting
        setup: [object] object storing ini inputs
        tweaks: [object]  parameters that can be changed from the main menu
        dataset: [array] weather profile for the entire search area, for use with weatherInterp
        consts: [object] Physical constants used throught the program
    """

    print('Data converted. Searching...')

    search = setup.search
    ref_pos = setup.ref_pos
    output_name = setup.output_name
    single_point = setup.manual_fragmentation_search
    ref_time = setup.start_datetime
    kotc = setup.restricted_time

    if len(single_point) != 0:
        print('Position from mean station location')

    # check if user defined occurrence time is used
    if kotc != None:
        kotc -= ref_time.hour*3600 + ref_time.minute*60 + ref_time.second + ref_time.microsecond/1e6

    # number of stations total
    n_stations = len(stns)

    # Station Location
    xstn = stns[0:n_stations, 0:3]
    tstn = stns[0:n_stations, 4]

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

    # Atmospheric resolution
    grid_size = setup.grid_size

    # Prevent search below stations
    if search[4] < max(xstn[:, 2]):

        # Must be just above the stations
        search[4] = max(xstn[:, 2]) + 0.0001

    # If automatic search
    if len(single_point) == 0:

        # Boundaries of search volume
        #  [x, y, z] local coordinates
        lb = [search[0], search[2], search[4]]
        ub = [search[1], search[3], search[5]]

        # Init pool of workers
        pool = multiprocessing.Pool(multiprocessing.cpu_count())

        # arguments to be passed to timeFunction()
        args = (stns, w, kotc, setup, [ref_pos.lat, ref_pos.lon, ref_pos.elev], dataset, pool)

        # Particle Swarm Optimization
        # x_opt - optimal supracenter location
        # f_opt - optimal supracenter error

        if setup.restrict_to_trajectory:
            #cons = [lineConstraintx, lineConstrainty, lineConstraintz]
            x_opt, f_opt = pso(timeFunction, lb, ub, f_ieqcons=lineConstraint, args=args, swarmsize=setup.swarmsize, maxiter=setup.maxiter, \
                        phip=setup.phip, phig=setup.phig, debug=False, omega=setup.omega, minfunc=setup.minfunc, minstep=setup.minstep) 
        else:
            x_opt, f_opt = pso(timeFunction, lb, ub, args=args, swarmsize=setup.swarmsize, maxiter=setup.maxiter, \
                        phip=setup.phip, phig=setup.phig, debug=False, omega=setup.omega, minfunc=setup.minfunc, minstep=setup.minstep) 

        pool.close()
        pool.join()
        print('Done Searching')
    
    # If manual search
    else:

        # set point to 'optimal' supracenter to check residuals, convert to geographic
        x_opt = np.asarray(geo2Loc(ref_pos.lat, ref_pos.lon, ref_pos.elev, single_point[0], single_point[1], single_point[2]))

    # Get results for current Supracenter
    time3D, az, tf, r, motc, sotc = outputWeather(n_stations, x_opt, stns, setup, consts, \
                                                            [ref_pos.lat, ref_pos.lon, ref_pos.elev], dataset, output_name, s_name, kotc, w)


    # Find error for manual searches
    if len(single_point) != 0:

        f_opt = np.dot(w, np.absolute(sotc - motc)**setup.fit_type)/nwn

    # x, y distance from Supracenter to each station
    horz_dist = np.zeros(n_stations)
    for i in range(n_stations):

        horz_dist[i] = np.sqrt((x_opt[0] - xstn[i, 0])**2 + (x_opt[1] - xstn[i, 1])**2)/1000

    # Optimal Supracenter Data
    # Convert to geographic
    x_opt = loc2Geo(ref_pos.lat, ref_pos.lon, ref_pos.elev, x_opt)

    # Calculate and Set the Occurrence Time into HH:MM:SS
    time_diff = motc + ref_time.microsecond/1e6 + ref_time.second + ref_time.minute*60 + ref_time.hour*3600

    otc = (datetime.datetime.min + datetime.timedelta(seconds=time_diff)).time()

    global sup
    global errors

    # Potential Supracenter locations
    sup = np.delete(sup, 0, 0)

    # Error in each potential Supracenter
    errors = np.delete(errors, 0)

    #Filter errors
    try:
        while max(errors) > setup.max_error/100:
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

    return results

    # # scatter plot(s)
    # min_search, max_search = scatterPlot(single_point, n_stations, xstn, s_name, r, x_opt, \
    #                                         reported_points, search, output_name, ref_pos, sup, errors, tweaks, dataset)

    # # residual plot
    # residPlot(x_opt, s_name, xstn, r, output_name, n_stations)

    # # output results
    # outputText(min_search, max_search, single_point, ref_time, otc, kotc, x_opt, f_opt, n_stations, tweaks, s_name, xstn, \
    #                                                                     r, w, az, tf, time3D, horz_dist, output_name, tstn)