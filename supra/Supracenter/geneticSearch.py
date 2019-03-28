"""Finds the optimal Supracenter using a grid search algorithm and graphs and exports the data"""

import copy
import datetime
import multiprocessing
import time
import os

import matplotlib.pyplot as plt
import numpy as np

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Supracenter.stringDecode import stringDecode, createString, reproduce, crossover, mutatestrings
from supra.Supracenter.defgrid import defgrid
from supra.Supracenter.cyscan import cyscan
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.angleConv import loc2Geo, geo2Loc
from supra.Supracenter.plot import residPlot, scatterPlot, outputText, outputWeather

sup = [0, 0, 0]
errors = [0]

def minIndex(a):
    ''' Helper function: returns the index of the smallest value in the list (ignoring nan values)

    Arguments:
        a: [list] list to search for min index

    Returns:
        [int] index of minimum value in list a
    '''
    return list(a).index(np.nanmin(a))


def geneticSearch(x, y, z, stns, tweaks, sounding, w, ref_pos, kotc):
    """ Optimizes the best solution of a supracenter out of a group specified by the user. 
        The final solution is returned after the last iteration

    Argurments:
        x, y, z: [ndarray] containing the positions potential supracenters in the desired search area                                
        stns: [ndarray] matrix containing the position coordinates (UTM) of the receivers in the data set 
              along with their arrival times (1:n,X,Y,Z,t), taken from SInfo (see Station Info File)  
              Stns[1] = [Easting, Northing, Elevation, Arrival Time + Reference Time]
              Stns[2] = [Easting, Northing, Elevation, Arrival Time + Reference Time]
              etc...
        tweaks: [object] parameters that can be changed from the main menu
        sounding: [array] weather profile for the entire search area, for use with weatherInterp   
        w: [list] list of weights for each station (0 < w < 1, 1 - default, 0 - ignore station)
        ref_pos: [list] [latitude, longitude, elevation] reference coordinate to be used for the local coordinate system
                 by default, is the average latitude and longitude between the stations, and an elevation of 0.
        kotc: [float] user defined occurrence time in seconds after refTime                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
            
    Returns:
        minx, miny, minz: [float] position of best fitting supracenter (x,y,z)   
        BigErr: [float] degree of fitness to the data                  
        R: [list] individual station residuals                   
        mOTC: [float] mean time of occurance, in seconds after reference time                         
        BigAz, BigTf, BigTimes: [list] individual station angles and traveltime parameters
    """

    # Set initial error to be large
    BigErr = 1e10

    # Number of stations to use in the solution
    n_stations = len(stns)

    # Station arrival times
    tobs = stns[0:n_stations, 3]

    # Station positions
    xstn = stns[0:n_stations, 0:3]

    # Number of points to search with
    n_supra = tweaks.swarm_size

    # Number of iterations
    n_iter = tweaks.max_iter

    # Total weight of stations
    nw = sum(w)

    # Conversion to binary for binary strings
    pSize = np.log2(len(x))

    sLength = 3*int(pSize)

    # Init of positions with binary strings 
    supracenters = np.array(createString(sLength, int(n_supra)))
    
    # Init arrays

    # Calculated time to each station
    time3D = np.zeros((n_supra, n_stations))
    time3D[:] = np.nan

    # Initial Azimuth to each station
    az = np.zeros((n_supra, n_stations))
    az[:] = np.nan

    # Initial Takeoff Angle to each station
    tf = np.zeros((n_supra, n_stations))
    tf[:] = np.nan

    # Difference in station occurence times and calculated times
    sotc = np.zeros_like(tobs)

    # Weighted mean station occurence times
    motc = np.zeros_like(sotc)

    BigCenter = [0]*sLength
    BigTimes = [0]*n_stations
    BigAz = [0]*n_stations
    BigTf = [0]*n_stations
    Err = [0]*n_supra

    # Init Supracenter position
    source = [0, 0, 0]

    for i in range(n_iter):

        for n in range(n_supra):

            # Find Supracenter position indices from binary strings
            Sp = stringDecode(supracenters[n, :], int(pSize), 3)

            # Find Supracenter position from indices
            source[0] = x[Sp[0] - 1]
            source[1] = y[Sp[1] - 1]
            source[2] = z[Sp[2] - 1]

            # weight of each station
            wn = w

            # total of weights
            nwn = nw
            
            # Reset station travel times
            time3D[n, :] = np.nan

            for j in range(n_stations):

                # Tolerance on station height
                ach = 0

                # Increase tolerance up to a point
                while (ach <= 1000) and np.isnan(time3D[n, j]):

                    # tolerance of cyscan
                    tol = 10

                    # Weather interpolation
                    zProfile, points = getWeather([source[0], source[1], source[2]], xstn[j, :], tweaks.weather_type, ref_pos, sounding)

                    # Increase tolerance up to a point
                    while ((tol <= 1000) and np.isnan(time3D[n, j])):

                        # Find travel time and angles
                        time3D[n, j], az[n, j], tf[n, j] = cyscan(np.array(source), np.array(xstn[j, :]), zProfile,\
                                                            tol=tol, wind=tweaks.enable_winds, n_theta=tweaks.n_theta, \
                                                            n_phi=tweaks.n_phi, precision=tweaks.precision)

                        # Increase tolerance by steps
                        if (np.isnan(time3D[n, j])):
                            if tol == 1000:
                                tol = 1001
                            elif tol == 500:
                                tol = 1000
                            elif tol == 100:    
                                tol = 500
                            if tol == 10:   
                                tol = 100

                    # Search for rays passing over station
                    if (np.isnan(time3D[n, j])):              
                        if ach == 0:
                            ach = 500
                        elif ach == 500:
                            ach = 1000
                        elif ach == 1000: 
                            ach = 1001

                # If traveltime still could not be found set time to zero
                if (np.isnan(time3D[n, j])):
                    time3D[n, j] = 0

                # Residual in arrival times
                sotc[j] = tobs[j] - time3D[n, j] 

            # Mean time of occurence
            motc = np.divide(np.dot(wn, sotc), nwn)

            # See page 20 of SUPRACENTER
            # User-defined occurence time
            if kotc != None:
                Err[n] = np.dot(wn, np.absolute(sotc - np.array([kotc]*n_stations)))/nwn
            
            # Unknown occurence time
            else:
                Err[n] = np.dot(wn, np.absolute(sotc - motc))/nwn

            # Save points and errors for plotting

            global sup 
            global errors

            # check: Err == nan
            if not Err[n] > 0: 
                Err[n] = tweaks.max_error

            # check: large err
            elif Err[n] > tweaks.max_error:
                errors = np.hstack((errors, tweaks.max_error))
                sup = np.vstack((sup, [source[0], source[1], source[2]]))

                # Answer should not converge to actual errors
                Err[n] = tweaks.max_error

            else:
                errors = np.hstack((errors, Err[n]))
                sup = np.vstack((sup, [source[0], source[1], source[2]]))

            if tweaks.debug:
            # print out current search location
                print("Supracenter: {:10.2f} m x {:10.2f} m y {:10.2f} m z  Error: {:10.2f}"\
                    .format(float(source[0]), float(source[1]), float(source[2]), Err[n]))

        # Save the Lowest Error Position for all iterations
        if (min(Err) < BigErr):
            BigErr, Bi = min(Err), minIndex(Err)
            BigCenter[0:sLength] = supracenters[Bi, 0:sLength]
            BigTimes[0:n_stations] = time3D[Bi, 0:n_stations]
            BigAz[0:n_stations] = az[Bi, 0:n_stations]
            BigTf[0:n_stations] = tf[Bi, 0:n_stations]
    
        if (i < n_iter): 

            # Reproduction Phase
            newSupra = reproduce(supracenters, Err)
    
            # Crossover Phase
            newSupra = crossover(newSupra)
    
            # Mutate Strings
            newSupra = mutatestrings(newSupra)
    
            # Reevaluate Supracenter
            supracenters = newSupra
    
    print("Done Searching")

    # Indicies of the optimal supracenter position
    minx, miny, minz = stringDecode(BigCenter, pSize, 3)

    # find time difference for each station
    for j in range(n_stations):

        sotc[j] = tobs[j] - BigTimes[j]

    # Find residuals
    R = sotc - motc

    return minx, miny, minz, BigErr, R, motc, BigAz, BigTf, BigTimes     


def gLocator(search_area, s_info, weights, ref_pos, s_name, reported_points, ref_time, kotc,
                     tweaks, sounding, output_name, single_point, consts):
    """ gLocator finds the supracenter using a grid search method, and prints out the results

    Arguments:
        search_area: [list] boundaries of the search area in local coordinates [xmin, xmax, ymin, ymax, zmin, zmax]
        s_info: [ndarray] station location and signal arrival times, given in local coordinates and seconds
        weights: [list] station weights (0 - ignore station, 1 - default)
        ref_pos: [ndarray] mean station location, used for conversion to local coordinates
        s_name: [list] list of station names
        reported_points: [list] list of extra points to be plotted
        ref_time: [Object] datetime object of the time the station times are in reference to
        kotc: [list] user defined occurence time, if the time of fragmentation is known
        tweaks: [Object] user defined tweaks and tolerances

    """

    print('Data converted. Searching...')

    if len(single_point) != 0:
        print('Position from mean station location')

    n_stations = len(s_info)
    sax, say, saz, size = defgrid(search_area, size=256)

    # Station positions
    xstn = s_info[0:n_stations, 0:3]
    tstn = s_info[0:n_stations, 3]

    # Prevent search below stations
    if search_area[4] < max(xstn[:, 2]):

        # Must be just above the stations
        search_area[4] = max(xstn[:, 2]) + 0.0001

    # if no user-defined occurence time
    if kotc != None:
        kotc -= ref_time.hour*3600 + ref_time.minute*60 + ref_time.second + ref_time.microsecond/1e6
        
    # Automatic search
    if len(single_point) == 0:
        mx, my, mz, Err, Resid, motc, Az, Tf, time3D = geneticSearch(sax, say, saz, s_info, tweaks, sounding, weights, ref_pos, kotc)

        # Get position of minimum & output current results to user for appraisal
        # Ask user for next step.
        LX = sax[mx]
        LY = say[my]
        LZ = saz[mz]

    # Manual search
    else:

        LX, LY, LZ = geo2Loc(ref_pos[0], ref_pos[1], ref_pos[2], single_point[0], single_point[1], single_point[2])

        LX = single_point[0]
        LY = single_point[1]
        LZ = single_point[2]

    # Get results for current Supracenter
    time3D, Az, Tf, Resid, motc, sotc = outputWeather(n_stations, [LX, LY, LZ], s_info, tweaks, consts, ref_pos,\
                                                        sounding, output_name, s_name, kotc, weights)

    # Find error for manual searches
    if single_point != []:

        Err = np.dot(weights, np.absolute(sotc - motc)**tweaks.fit_type)/sum(weights)


    global sup
    global errors

    # Potential Supracenter locations
    sup = np.delete(sup, 0, 0)

    # Error in each potential Supracenter
    errors = np.delete(errors, 0)

    for i in range(1):
        a = []
        std_error = np.std(errors)
        lim = np.mean(errors) + 0*std_error

        for i in range(len(errors)):
            if errors[i] >= lim:
                a.append(i)

        errors = np.delete(errors, (a), axis=0)
        sup = np.delete(sup, (a), axis=0)

    # Calculate and Set the Occurance Time into hh:mm:ss
    time_diff = motc + ref_time.microsecond/1e6 + ref_time.second + ref_time.minute*60 + ref_time.hour*3600

    otc = (datetime.datetime.min + datetime.timedelta(seconds=time_diff)).time()
    
    x_opt = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], [LX, LY, LZ])

    # x, y distance from Supracenter to each station
    horz_dist = np.zeros(n_stations)
    for i in range(n_stations):
        horz_dist[i] = np.sqrt((x_opt[0] - xstn[i, 0])**2 + (x_opt[1] - xstn[i, 1])**2)/1000

    # scatter plot(s)
    min_search, max_search = scatterPlot(single_point, n_stations, xstn, s_name, Resid, x_opt, \
                                            reported_points, search_area, output_name, ref_pos, sup, errors, tweaks, sounding)
    
    # residual plot
    residPlot(x_opt, s_name, xstn, Resid, output_name, n_stations)

    # output results
    outputText(min_search, max_search, single_point, ref_time, otc, kotc, x_opt, Err, n_stations, tweaks, s_name, xstn, \
                                                                        Resid, weights, Az, Tf, time3D, horz_dist, output_name, tstn)