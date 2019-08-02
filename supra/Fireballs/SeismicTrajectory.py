""" Determine the fireball trajectory from seismic data.

Modified method of Pujol et al. (2005).

"""

from __future__ import print_function, division, absolute_import

import sys
import os
import multiprocessing
import threading
import traceback
import time
import copy
import math
import numpy as np
import scipy.signal
import scipy.optimize
import scipy.interpolate as interp
import argparse
import pyswarm
import datetime
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects
from functools import partial

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

try:
    # Python 2  
    import ConfigParser as configparser

except:
    # Python 3
    import configparser


from supra.Supracenter.cyscan2 import cyscan
from supra.Supracenter.cyzInteg import zInterp
import supra.Supracenter.cyweatherInterp
from supra.Supracenter.netCDFconv import storeHDF, storeNetCDFECMWF, storeNetCDFUKMO, readCustAtm, storeAus
from supra.Supracenter.fetchECMWF import fetchECMWF
from supra.Supracenter.fetchMERRA import fetchMERRA
from supra.Utils.AngleConv import loc2Geo, angle2NDE, geo2Loc
from supra.Supracenter.cyzInteg import zInteg
from supra.Supracenter.stationDat import readTimes
from supra.Utils.Classes import Position, Constants
from supra.Utils.pso import pso
from wmpl.Formats.CSSseismic import loadCSSseismicData
from wmpl.Utils.TrajConversions import date2JD, jd2Date, raDec2ECI, geo2Cartesian, cartesian2Geo, raDec2AltAz, eci2RaDec, latLonAlt2ECEF, ecef2ENU, enu2ECEF, ecef2LatLonAlt
from wmpl.Utils.Math import vectMag, vectNorm, rotateVector, meanAngle
from wmpl.Utils.Plotting import Arrow3D, set3DEqualAxes
from wmpl.Utils.PlotMap import GroundMap
from wmpl.Utils.PlotCelestial import CelestialPlot

def findPoints(setup):    
    GRID_SPACE = 50
    MIN_HEIGHT = 15000
    MAX_HEIGHT = 60000

    u = setup.trajectory.vector.xyz
    ground_point = setup.trajectory.pos_f.xyz
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

def parseWeather(setup, t=0):
    consts = Constants()
     # Parse weather type
    setup.weather_type = setup.weather_type.lower()
    if setup.weather_type not in ['custom', 'none', 'ecmwf', 'merra', 'ukmo', 'binary']:
        print('INI ERROR: [Atmospheric] weather_type must be one of: custom, none, ecmwf, merra, ukmo, binary')
        sys.exit()

    # Custom weather data    
    if setup.weather_type == 'custom':

        # try:
        sounding = readCustAtm(setup.sounding_file, consts)
        # except:
        #     print('ERROR: Could not read custom atmosphere profile!')
        #     exit()

        # Check file name
        if len(sounding) < 2:
            print('FILE ERROR: file must have at least two data points!')
            sys.exit()


    # MERRA-2
    elif setup.weather_type == 'merra':

        # Check file type
        if '.nc' not in setup.sounding_file:
            print("FILE ERROR: custom data set must be a .nc file!")
            sys.exit()

        # # Run data fetching script
        # if setup.get_data == True:
        #     print('Getting data...')
        #     fetchMERRA(setup)

        try:
            sounding = storeHDF(setup.sounding_file, consts)
        except:
            print('ERROR: could not read merra atmosphere profile!')
            exit()

    # ECMWF
    elif setup.weather_type == 'ecmwf':

        # Check file type
        if '.nc' not in setup.sounding_file:
            print("FILE ERROR: custom data set must be a .nc file from ECMWF!")
            sys.exit()

        # # Run data fetching script
        # if setup.get_data == True:
        #     fetchECMWF(setup, setup.sounding_file)

        # GetIRISData/MakeIRISPicks

        try:

            #Get closest hour
            start_time = (setup.fireball_datetime.hour + np.round(setup.fireball_datetime.minute/60) + t)%24
            sounding = storeNetCDFECMWF(setup.sounding_file, setup.lat_centre, setup.lon_centre, consts, start_time=start_time)

        # SeismicTrajectory
        except:

            try:
                start_time = (setup.fireball_datetime.hour + t)%24
                sounding = storeNetCDFECMWF(setup.sounding_file, setup.x0, setup.y0, consts, start_time=start_time)
            except:
                print("ERROR: Unable to use weather file, or setup.start_datetime/setup.atm_hour is set up incorrectly. Try checking if sounding_file exists")
                print(traceback.format_exc())
                exit()

        

    # UKMO
    elif setup.weather_type == 'ukmo':
        # Check file type
        if '.nc' not in setup.sounding_file:
            print("FILE ERROR: custom data set must be a .nc file from UKMO!")
            sys.exit()
        
        try:
            sounding = storeNetCDFUKMO(setup.sounding_file, setup.search_area, consts)
        except:
            print('ERROR: Could not read UKMO atmosphere profile!')
            exit()

    elif setup.weather_type == 'binary':

        
        sounding = storeAus(consts)

    else:

        # Sample fake weather profile, the absolute minimum that can be passed
        sounding = np.array([[    0.0, setup.v_sound, 0.0, 0.0],
                             [    0.0, setup.v_sound, 0.0, 0.0],
                             [99999.0, setup.v_sound, 0.0, 0.0]])

    return sounding

def mergeChannels(seismic_data):
    """ Merge seismic data from the same channels which are fragmented into several chunks which start at a
        different time

    """

    merged_data = []
    merged_indices = []

    for i, entry1 in enumerate(seismic_data):

        # Skip the entry if it was already merged
        if i in merged_indices:
            continue

        # Unpack the loaded seismic data
        site1, w1, time_data1, waveform_data1 = entry1

        merged_sites = [site1]
        merged_wfdisc = [w1]
        merged_time = [time_data1]
        merged_waveform = [waveform_data1]


        for j, entry2 in enumerate(seismic_data):

            # Skip the same entries
            if i == j:
                continue

            # Skip the entry if it was already merged
            if j in merged_indices:
                continue

            # Unpack the seismic data
            site2, w2, time_data2, waveform_data2 = entry2


            # Check if the channels are the same
            if (w1.sta == w2.sta) and (w1.chan == w2.chan):

                merged_sites.append(site2)
                merged_wfdisc.append(w2)
                merged_time.append(time_data2)
                merged_waveform.append(waveform_data2)

                merged_indices.append(j)


        # Add all merged data to the list, but keep the first site and wfdisc info as the reference one
        merged_data.append([merged_sites, merged_wfdisc, merged_time, merged_waveform])


    return merged_data



def wrapRaDec(ra, dec):
    """ Wraps RA and Dec into their limits. """

    # Wrap RA into [0, 2pi) range
    ra = ra%(2*np.pi)

    # Wrap Dec into [-pi/2, pi/2]
    dec = (dec + np.pi/2)%np.pi - np.pi/2

    return ra, dec


def timeOfArrival(stat_coord, x0, y0, t0, v, azim, zangle, setup, points, u, sounding=[], ref_loc=Position(0, 0, 0), \
                        travel=False, fast=False, theo=False):
    """ Calculate the time of arrival at given coordinates in the local coordinate system for the given
        parameters of the fireball.

    Arguments:
        stat_coord: [3 element ndarray] Coordinates of the station in the local coordinate system.
        x0: [float] Intersection with the X axis in the local coordinate system (meters).
        y0: [float] Intersection with the Y axis in the local coordinate system (meters).
        t0: [float] Time when the trajectory intersected the reference XY plane (seconds), offset from 
            some reference time.
        v: [float] Velocity of the fireball (m/s).
        azim: [float] Fireball azimuth (+E of due N). (radians)
        zangle: [float] Zenith angle. (radians)
        setup: [Object] Object containing all user-defined parameters
        sounding: [ndarray] atmospheric profile of the search area
        travel: [boolean] switch to only return the travel time

    Return:
        ti: [float] Balistic shock time of arrival to the given coordinates (seconds).

    """    

    #azim = (np.pi - azim)%(2*np.pi)

    # Calculate the mach angle
    #cos(arcsin(x)) = sqrt(1 - x^2)

    beta = math.sqrt(1 - (setup.v_sound/v/1000)**2)

    # Difference from the reference point on the trajectory and the station
    b = stat_coord - np.array([x0, y0, 0])

    # Calculate the distance along the trajectory
    dt = abs(np.dot(b, -u))

    # Calculate the distance perpendicular to the trajectory
    dp = math.sqrt(abs(vectMag(b)**2 - dt**2))

    # No winds
    if setup.weather_type == 'none':

        if travel:
            # travel from trajectory only
            ti = (dp*beta)/setup.v_sound
        else:
            if theo:
                # Calculate the time of arrival
                ti = t0 + dt/v + (dp*beta)/setup.v_sound
            else:
                # Calculate the time of arrival
                ti = t0 - dt/v + (dp*beta)/setup.v_sound   

    # Winds
    else:

        R = waveReleasePointWinds(stat_coord, setup, sounding, ref_loc, points, u)

        if travel:
            # travel from trajectory only
            ti = R[3]*beta

        else:
            if theo:
                # Calculate time of arrival
                ti = t0 +dt/v + R[3]*beta
            else:
                # Calculate time of arrival
                ti = t0 -dt/v + R[3]*beta

    return ti


# methods for finding angle between two vectors
# import math
# def dotproduct(v1, v2):
#   return sum((a*b) for a, b in zip(v1, v2))

# def length(v):
#   return math.sqrt(dotproduct(v, v))

# def angle(v1, v2):
#   return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))

def waveReleasePointWinds(stat_coord, setup, sounding, ref_loc, points, u):
    #azim = (np.pi - azim)%(2*np.pi)
    # Break up the trajectory into points

    #ANGLE_TOL = 25 #deg


    # Trajectory vector
    #u = np.array([-np.cos(azim)*np.sin(zangle), np.sin(azim)*np.sin(zangle), -np.cos(zangle)])
    a = (len(points))
    T = [0]*a
    az = [0]*a
    tf = [0]*a
    prop_mag = [0]*a

    # Cut down atmospheric profile to the correct heights, and interp
    z_profile, _ = supra.Supracenter.cyweatherInterp.getWeather(points[0], stat_coord, setup.weather_type, \
                     ref_loc, copy.copy(sounding), convert=True)
    
    D = np.array(stat_coord)

    cyscan_res = []

    # Compute time of flight residuals for all stations
    for i in range(a):
        S = np.array(points[i])

        zProfile = zInteg(D[2], S[2], z_profile)
        zProfile = zInterp(D[2], S[2], zProfile, div=37)
        cyscan_res.append(cyscan(S, D, zProfile, wind=setup.enable_winds, n_theta=setup.n_theta, n_phi=setup.n_phi, h_tol=setup.h_tol, v_tol=setup.v_tol))

    cyscan_res = np.array(cyscan_res)
    T = cyscan_res[:, 0]
    az = cyscan_res[:, 1]
    tf = cyscan_res[:, 2]

    az = np.radians(az)
    tf = np.radians(180 - tf)

    v = [0]*a
    for ii in range(a):
        v[ii] = np.array([np.sin(az[ii])*np.sin(tf[ii]), np.cos(az[ii])*np.sin(tf[ii]), -np.cos(tf[ii])])
        prop_mag[ii] = np.dot(u, v[ii])

    # Minimize angle to 90 degrees from trajectory

    # plt.plot(points[:, 2], np.degrees(np.arccos(prop_mag)))
    # plt.show()
    #np.save(setup.working_directory + str(stat_coord), prop_mag_raw)
    prop_mag = np.absolute(prop_mag)

    # for ii, element in enumerate(prop_mag):
    #     if element > np.cos(np.radians(ANGLE_TOL)):
    #         prop_mag[ii] = np.nan

    try:
        best_indx = np.nanargmin(prop_mag)
        #best_indx = np.nanargmax(diff)
    except:
        return np.array([np.nan, np.nan, np.nan, np.nan])

    S = np.array(points[best_indx])
    T = np.array(T[best_indx])
    # Return point
    return np.append(S, T)

def waveReleasePoint(stat_coord, x0, y0, t0, v, azim, zangle, v_sound):
    """ Calculate the point on the trajectory from which the balistic wave was released and heard by the given
        station.

    Arguments:
        stat_coord: [3 element ndarray] Coordinates of the station in the local coordinate system.
        x0: [float] Intersection with the X axis in the local coordinate system (meters).
        y0: [float] Intersection with the Y axis in the local coordinate system (meters).
        t0: [float] Time when the trajectory intersected the reference XY plane (seconds), offset from 
            some reference time.
        v: [float] Velocity of the fireball (m/s).
        azim: [float] Fireball azimuth (+E of due N).
        zangle: [float] Zenith angle.
        v_sound: [float] Average speed of sound (m/s).

    Return:
        traj_point: [3 element ndarray] Location of the release point in the local coordinate system.

    """

    # back azimuth
    #azim = (np.pi - azim)%(2*np.pi)


    # Calculate the mach angle
    beta = np.arcsin(v_sound/v)

    # Trajectory vector
    u = np.array([np.sin(azim)*np.sin(zangle), np.cos(azim)*np.sin(zangle), -np.cos(zangle)])

    # Difference from the reference point on the trajectory and the station
    b = stat_coord - np.array([x0, y0, 0])

    # Calculate the distance along the trajectory
    dt = np.abs(np.dot(b, -u))

    # Calculate the distance perpendicular to the trajectory
    dp = np.sqrt(vectMag(b)**2 - dt**2)

    # Vector from (x0, y0) to the point of wave release
    r = -u*(dt + dp*np.tan(beta))
    #r = -u*dt

    # Position of the wave release in the local coordinate system
    traj_point = np.array([x0, y0, 0]) + r


    return traj_point



def timeResidualsAzimuth(params, stat_coord_list, arrival_times, setup, sounding, v_fixed=None, \
        print_residuals=False, pool=[]):
    """ Cost function for seismic fireball trajectory optimization. The function uses 

    Arguments:
        params: [list] Estimated parameters: x0, t0, t0, v, azim, elev.
        stat_coord_list: [list of ndarrays] A list of station coordinates (x, y, z) in the reference coordinate system.
        arrival_times: [list] A list of arrival times of the sound wave to the seismic station (in seconds 
            from some reference time).
        setup: [Object] Object containing all user-defined parameters
        sounding: [ndarray] atmospheric profile of the search area

    Keyword arguments:
        azim_off: [float] Azimuth around which the given values are centred. If None (default), it is assumed 
            that the azimuth is calculated +E of due North.
        v_fixed: [float] Use a fixed velocity. Set to None to ignore

    """
    # Unpack estimated parameters
    x0, y0, t0, v, azim, zangle = params
    
    # if v_fixed != '':

    #     # Keep the fireball velocity fixed
    #     v = v_fixed
    if v_fixed != None:
        v = float(v_fixed)

    ref_pos = Position(setup.lat_centre, setup.lon_centre, 0)
    cost_value = 0
    
    points = findPoints(setup)
    u = setup.trajectory.vector.xyz
    
    # Go through all arrival times
    for i, (t_obs, stat_coord) in enumerate(zip(arrival_times, stat_coord_list)):
        ### Calculate the difference between the observed and the prediced arrival times ###
        ######################################################################################################

        # Calculate the time of arrival
        ti = timeOfArrival(stat_coord, x0, y0, t0, v, np.radians(azim), np.radians(zangle), setup, points, u, sounding=sounding, ref_loc=ref_pos, theo=True)

        # Smooth approximation of l1 (absolute value) loss

        cost_value += 2*((1 + (t_obs - ti)**2)**0.5 - 1)

        if np.isnan(ti): 
            continue     
    
        ######################################################################################################
    # Save points for plotting later

    # Save points and errors for plotting

    return cost_value



def convertStationCoordinates(station_list, ref_indx):
    """ Converts the coordinates of stations into the local coordinate system, the origin being one of the
        given stations.
    
    Arguments:
        station_list: [list] A list of stations and arrival times, each entry is a tuple of:
            (name, lat, lon, elev, arrival_time_jd), where latitude and longitude are in radians, the 
            zangle is in meters and Julian date of the arrival time.
        ref_indx: [int] Index of the reference station which will be in the origin of the local coordinate
            system.
    """


    # Select the location of the reference station
    _, _, _, lat0, lon0, elev0, _,  _, _ = station_list[ref_indx]

    stat_coord_list = []

    # Calculate the local coordinates of stations at the given time
    for i, entry in enumerate(station_list):

        _, _, _, lat, lon, elev,_,_, _ = entry

        # Convert geographical to local coordinates
        stat_coord = geo2Loc(lat0, lon0, elev0, lat, lon, elev)

        ######

        stat_coord_list.append(stat_coord)


    return stat_coord_list



def estimateSeismicTrajectoryAzimuth(station_list, setup, sounding, p0=None, azim_range=None, elev_range=None, \
        v_fixed=None, allTimes=[], ax=None):
    """ Estimate the trajectory of a fireball from seismic/infrasound data by modelling the arrival times of
        the balistic shock at given stations.

    The method is implemented according to Pujol et al. (2005) and Ishihara et al. (2003).

    Arguments:
        station_list: [list] A list of stations and arrival times, each entry is a tuple of:
            (name, lat, lon, elev, arrival_time_jd), where latitude and longitude are in radians, the 
            zangle is in meters and Julian date of the arrival time.
        setup: [Object] Object containing all user-defined parameters
        sounding: [ndarray] atmospheric profile of the search area

    Keyword arguments:
        p0: [6 element ndarray] Initial parameters for trajectory estimation:
            p0[0]: [float] p0, north-south offset of the trajectory intersection with the ground (in km, +S).
            p0[1]: [float] y0, east-west offset of the trajectory intersection with the ground (in km, +E).
            p0[2]: [float] Time when the trajectory was at point (p0, y0), reference to the reference time 
                (seconds).
            p0[3]: [float] Velocity of the fireball (km/s).
            p0[4]: [float] Initial azimuth (+E of due north) of the fireball (radians).
            p0[5]: [float] Initial zenith angle of the fireball (radians).
        azim_range: [list of floats] (min, max) azimuths for the search, azimuths should be +E of due North
            in degrees. If the range of azmiuths traverses the 0/360 discontinuity, please use negative
            values for azimuths > 180 deg!
        elev_range: [list of floats] (min, max) zangles for the search, in degrees. The zangle is 
            measured from the horizon up.
        v_fixed: [float or None] value to restrict the velocity of the meteor. If set to None, there is no restriction

    """

    if ax is None:
        ax = plt.gca()

    # Extract all Julian dates
    jd_list = [entry[6] for entry in station_list]
    pick_time = [float(entry[7]) for entry in station_list]
    station_no = [int(float(entry[8])) for entry in station_list]

    # Calculate the arrival times as the time in seconds from the earliest JD

    jd_ref = min(jd_list)
    jd_list = np.array(jd_list)

    # Get the index of the first arrival station
    first_arrival_indx = np.argwhere(jd_list == jd_ref)[0][0]

    ref_pick_time = pick_time[first_arrival_indx]

    # Convert station coordiantes to local coordinates, with the station of the first arrival being the
    # origin of the coordinate system
    stat_coord_list = convertStationCoordinates(station_list, first_arrival_indx)

    ref_pos = Position(setup.lat_centre, setup.lon_centre, 0)

    setup.pos_min.pos_loc(ref_pos)
    setup.pos_max.pos_loc(ref_pos)

    bounds = [
        (setup.pos_min.x, setup.pos_max.x), # X0
        (setup.pos_min.y, setup.pos_max.y), # Y0
        (setup.t_min, setup.t_max), # t0
        (setup.v_min, setup.v_max), # Velocity (m/s)
        (setup.azimuth_min.deg, setup.azimuth_max.deg),     # Azimuth
        (setup.zenith_min.deg, setup.zenith_max.deg)  # Zenith angle
        ]

    print("Bounds:")
    print("x      : {:+8.2f} - {:+8.2f} km".format(setup.pos_min.x/1000, setup.pos_max.x/1000))
    print("y      : {:+8.2f} - {:+8.2f} km".format(setup.pos_min.y/1000, setup.pos_max.y/1000))
    print("t      : {:8.2f} - {:8.2f} s".format(setup.t_min, setup.t_max))
    print("v      : {:8.2f} - {:8.2f} km/s".format(setup.v_min/1000, setup.v_max/1000))
    print("Azimuth: {:8.2f} - {:8.2f} deg fN".format(setup.azimuth_min.deg, setup.azimuth_max.deg))
    print("Zenith : {:8.2f} - {:8.2f} deg".format(setup.zenith_min.deg, setup.zenith_max.deg))


    # Extract lower and upper bounds
    lower_bounds = [bound[0] for bound in bounds]
    upper_bounds = [bound[1] for bound in bounds]

    class MiniResults():
        def __init__(self):
            self.x = None

    # Run PSO several times and choose the best solution
    solutions = []
    for i in range(setup.run_times):

        # Use PSO for minimization
        x, fopt, particles, errors = pso(timeResidualsAzimuth, lower_bounds, upper_bounds, args=(stat_coord_list, \
            pick_time, setup, sounding, v_fixed), maxiter=setup.maxiter, swarmsize=setup.swarmsize, \
            phip=setup.phip, phig=setup.phig, debug=False, omega=setup.omega, \
            processes=multiprocessing.cpu_count(), particle_output=True)


        solutions.append([x, fopt])

    if setup.perturb == True:

        #allTimes = [perturb_no, station_no, ball/frag, frag_no]
        perturb_times = allTimes.shape[0]

        p_arrival_times = allTimes[:, station_no, 0, 0] - float(ref_pick_time)

        x_perturb = [0]*perturb_times
        fopt_perturb = [0]*perturb_times

        for i in range(perturb_times):
            

            #remove nan
            #p_arrival_times[i] = [j for j in p_arrival_times[i] if j != j]

            # Use PSO for minimization, with perturbed arrival times
            x_perturb[i], fopt_perturb[i] = pso(timeResidualsAzimuth, lower_bounds, upper_bounds, args=(stat_coord_list, \
                p_arrival_times[i], setup, sounding, v_fixed), maxiter=setup.maxiter, swarmsize=setup.swarmsize, \
                phip=setup.phip, phig=setup.phig, debug=False, omega=setup.omega, processes=multiprocessing.cpu_count())
            
            if i == 0:
                print('Computational  best estimation', fopt_perturb[i])
            else:
                print('Perturbation', i, 'best estimation', fopt_perturb[i])
    else:
        x_perturb, fopt_perturb = [], []

    # Choose the solution with the smallest residuals
    fopt_array = np.array([fopt for x_val, fopt_val in solutions])
    best_indx = np.argmin(fopt_array)

    x, fopt = solutions[best_indx]

    res = MiniResults()
    res.x = x


    # Extract estimated parameters
    x0, y0 = res.x[:2]
    t0 = res.x[2]
    v_est = res.x[3]
    azim, zangle = res.x[4:]


    ref_pos = Position(setup.lat_centre, setup.lon_centre, 0)


    lat_fin, lon_fin, _ = loc2Geo(ref_pos.lat, ref_pos.lon, ref_pos.elev, [x0, y0, 0])

    # Print the time residuals per every station
    timeResidualsAzimuth(res.x, stat_coord_list, pick_time, setup, sounding, v_fixed=v_fixed, 
        print_residuals=True)

    # Plot the stations and the estimated trajectory
    residuals = plotStationsAndTrajectory(station_list, [x0, y0, t0, v_est/1000, azim, zangle], setup, sounding, x_perturb=x_perturb, ax=ax)


    final_pos = Position(lat_fin, lon_fin, 0)
    results = [final_pos, ref_pos, t0, v_est, azim, zangle, residuals]

    results = []
    print("Results:")
    print("Trajectory (Nominal)           | Lat {:+10.4f} N Lon {:+10.4f} E t {:5.2f} s v {:7.4f} km/s Azimuth {:6.2f} deg fN Zenith {:5.2f} Error {:10.4f}"\
                        .format(lat_fin, lon_fin, t0, v_est/1000, azim, zangle, fopt))
    results.append([lat_fin, lon_fin, t0, v_est/1000, azim, zangle, fopt])
    for i in range(perturb_times):
        lat_fin, lon_fin, _ = loc2Geo(ref_pos.lat, ref_pos.lon, ref_pos.elev, [x_perturb[i][0], x_perturb[i][1], 0])
        print("Trajectory (Perturbation {:4d}) | Lat {:+10.4f} N Lon {:+10.4f} E t {:5.2f} s v {:7.4f} km/s Azimuth {:6.2f} deg fN Zenith {:5.2f} Error {:10.4f}"\
                        .format(i, lat_fin, lon_fin, x_perturb[i][2], x_perturb[i][3]/1000,\
                         x_perturb[i][4], x_perturb[i][5], fopt_perturb[i]))   
        results.append([lat_fin, lon_fin, x_perturb[i][2], x_perturb[i][3]/1000, x_perturb[i][4], x_perturb[i][5], fopt_perturb[i]])

    particles = np.array(particles)
    errors = np.array(particles)

    return particles, errors, results

def addPressure(sounding):
    """ Helper function: adds a model pressure based off of altitude to the weather data so that a darkflight.c
    ini file can be exported. It takes both lists of heights and compressed/stretches the model list to fit the
    sounding profile.

    Arguments:
        sounding [ndarray]: sounding profile to add pressure levels to.

    Returns:
        p_compress [ndarray]: list of model pressures which align with the heights from the sounding profile.
    """

    # model pressure file, taken from here: https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels
    model_pressure = 'wmpl/Fireballs/pressure_conv.txt'

    with open(model_pressure) as f:

        data = np.array([0, 0])

        # Parse file contents
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '')

            # Split the line by the delimiter
            line = line.split()

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = float(entry.strip())

            # Add the contents of the line to the data list
            data = np.vstack((data, line))

        # First row was all zeroes
        data = np.delete(data, 0, 0)

    data[:, 0], data[:, 1] = data[:, 1], data[:, 0].copy()
    h = sounding[:, 0]

    pressure_list = zInteg(h[0], h[-1], data)

    # stretch/compress model list to fit sounding list
    p_interp = interp.interp1d(np.arange(pressure_list[:, 0].size), pressure_list[:, 1])
    p_compress = p_interp(np.linspace(0, pressure_list[:, 0].size - 1, h[:].size))

    return p_compress

def getTrajectoryVector(azim, zangle):
        
    # Calculate the trajectory vector
    u = np.array([np.sin(azim)*np.sin(zangle), np.cos(azim)*np.sin(zangle), -np.cos(zangle)])

    return u

def plotStationsAndTrajectory(station_list, params, setup, sounding, x_perturb=[], ax=None):
    """ Plots the stations in the local coordinate system and the fitted trajectory. 
        
    Arguments:
        station_list: [list] contains station names, locations, and travel times of the signal
        params: [list] list of the six parameters of the fireball, x0, y0, t0, v, azim, zangle  
        setup: [Object] Object containing all user-defined parameters
        sounds: [ndarray] atmospheric profile of the search area

    """
    
    if ax is None:
        ax = plt.gca()

    # Unpack estimated parameters
    x0, y0, t0, v, azim, zangle = params

    perturb_times = len(x_perturb)

    x_p = [0]*perturb_times
    y_p = [0]*perturb_times
    t_p = [0]*perturb_times
    v_p = [0]*perturb_times
    az_p = [0]*perturb_times 
    ze_p = [0]*perturb_times 
    
    if len(x_perturb) != 0:

        for i in range(perturb_times):
            x_p[i] = x_perturb[i][0] 
            y_p[i] = x_perturb[i][1] 
            t_p[i] = x_perturb[i][2]
            v_p[i] = x_perturb[i][3]
            az_p[i] = x_perturb[i][4]
            ze_p[i] = x_perturb[i][5]

    # Extract all Julian dates
    jd_list = [entry[6] for entry in station_list]
    pick_time = [float(entry[7]) for entry in station_list]
    # Calculate the arrival times as the time in seconds from the earliest JD
    ref_indx = np.argmin(jd_list)
    jd_list = np.array(jd_list)
    stat_obs_times_of_arrival = pick_time#(jd_list - jd_ref)*86400.0

    # Convert station coordiantes to local coordinates, with the station of the first arrival being the
    # origin of the coordinate system
    stat_coord_list = convertStationCoordinates(station_list, ref_indx)

    # Extract coordinates of the reference station
    #lat0, lon0, elev0 = station_list[ref_indx][3:6]
    lat0, lon0, elev0 = setup.lat_centre, setup.lon_centre, 0

    azim = np.radians(azim)
    zangle = np.radians(zangle)

    # Calculate the trajectory vector
    u = np.array([np.sin(azim)*np.sin(zangle), np.cos(azim)*np.sin(zangle), -np.cos(zangle)])

    # Calculate modelled times of arrival and points of wave releasefor every station
    stat_model_times_of_arrival = []
    wave_release_points = []

    points = findPoints(setup)
    u = setup.trajectory.vector.xyz
    for stat_coord in stat_coord_list:

        # Calculate time of arrival
        ti = timeOfArrival(stat_coord, x0, y0, t0, v, azim, zangle, setup, points, u, sounding=sounding, ref_loc=Position(lat0, lon0, elev0))
        stat_model_times_of_arrival.append(ti)

        # Calculate point of wave release
        traj_point = waveReleasePoint(stat_coord, x0, y0, t0, v, azim, zangle, setup.v_sound)
        wave_release_points.append(traj_point)

    stat_model_times_of_arrival = np.array(stat_model_times_of_arrival)
    wave_release_points = np.array(wave_release_points)

    # Get the release points with the highest and the lowest release height
    high_point = np.argmax(wave_release_points[:, 2])
    low_point = np.argmin(wave_release_points[:, 2])

    # Terminal output

    print('Wave release:')
    print(wave_release_points[:, 2])
    print(' - top:', wave_release_points[high_point, 2], 'km')
    print(' - bottom:', wave_release_points[low_point, 2], 'km')
    print()
    print('Residuals:')
    print('{:>10s}: {:6s}'.format('Station', 'res (s)'))
    
    sqr_res_acc = 0
    res = []
    residuals = []

    # Print the residuals per station
    for stat_entry, toa_obs, toa_model in zip(station_list, stat_obs_times_of_arrival, \
        stat_model_times_of_arrival):

        # Station network code and name
        station_name = stat_entry[1] + '-' + stat_entry[2]

        res.append(toa_obs - toa_model)
        print('{:>10s}: {:+06.2f}s'.format(station_name, toa_obs - toa_model))

        sqr_res_acc += (toa_obs - toa_model)**2
        residuals.append([station_name, toa_obs - toa_model])

    print('RMS:', np.sqrt(sqr_res_acc/len(station_list)), 's')

    # Corner locations as boundaries for the contour plot
    corner_times = [0, 0, 0, 0]
    corner_coords = np.array([[setup.x_min, setup.y_min, 0.0],
                     [setup.x_min, setup.y_max, 0.0],
                     [setup.x_max, setup.y_min, 0.0],
                     [setup.x_max, setup.y_max, 0.0]])

    # The maximum time in the contour plot must be in one of the corners
    for i, coor in enumerate(corner_coords):
        ti = timeOfArrival(coor, x0, y0, t0, v, azim, zangle, setup, points, u, ref_loc=Position(lat0, lon0, elev0), sounding=sounding)
        corner_times[i] = ti
    
    # Determine the maximum absolute time of arrival (either positive of negative)
    # toa_abs_max = np.max([np.abs(np.min(stat_model_times_of_arrival)), np.max(stat_model_times_of_arrival)])
    toa_abs_max = np.max(np.abs(corner_times))

    # Calcualte absolute observed - calculated residuals
    toa_residuals = np.abs(stat_obs_times_of_arrival - stat_model_times_of_arrival)

    # Determine the maximum absolute residual
    toa_res_max = np.max(toa_residuals)


    ### PLOT 3D ###
    # ##########################################################################################################

    # # # Setup 3D plot
    # # fig = plt.figure()
    # # ax = fig.gca(projection='3d')


    # # # Plot the stations except the reference
    # # station_mask = np.ones(len(x), dtype=bool)
    # # station_mask[ref_indx] = 0
    # # ax.scatter(x[station_mask], y[station_mask], z[station_mask], c=stat_model_times_of_arrival[station_mask], \
    # #     depthshade=0, cmap='viridis', vmin=0.00, vmax=toa_abs_max, edgecolor='k')

    # # # Plot the reference station
    # # ax.scatter(x[ref_indx], y[ref_indx], z[ref_indx], c='k', zorder=5)

    # lats = [np.degrees(float(i)) for i in [row[3] for row in station_list]]
    # lons = [np.degrees(float(i)) for i in [row[4] for row in station_list]]
    # elev = [(float(i)) for i in [row[5] for row in station_list]]

    # ############### PLOT STATIONS TOGGLE
    # PLOT = False
    # # Plot the stations (the colors are observed - calculated residuals)
    # if PLOT:
    #     ax.scatter(lats, lons, elev, c=toa_residuals, cmap='inferno_r', \
    #         edgecolor='0.5', linewidths=1, vmin=0, vmax=toa_res_max)

    # # plt.colorbar(stat_scat, label='abs(O - C) (s)')

    # # Plot the trajectory intersection with the ground
    # # ax.scatter(x0, y0, 0, c='g')

    # # Plot the lowest and highest release points
    # wrph_x, wrph_y, wrph_z = wave_release_points[high_point]
    # wrpl_x, wrpl_y, wrpl_z = wave_release_points[low_point]
    # # ax.scatter(wrph_x, wrph_y, wrph_z)
    # # ax.scatter(wrpl_x, wrpl_y, wrpl_z)
  
    # print('Trajectory vector:', u)

    # # Get the limits of the plot
    # x_min, x_max = ax.get_xlim()
    # y_min, y_max = ax.get_ylim()
    # z_min, z_max = ax.get_zlim()

    # x_min = setup.x_min
    # x_max = setup.x_max
    # y_min = setup.y_min
    # y_max = setup.y_max

    # # Get the maximum range of axes
    # x_range = abs(x_max - x_min)
    # y_range = abs(y_max - y_min)
    # z_range = abs(z_max - z_min)

    # traj_len = 0.10*max([x_range, y_range, z_range])

    # # Calculate the beginning of the trajectory
    # x_beg = x0 - traj_len*u[0]
    # y_beg = y0 - traj_len*u[1]
    # z_beg = -traj_len*u[2]


    # top_lat, top_lon, top_elev = loc2Geo(lat0, lon0, elev0, [x_beg, y_beg, z_beg])
    # bot_lat, bot_lon, bot_elev = loc2Geo(lat0, lon0, elev0, [x0, y0, 0])
    # # Plot the trajectory
    # ax.plot([top_lat, bot_lat], [top_lon, bot_lon], [top_elev, bot_elev], c='b')

    # if len(x_perturb) != 0:
    #     u_p = [0]*perturb_times

    #     for i in range(perturb_times):
    #         u_p[i] = getTrajectoryVector(np.radians(az_p[i]), np.radians(ze_p[i]))

    #         x_beg_p = x_p[i] - traj_len*u_p[i][0]
    #         y_beg_p = y_p[i] - traj_len*u_p[i][1]
    #         z_beg_p = -traj_len*u_p[i][2]

    #         top_lat, top_lon, top_elev = loc2Geo(lat0, lon0, elev0, [x_beg_p, y_beg_p, z_beg_p])
    #         bot_lat, bot_lon, bot_elev = loc2Geo(lat0, lon0, elev0, [x_p[i], y_p[i], 0])
    #         ax.plot([top_lat, bot_lat], [top_lon, bot_lon], [top_elev, bot_elev], c='b', alpha=0.4)

    # Plot wave release trajectory segment
    #ax.plot([wrph_x, wrpl_x], [wrph_y, wrpl_y], [wrph_z, wrpl_z], color='red', linewidth=2)

    # ### Plot the boom corridor ###
    # x_data = np.linspace(x_min, x_max, setup.contour_res)
    # y_data = np.linspace(y_min, y_max, setup.contour_res)
    # xx, yy = np.meshgrid(x_data, y_data)

    # # Make an array of all plane coordinates
    # plane_coordinates = np.c_[xx.ravel(), yy.ravel(), np.zeros_like(xx.ravel())]

    # times_of_arrival = np.zeros_like(xx.ravel())
    # print("Creating contour plot...")

    # # Calculate times of arrival for each point on the reference plane
    # for i, plane_coords in enumerate(plane_coordinates):

    #     ti = timeOfArrival(plane_coords, x0, y0, t0, v, azim, zangle, setup, points, u, sounding=sounding, ref_loc=Position(lat0, lon0, elev0))

    #     times_of_arrival[i] = ti

    # times_of_arrival = times_of_arrival.reshape(setup.contour_res, setup.contour_res)

    # # Determine range and number of contour levels, so they are always centred around 0
    # # toa_abs_max = np.max([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
    # levels = np.linspace(0.00, toa_abs_max, 50)

    # # Plot colorcoded times of arrival on the surface
    # toa_conture = ax.contourf(xx, yy, times_of_arrival, levels, zdir='z', offset=np.min(z_stat), \
    #     cmap='viridis', alpha=1.0)
    # Add a color bar which maps values to colors
    # plt.colorbar(toa_conture, label='Time of arrival (s)')

    ######

    # Constrain the plot to the initial bounds
    # ax.set_xlim(x_min, x_max)
    # ax.set_ylim(y_min, y_max)

    # Set a constant aspect ratio
    #ax.set_aspect('equal', adjustable='datalim')

    # Set a constant and equal aspect ratio
    #ax.set_aspect('equal', adjustable='box-forced')
    #set3DEqualAxes(ax)

    ax.set_xlabel('Latitude')
    ax.set_ylabel('Longitude')
    ax.set_zlabel('Elevation')


    # # plt.savefig(os.path.join(setup.working_directory, setup.fireball_name) + '_3d.png', dpi=300)


    ##########################################################################################################
    return residuals
def getStationList(file_name): 
    """ Reads station .csv file and produces a list containing the station's position and signal time. 
        Accepts files exported from MakeIRISPicks.py. A custom file can be made in the following form:
        *header
        pick_group, station_network, station_code, latitude (deg), longitude (deg), zangle (m), time_of_arrival (Julian Date)

    Arguments:
        file_name: [string] location of the station.csv file

    Returns:
        data: [ndarray] parsed station location and times
    """   
    with open(file_name) as f:

        # Skip the header
        for i in range(1):
            next(f)

        data = []
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '').replace('\t', '')

            # Split the line by the delimiter
            line = line.split(',')

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = (entry.strip())
                if i in [3, 4, 5, 6]:
                    line[i] = float(line[i])
                if i in [3, 4]:
                    line[i] = np.radians(line[i])
            # Add the contents of the line to the data list
            data.append(line)

        return data


def readPoints(output, header=0):
    """ Reads points from a SeismicTrajectory output file so that past graphs can be replotted

    Arguments: 
        output: [string] file containg data that is to be graphed
        header: [integer] number of headers contained in the file

    Returns:
        sup: [ndarray] search data that is to be graphed (global variable)

    """

    with open(output) as f:

        for i in range(header):
            next(f)

        data = np.array([0]*9)
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '')

            # Split the line by the delimiter
            line = line.split(',')

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = float(entry.strip())

            # Add the contents of the line to the data list
            data = np.vstack((data, np.array(line)))

        # First row was all zeroes
        data = np.delete(data, 0, 0)

    particles = data[:, 2:8]
    errors = data[:, 8]

    return sup

if __name__ == "__main__":

    ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description="""
            ~~Seismic Trajectory~~ 
    Search for trajectory of fireball using 
    seismic/infrasound data.

    Denis Vida
    """,
        formatter_class=argparse.RawTextHelpFormatter)

    arg_parser.add_argument('input_file', type=str, help='Path to Supracenter input file.')

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    #################
    # get ini values
    setup = configRead(cml_args.input_file)
    configParse(setup, 'trajectory')

    # Read station file
    station_list = getStationList(setup.station_picks_file)

    # Extract all Julian dates
    jd_list = [entry[6] for entry in station_list]
    # Calculate the arrival times as the time in seconds from the earliest JD
    jd_ref = min(jd_list)
    ref_indx = np.argmin(jd_list)

    try:
        _, _, _, lat0, lon0, elev0, ref_time, pick_time, station_no = station_list[ref_indx]
        lat, lon0, elev0 = setup.lat_centre, setup.lon_centre, 0

    except:
        print("ERROR: data_picks.csv files created previous to Jan. 8, 2019 are lacking a channel tag added. Redownloading the waveform files will likely fix this")

    # Date for weather
    ref_time = jd2Date(jd_ref)
    setup.ref_time = datetime.datetime(*(map(int, ref_time)))

    # Find search area in lat/lon (weather area)

    # Init the constants
    consts = Constants()

    # Convert switches to booleans
    # setup.mode =            setup.mode.lower()   
    # setup.debug =           (setup.debug.lower() == 'true')  
    # setup.enable_winds =    (setup.enable_winds.lower() == 'true')
    # setup.get_data =        (setup.get_data.lower() == 'true')
    # setup.perturb =         (setup.perturb.lower() == 'true')

    sounding = parseWeather(setup, consts)

    # # input is given as latitude/longitude
    # if setup.geomode:


    # Set up search parameters
    p0 = [setup.lat_f, setup.lon_f, setup.t0, setup.v, setup.azim, setup.zangle]

    # Set up perturbed array
    if setup.perturb == True:
        try:

            allTimes = readTimes(setup.arrival_times_file)

        except:
            print("ERROR: unable to find perturbed times. Place a file 'all_pick_times.npy' generated by MakeIRISPicks into the file_name indicated in the SeismicTrajectory.ini file")
        print("Status: Loaded picks from perturbations")
    else:
        allTimes = []

    # If searching
    if setup.run_mode == 'search':                
        estimateSeismicTrajectoryAzimuth(station_list, setup, sounding, p0=p0, azim_range=[setup.azimuth_min, setup.azimuth_max],
            elev_range=[setup.zangle_min, setup.zangle_max], v_fixed=setup.v_fixed, allTimes=allTimes)
    
    # Replot 
    elif setup.run_mode == 'replot':
        dat = readPoints(setup.points_name, header=1)
        p0 = [dat[-1, 0], dat[-1, 1], dat[-1, 2], dat[-1, 3], dat[-1, 4], dat[-1, 5]]
        plotStationsAndTrajectory(station_list, p0, setup, sounding)
 
    elif setup.run_mode == 'manual':
        plotStationsAndTrajectory(station_list, p0, setup, sounding)

    else:
        print('Invalid mode! Use search, replot, or manual.')
        exit()
