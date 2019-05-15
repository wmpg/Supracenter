""" Determine the fireball trajectory from seismic data.

Modified method of Pujol et al. (2005).

"""

from __future__ import print_function, division, absolute_import

import sys
import os
import multiprocessing

import copy
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


from supra.Supracenter.cyscan import cyscan
import supra.Supracenter.cyweatherInterp
from supra.Supracenter.netCDFconv import storeHDF, storeNetCDFECMWF, storeNetCDFUKMO, readCustAtm
from supra.Supracenter.fetchECMWF import fetchECMWF
from supra.Supracenter.fetchMERRA import fetchMERRA
from supra.Supracenter.angleConv import loc2Geo, latLon2Local, local2LatLon
from supra.Supracenter.cyzInteg import zInteg
from supra.Supracenter.stationDat import readTimes
from supra.Fireballs.Program import configParse, configRead
from wmpl.Formats.CSSseismic import loadCSSseismicData
from wmpl.Utils.TrajConversions import date2JD, jd2Date, raDec2ECI, geo2Cartesian, cartesian2Geo, raDec2AltAz, eci2RaDec, latLonAlt2ECEF, ecef2ENU, enu2ECEF, ecef2LatLonAlt
from wmpl.Utils.Math import vectMag, vectNorm, rotateVector, meanAngle
from wmpl.Utils.Plotting import Arrow3D, set3DEqualAxes
from wmpl.Utils.PlotMap import GroundMap
from wmpl.Utils.PlotCelestial import CelestialPlot

sup = [0, 0, 0, 0, 0, 0]
errors = [0]


class Constants:

    def __init__(self):

        # Ideal gas constant
        self.R = 8.31432 # J/K mol

        # Kelvins to C
        self.K = 273.15

        # Heat capacity ratio
        self.GAMMA = 1.40

        # Molar mass of the air
        self.M_0 = 0.0289644 #kg/mol


def parseWeather(setup, consts, time=0):

     # Parse weather type
    setup.weather_type = setup.weather_type.lower()
    if setup.weather_type not in ['custom', 'none', 'ecmwf', 'merra', 'ukmo']:
        print('INI ERROR: [Atmospheric] weather_type must be one of: custom, none, ecmwf, merra, ukmo')
        sys.exit()

    # Custom weather data    
    if setup.weather_type == 'custom':

        try:
            sounding = readCustAtm(setup.sounding_file, consts)
        except:
            print('ERROR: Could not read custom atmosphere profile!')
            exit()

        # Check file name
        if len(sounding) < 2:
            print('FILE ERROR: file must have at least two data points!')
            sys.exit()

        # Check file type
        if '.txt' not in setup.sounding_file:
            print("FILE ERROR: custom data set must be a .txt file!")
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
            start_time = (setup.start_datetime.hour + np.round(setup.start_datetime.minute/60) + time)%24
            sounding = storeNetCDFECMWF(setup.sounding_file, setup.lat_centre, setup.lon_centre, consts, start_time=start_time)

        # SeismicTrajectory
        except:

            try:
                start_time = (setup.atm_hour + time)%24
                sounding = storeNetCDFECMWF(setup.sounding_file, setup.x0, setup.y0, consts, start_time=start_time)
            except:
                print("ERROR: Unable to use weather file, or setup.start_datetime/setup.atm_hour is set up incorrectly. Try checking if sounding_file exists")
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

    else:

        # Sample fake weather profile, the absolute minimum that can be passed
        sounding = np.array([[    0.0, setup.v_sound, 0.0, 0.0],
                             [10000.0, setup.v_sound, 0.0, 0.0],
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




def timeOfArrival(stat_coord, x0, y0, t0, v, azim, zangle, setup, sounding=[], ref_loc=[0, 0, 0], travel=False, fast=False):
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
    beta = np.arcsin(setup.v_sound/v)

    # Trajectory vector
    u = np.array([np.sin(azim)*np.sin(zangle), np.cos(azim)*np.sin(zangle), -np.cos(zangle)])

    # Difference from the reference point on the trajectory and the station
    b = stat_coord - np.array([x0, y0, 0])

    # Calculate the distance along the trajectory
    dt = np.abs(np.dot(b, -u))

    # Calculate the distance perpendicular to the trajectory
    dp = np.sqrt(abs(vectMag(b)**2 - dt**2))

    # No winds
    if setup.weather_type == 'none':

        if travel:
            # travel from trajectory only
            ti = (dp*np.cos(beta))/setup.v_sound
        else:
            # Calculate the time of arrival
            ti = t0 - dt/v + (dp*np.cos(beta))/setup.v_sound

    # Winds
    else:

        # Source location
        try: 
            setup.fast_ballistic = setup.fast_ballistic
        except:
            print("CODE WARNING: setup.fast_ballistic not found, setting to false")
            setup.fast_ballistic = False

        if fast:
            # Faster, no-winds solution
            S = waveReleasePoint(stat_coord, x0, y0, t0, v, azim, zangle, setup.v_sound)
        else:
            if setup.fast_ballistic:
                S = waveReleasePoint(stat_coord, x0, y0, t0, v, azim, zangle, setup.v_sound)
            else:
                # Slow, winds solution
                S = waveReleasePointWinds(stat_coord, x0, y0, t0, v, azim, zangle, setup, sounding, ref_loc)

        # Detector location
        D = stat_coord

        # Cut down atmospheric profile to the correct heights, and interp
        zProfile, _ = supra.Supracenter.cyweatherInterp.getWeather(S, D, setup.weather_type, \
                         [ref_loc[0], ref_loc[1], ref_loc[2]], copy.copy(sounding), convert=True)

        # Time of arrival between the points with atmospheric profile
        T, _, _ = cyscan(np.array(S), np.array(D), zProfile, wind=setup.enable_winds)

        if travel:
            # travel from trajectory only
            ti = T*np.cos(beta)

        else:
            # Calculate time of arrival
            ti = t0 -dt/v + T*np.cos(beta)

    return ti

# methods for finding angle between two vectors
# import math
# def dotproduct(v1, v2):
#   return sum((a*b) for a, b in zip(v1, v2))

# def length(v):
#   return math.sqrt(dotproduct(v, v))

# def angle(v1, v2):
#   return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))
def loop(D, points, weather_type, ref_loc, sounding, enable_winds, u, T, az, tf, prop_mag, i):
    
    S = points[i]

    # Cut down atmospheric profile to the correct heights, and interp
    zProfile, _ = supra.Supracenter.cyweatherInterp.getWeather(S, D, weather_type, \
                     [np.degrees(ref_loc[0]), np.degrees(ref_loc[1]), ref_loc[2]], sounding, convert=True)

    # Time of arrival between the points with atmospheric profile
    T[i], az[i], tf[i] = cyscan(np.array(S), np.array(D), zProfile, wind=enable_winds)

    #az[i] = angle2NDE(az[i])
    az[i] = np.radians(az[i])
    tf[i] = np.radians(tf[i])

    #v = np.array([-np.cos(az[i])*np.sin(tf[i]), np.sin(az[i])*np.sin(tf[i]), -np.cos(tf[i])])
    v = np.array([np.sin(az[i])*np.sin(tf[i]), np.cos(az[i])*np.sin(tf[i]), -np.cos(tf[i])])

    prop_mag[i] = np.dot(u, v)

    return prop_mag[i]

def waveReleasePointWinds(stat_coord, x0, y0, t0, v, azim, zangle, setup, sounding, ref_loc):
    #azim = (np.pi - azim)%(2*np.pi)
    # Break up the trajectory into points
    GRID_SPACE = 100

    pool = multiprocessing.Pool(multiprocessing.cpu_count()) 

    # Trajectory vector
    #u = np.array([-np.cos(azim)*np.sin(zangle), np.sin(azim)*np.sin(zangle), -np.cos(zangle)])
    u = np.array([np.sin(azim)*np.sin(zangle), np.cos(azim)*np.sin(zangle), -np.cos(zangle)])

    # define line bottom boundary
    ground_point = (x0, y0, stat_coord[2])
    
    # find top boundary of line given maximum elevation of trajectory
    if setup.elev_i != 0:
        scale = -setup.elev_i*1000/u[2]

    else:
        scale = -100*1000/u[2]

    # define line top boundary
    top_point = ground_point - scale*u

    ds = scale / GRID_SPACE
    points = []

    for i in range(GRID_SPACE + 1):
        points.append(top_point + i*ds*u)

    points = np.array(points)

    D = stat_coord

    T = [0]*(GRID_SPACE + 1)
    az = [0]*(GRID_SPACE + 1)
    tf = [0]*(GRID_SPACE + 1)
    prop_mag = [0]*(GRID_SPACE + 1)

    # Pass data through to multiprocess loop
    iterable = range(GRID_SPACE + 1)

    weather_type = setup.weather_type
    enable_winds = setup.enable_winds

    # Store functions to be used with multiprocessing
    func = partial(loop, D, points, weather_type, ref_loc, sounding, enable_winds, u, T, az, tf, prop_mag)

    # Compute time of flight residuals for all stations
    prop_mag = pool.map(func, iterable)

    pool.close()
    pool.join()

    # Minimize angle to 90 degrees from trajectory
    prop_mag = np.absolute(prop_mag)
    try:
        best_indx = np.nanargmin(prop_mag)
    except:
        print('not good')
        best_indx = 0

    # Return point
    return (top_point + best_indx*ds*u)


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



def timeResidualsAzimuth(params, stat_coord_list, arrival_times, setup, sounding, azim_off=None, v_fixed=None, \
        print_residuals=False):
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

    # Convert the values from km to m
    x0 *= 1000
    y0 *= 1000
    v  *= 1000

    # if v_fixed != '':

    #     # Keep the fireball velocity fixed
    #     v = v_fixed
    try:
        v = float(v_fixed)
    except:
        pass


    # If the azimuth offset is given, recalculate the proper values
    if azim_off is not None:
        azim += azim_off

    # Wrap azimuth and zenith angle to the allowed range
    azim = azim%(2*np.pi)
    zangle = zangle%(2*np.pi)

    cost_value = 0

    # Go through all arrival times
    for i, (t_obs, stat_coord) in enumerate(zip(arrival_times, stat_coord_list)):

        ### Calculate the difference between the observed and the prediced arrival times ###
        ######################################################################################################

        # Calculate the time of arrival
        ti = timeOfArrival(stat_coord, x0, y0, t0, v, azim, zangle, setup, sounding=sounding)

        # Smooth approximation of l1 (absolute value) loss
        z = (t_obs - ti)**2
        cost = 2*((1 + z)**0.5 - 1)

        cost_value += cost

        if print_residuals:
            print("{:>3d}, {:<.3f}".format(i, t_obs - ti))

        ######################################################################################################

    # Save points for plotting later
    global sup 
    global errors

    # Save points and errors for plotting

    # Handling of large errors (such as np.nan values from weather, etc.)
    MAX_ERROR = 1000000

    # check: large err
    if cost_value > MAX_ERROR:
        errors = np.hstack((errors, MAX_ERROR))
        sup = np.vstack((sup, [x0, y0, t0, v, azim, zangle]))
        #cost_value = MAX_ERROR

    # If no errors:
    else:
        errors = np.hstack((errors, cost_value))
        sup = np.vstack((sup, [x0, y0, t0, v, azim, zangle]))

    # Optional print status of searches
    if setup.debug:
        print('Loc: x0 = {:10.2f}, y0 = {:10.2f}, t = {:8.2f}, v = {:4.1f}, azim = {:7.2f}, zangle = {:5.2f}, Err: {:8.2f}'\
            .format(x0, y0, t0, v/1000, np.degrees(azim), np.degrees(zangle), cost_value))\

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
        stat_coord = latLon2Local(lat0, lon0, elev0, lat, lon, elev)

        ######

        stat_coord_list.append(stat_coord)


    return stat_coord_list



def estimateSeismicTrajectoryAzimuth(station_list, setup, sounding, p0=None, azim_range=None, elev_range=None, \
        v_fixed=None, allTimes=[]):
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

    # Extract all Julian dates
    jd_list = [entry[6] for entry in station_list]
    pick_time = [float(entry[7]) for entry in station_list]
    station_no = [int(float(entry[8])) for entry in station_list]

    # Calculate the arrival times as the time in seconds from the earliest JD

    jd_ref = min(jd_list)
    jd_list = np.array(jd_list)
    arrival_times = (jd_list - jd_ref)*86400.0

    # Get the index of the first arrival station
    first_arrival_indx = np.argwhere(jd_list == jd_ref)[0][0]

    ref_pick_time = pick_time[first_arrival_indx]

    # Convert station coordiantes to local coordinates, with the station of the first arrival being the
    # origin of the coordinate system
    stat_coord_list = convertStationCoordinates(station_list, first_arrival_indx)

        
    if p0 is None:

        # Initial parameters:
        p0 = np.zeros(6)

        # Initial parameter
        # Initial point (km)
        # X direction (south positive) from the reference station
        p0[0] = 0
        # Y direction (east positive) from the reference station
        p0[1] = 0

        # Set the time of the wave release to 1 minute (about 20km at 320 m/s)
        p0[2] = -60

        # Set fireball velocity (km/s)
        p0[3] = 20

        # Set initial direction
        # Azim
        p0[4] = np.radians(0)
        # Zangle
        p0[5] = np.radians(45)


    # Find the trajectory by minimization
    # The minimization will probably fail as t0 and the velocity of the fireball are dependant, but the
    # radiant should converge


    # Set azimuth bounds
    # if azim_range is not None:
    #azim_min, azim_max = azim_range

    # Check if the azimuths traverse more than 180 degrees, meaning they are traversing 0/360
    if abs(setup.azimuth_max - setup.azimuth_min) > 180:
        setup.azimuth_min -= 360

    # # Convert the azimuths to +E of due S
    # setup.azimuth_min = np.radians(180 - setup.azimuth_min)
    # setup.azimuth_max = np.radians(180 - setup.azimuth_max)

    # Center the search around the mean azimuth
    azim_avg = (setup.azimuth_min + setup.azimuth_max)/2

    # Calculate the upper and lower azimuth boundaries from the mean
    azim_down = azim_avg - setup.azimuth_max
    azim_up =   azim_avg - setup.azimuth_min

    # Fix order, since min/max of angle angle is counter-intuitive
    if setup.zangle_min < setup.zangle_max:
        setup.zangle_min, setup.zangle_max = setup.zangle_max, setup.zangle_min

    # Set the zangle bounds
    # Calculate zenith angles
    zangle_min = np.min((np.radians(setup.zangle_max), np.radians(setup.zangle_min)))
    zangle_max = np.max((np.radians(setup.zangle_max), np.radians(setup.zangle_min)))

        # Fix order, since min/max of angle angle is counter-intuitive
    if setup.azimuth_min < setup.azimuth_max:
        setup.azimuth_min, setup.azimuth_max = setup.azimuth_max, setup.azimuth_min

    # Set the zangle bounds
    # Calculate zenith angles
    azim_down = np.min((np.radians(setup.azimuth_max), np.radians(setup.azimuth_min)))
    azim_up   = np.max((np.radians(setup.azimuth_max), np.radians(setup.azimuth_min)))

    # Set the bounds for every parameters
    bounds = [
        (setup.x_min, setup.x_max), # X0
        (setup.y_min, setup.y_max), # Y0
        (setup.t_min, setup.t_max), # t0
        (setup.v_min, setup.v_max), # Velocity (km/s)
        (azim_down, azim_up),     # Azimuth
        (zangle_min, zangle_max)  # Zenith angle
        ]

    print('Bounds:', bounds)

    # Extract lower and upper bounds
    lower_bounds = [bound[0] for bound in bounds]
    upper_bounds = [bound[1] for bound in bounds]

    class MiniResults():
        def __init__(self):
            self.x = None

    # Run PSO several times and choose the best solution
    solutions = []

    for i in range(setup.run_times):

        print('Running PSO, run', i)
        print()
        # Use PSO for minimization
        x, fopt = pyswarm.pso(timeResidualsAzimuth, lower_bounds, upper_bounds, args=(stat_coord_list, \
            pick_time, setup, sounding, np.radians(azim_avg), v_fixed), maxiter=setup.maxiter, swarmsize=setup.swarmsize, \
            phip=setup.phip, phig=setup.phig, debug=setup.debug, omega=setup.omega)

        print('Run', i, 'best estimation', fopt)

        solutions.append([x, fopt])

    if setup.perturb == True:

        #allTimes = [perturb_no, station_no, ball/frag, frag_no]
        perturb_times = allTimes.shape[0]

        p_arrival_times = allTimes[:, station_no, 0, 0] - float(ref_pick_time)

        x_perturb = [0]*perturb_times
        fopt_perturb = [0]*perturb_times
        for i in range(perturb_times):
            
            if i == 0:
                print('Computational Picks')
            else:    
                print('Perturbation', i)

            #remove nan
            #p_arrival_times[i] = [j for j in p_arrival_times[i] if j != j]

            # Use PSO for minimization, with perturbed arrival times
            x_perturb[i], fopt_perturb[i] = pyswarm.pso(timeResidualsAzimuth, lower_bounds, upper_bounds, args=(stat_coord_list, \
                p_arrival_times[i], setup, sounding, azim_avg, v_fixed), maxiter=setup.maxiter, swarmsize=setup.swarmsize, \
                phip=setup.phip, phig=setup.phig, debug=setup.debug, omega=setup.omega)
            
            if i == 0:
                print('Computational picks, best estimation', fopt_perturb[i])
            else:
                print('Perturbation', i, 'best estimation', fopt_perturb[i])
    else:
        x_perturb, fopt_perturb = [], []

    # Choose the solution with the smallest residuals
    fopt_array = np.array([fopt for x, fopt in solutions])
    best_indx = np.argmin(fopt_array)

    x, fopt = solutions[best_indx]
    print("X", x)

    res = MiniResults()
    res.x = x

    print("X:", res.x)
    print('Final function value:', fopt)

    # Extract estimated parameters
    x0, y0 = res.x[:2]
    t0 = res.x[2]
    v_est = res.x[3]
    azim, zangle = res.x[4:]

    # Recalculate the proper azimuth if the range was given
    if azim_range:
        azim += azim_avg

    azim = np.radians(azim)
    # Wrap azimuth and zangle to the allowed range
    azim = azim%(2*np.pi)
    zangle = zangle%(np.pi/2)

    lat_fin, lon_fin, _ = local2LatLon((float(station_list[first_arrival_indx][3])), (float(station_list[first_arrival_indx][4])), 0, [x0*1000, y0*1000, 0])

    print('--------------------')
    print('RESULTS:')
    print('x0:', x0, 'km')
    print('y0:', y0, 'km')
    print('lat:', np.degrees(lat_fin))
    print('lon:', np.degrees(lon_fin))
    print('t0:', t0, 's')
    print('v:', v_est, 'm/s')
    print('Azim (+E of due N) : ', (np.degrees(azim))%360)
    print('Elev (from horizon):', 90 - np.degrees(zangle))

    # Print the time residuals per every station
    timeResidualsAzimuth(res.x, stat_coord_list, pick_time, setup, sounding, azim_off=azim_avg, v_fixed=v_fixed, 
        print_residuals=True)

    # Plot the stations and the estimated trajectory
    plotStationsAndTrajectory(station_list, [1000*x0, 1000*y0, t0, 1000*v_est, azim, zangle], setup, sounding, x_perturb=x_perturb)

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

def plotStationsAndTrajectory(station_list, params, setup, sounding, x_perturb=[]):
    """ Plots the stations in the local coordinate system and the fitted trajectory. 
        
    Arguments:
        station_list: [list] contains station names, locations, and travel times of the signal
        params: [list] list of the six parameters of the fireball, x0, y0, t0, v, azim, zangle  
        setup: [Object] Object containing all user-defined parameters
        sounds: [ndarray] atmospheric profile of the search area

    """

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
            x_p[i] = x_perturb[i][0] * 1000
            y_p[i] = x_perturb[i][1] * 1000
            t_p[i] = x_perturb[i][2]
            v_p[i] = x_perturb[i][3] * 1000
            az_p[i] = x_perturb[i][4]
            ze_p[i] = x_perturb[i][5]

    # print("############")
    # print('x', x0, x_p)
    # print('y', y0, y_p)
    # print('t', t0, t_p)
    # print('v', v, v_p)
    # print('azimuth', azim, az_p)
    # print('zenith', zangle, ze_p)
    # print("#############")

    # Extract all Julian dates
    jd_list = [entry[6] for entry in station_list]
    pick_time = [float(entry[7]) for entry in station_list]
    # Calculate the arrival times as the time in seconds from the earliest JD
    jd_ref = min(jd_list)
    ref_indx = np.argmin(jd_list)
    jd_list = np.array(jd_list)
    stat_obs_times_of_arrival = pick_time#(jd_list - jd_ref)*86400.0


    
    # Convert station coordiantes to local coordinates, with the station of the first arrival being the
    # origin of the coordinate system
    stat_coord_list = convertStationCoordinates(station_list, ref_indx)

    # Extract station coordinates
    x_stat, y_stat, z_stat = np.array(stat_coord_list).T

    # Convert station coordinates to km
    x_stat /= 1000
    y_stat /= 1000
    z_stat /= 1000

    x0 /= 1000
    y0 /= 1000

    # Extract coordinates of the reference station
    lat0, lon0, elev0 = station_list[ref_indx][3:6]


    # Calculate the trajectory vector
    u = np.array([np.sin(azim)*np.sin(zangle), np.cos(azim)*np.sin(zangle), -np.cos(zangle)])

    # Calculate modelled times of arrival and points of wave releasefor every station
    stat_model_times_of_arrival = []
    wave_release_points = []

    for stat_coord in stat_coord_list:

        # Calculate time of arrival
        ti = timeOfArrival(stat_coord, 1000*x0, 1000*y0, t0, v, azim, zangle, setup, sounding=sounding, ref_loc=[lat0, lon0, elev0])
        stat_model_times_of_arrival.append(ti)

        # Calculate point of wave release
        traj_point = waveReleasePoint(stat_coord, 1000*x0, 1000*y0, t0, v, azim, zangle, setup.v_sound)/1000
        wave_release_points.append(traj_point)

    stat_model_times_of_arrival = np.array(stat_model_times_of_arrival)
    wave_release_points = np.array(wave_release_points)

    # Get the release points with the highest and the lowest release height
    high_point = np.argmax(wave_release_points[:, 2])
    low_point = np.argmin(wave_release_points[:, 2])

    # Terminal output
    print('x0:', x0, 'km')
    print('y0:', y0, 'km')
    print('t0:', t0, 's')
    print('v:', v, 'm/s')
    print('Azim (+E of due N) : ', np.degrees(azim))
    print('Elev (from horizon):', 90 - np.degrees(zangle))
    print('Wave release:')
    print(wave_release_points[:, 2])
    print(' - top:', wave_release_points[high_point, 2], 'km')
    print(' - bottom:', wave_release_points[low_point, 2], 'km')
    print()
    print('Residuals:')
    print('{:>10s}: {:6s}'.format('Station', 'res (s)'))
    
    sqr_res_acc = 0
    res = []

    # Print the residuals per station
    for stat_entry, toa_obs, toa_model in zip(station_list, stat_obs_times_of_arrival, \
        stat_model_times_of_arrival):

        # Station network code and name
        station_name = stat_entry[1] + '-' + stat_entry[2]

        res.append(toa_obs - toa_model)

        print('{:>10s}: {:+06.2f}s'.format(station_name, toa_obs - toa_model))

        sqr_res_acc += (toa_obs - toa_model)**2


    print('RMS:', np.sqrt(sqr_res_acc/len(station_list)), 's')

    # Corner locations as boundaries for the contour plot
    corner_times = [0, 0, 0, 0]
    corner_coords = [[1000*setup.x_min, 1000*setup.y_min, 0], \
                     [1000*setup.x_min, 1000*setup.y_max, 0], \
                     [1000*setup.x_max, 1000*setup.y_min, 0], \
                     [1000*setup.x_max, 1000*setup.y_max, 0]]

    # The maximum time in the contour plot must be in one of the corners
    for i, coor in enumerate(corner_coords):
        ti = timeOfArrival(coor, 1000*x0, 1000*y0, t0, v, azim, zangle, setup, sounding=sounding)
        corner_times[i] = ti
    
    # Determine the maximum absolute time of arrival (either positive of negative)
    # toa_abs_max = np.max([np.abs(np.min(stat_model_times_of_arrival)), np.max(stat_model_times_of_arrival)])
    toa_abs_max = np.max(np.abs(corner_times))

    # Calcualte absolute observed - calculated residuals
    toa_residuals = np.abs(stat_obs_times_of_arrival - stat_model_times_of_arrival)

    # Determine the maximum absolute residual
    toa_res_max = np.max(toa_residuals)


    ### PLOT 3D ###
    ##########################################################################################################

    # Setup 3D plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')


    # # Plot the stations except the reference
    # station_mask = np.ones(len(x), dtype=bool)
    # station_mask[ref_indx] = 0
    # ax.scatter(x[station_mask], y[station_mask], z[station_mask], c=stat_model_times_of_arrival[station_mask], \
    #     depthshade=0, cmap='viridis', vmin=0.00, vmax=toa_abs_max, edgecolor='k')

    # # Plot the reference station
    # ax.scatter(x[ref_indx], y[ref_indx], z[ref_indx], c='k', zorder=5)

    # Plot the stations (the colors are observed - calculated residuasl)
    stat_scat = ax.scatter(x_stat, y_stat, z_stat, c=toa_residuals, depthshade=0, cmap='inferno_r', \
        edgecolor='0.5', linewidths=1, vmin=0, vmax=toa_res_max)

    plt.colorbar(stat_scat, label='abs(O - C) (s)')

    # Plot the trajectory intersection with the ground
    ax.scatter(x0, y0, 0, c='g')

    # Plot the lowest and highest release points
    wrph_x, wrph_y, wrph_z = wave_release_points[high_point]
    wrpl_x, wrpl_y, wrpl_z = wave_release_points[low_point]
    ax.scatter(wrph_x, wrph_y, wrph_z)
    ax.scatter(wrpl_x, wrpl_y, wrpl_z)
  
    print('Trajectory vector:', u)

    # Get the limits of the plot
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    z_min, z_max = ax.get_zlim()

    x_min = setup.x_min
    x_max = setup.x_max
    y_min = setup.y_min
    y_max = setup.y_max

    # Get the maximum range of axes
    x_range = abs(x_max - x_min)
    y_range = abs(y_max - y_min)
    z_range = abs(z_max - z_min)

    traj_len = 1.5*max([x_range, y_range, z_range])

    # Calculate the beginning of the trajectory
    x_beg = x0 - traj_len*u[0]
    y_beg = y0 - traj_len*u[1]
    z_beg = -traj_len*u[2]

    # Plot the trajectory
    ax.plot([x0, x_beg], [y0, y_beg], [0, z_beg], c='k')

    if len(x_perturb) != 0:
        u_p = [0]*perturb_times

        for i in range(perturb_times):
            u_p[i] = getTrajectoryVector(az_p[i], ze_p[i])

            x_beg_p = x_p[i]/1000 - traj_len*u_p[i][0]
            y_beg_p = y_p[i]/1000 - traj_len*u_p[i][1]
            z_beg_p = -traj_len*u_p[i][2]

            ax.plot([x_p[i]/1000, x_beg_p], [y_p[i]/1000, y_beg_p], [0, z_beg_p], c='b', alpha=0.4)

    # Plot wave release trajectory segment
    ax.plot([wrph_x, wrpl_x], [wrph_y, wrpl_y], [wrph_z, wrpl_z], color='red', linewidth=2)

    ### Plot the boom corridor ###
    img_dim = setup.img_dim
    x_data = np.linspace(x_min, x_max, img_dim)
    y_data = np.linspace(y_min, y_max, img_dim)
    xx, yy = np.meshgrid(x_data, y_data)

    # Make an array of all plane coordinates
    plane_coordinates = np.c_[xx.ravel(), yy.ravel(), np.zeros_like(xx.ravel())]

    times_of_arrival = np.zeros_like(xx.ravel())
    print("Creating contour plot...")

    # Calculate times of arrival for each point on the reference plane
    for i, plane_coords in enumerate(plane_coordinates):

        ti = timeOfArrival(1000*plane_coords, 1000*x0, 1000*y0, t0, v, azim, zangle, setup, sounding=sounding, ref_loc=[lat0, lon0, elev0])

        times_of_arrival[i] = ti

    times_of_arrival = times_of_arrival.reshape(img_dim, img_dim)

    # Determine range and number of contour levels, so they are always centred around 0
    # toa_abs_max = np.max([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
    levels = np.linspace(0.00, toa_abs_max, 50)

    # Plot colorcoded times of arrival on the surface
    toa_conture = ax.contourf(xx, yy, times_of_arrival, levels, zdir='z', offset=np.min(z_stat), \
        cmap='viridis', alpha=1.0)
    # Add a color bar which maps values to colors
    fig.colorbar(toa_conture, label='Time of arrival (s)')

    ######

    # Constrain the plot to the initial bounds
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

    # Set a constant aspect ratio
    #ax.set_aspect('equal', adjustable='datalim')

    # Set a constant and equal aspect ratio
    ax.set_aspect('equal', adjustable='box-forced')
    set3DEqualAxes(ax)

    ax.set_xlabel('X (+ south)')
    ax.set_ylabel('Y (+ east)')
    ax.set_zlabel('Z (+ zenith)')


    plt.savefig(os.path.join(setup.working_directory, setup.fireball_name) + '_3d.png', dpi=300)

    plt.show()


    ### PLOT THE MAP ###
    ##########################################################################################################

    # Calculate the coordinates of the trajectory intersection with the ground
    lat_i, lon_i, elev_i = local2LatLon(lat0, lon0, elev0, 1000*np.array([x0, y0, 0]))

    # Calculate the coordinate of the beginning of the trajectory
    lat_beg, lon_beg, elev_beg = local2LatLon(lat0, lon0, elev0, 1000*np.array([x_beg, y_beg, z_beg]))


    # Calculate coordinates of trajectory release points
    wrph_lat, wrph_lon, wrph_elev = local2LatLon(lat0, lon0, elev0, 1000*np.array([wrph_x, wrph_y, wrph_z]))
    wrpl_lat, wrpl_lon, wrpl_elev = local2LatLon(lat0, lon0, elev0, 1000*np.array([wrpl_x, wrpl_y, wrpl_z]))


    ### Calculate pointing vectors from stations to release points

    wr_points = []
    for wr_point, stat_point in zip(wave_release_points, np.c_[x_stat, y_stat, z_stat]):

        # Calculate vector pointing from station to release point
        wr_point = wr_point - stat_point

        wr_points.append(wr_point)

    
    # Normalize by the longest vector, which will be 1/10 of the trajectory length

    wr_len_max = max([vectMag(wr_point) for wr_point in wr_points])

    wr_vect_x = []
    wr_vect_y = []
    wr_vect_z = []

    for wr_point in wr_points:

        wr_x, wr_y, wr_z = (wr_point/wr_len_max)*traj_len/10

        wr_vect_x.append(wr_x)
        wr_vect_y.append(wr_y)
        wr_vect_z.append(wr_z)

    ###

    # Extract station coordinates
    stat_lat = []
    stat_lon = []
    stat_elev = []

    for entry in station_list:

        lat_t, lon_t, elev_t = entry[3:6]

        stat_lat.append(lat_t)
        stat_lon.append(lon_t)
        stat_elev.append(elev_t)

    stat_lat = np.array(stat_lat)
    stat_lon = np.array(stat_lon)
    stat_elev = np.array(stat_elev)

    # Init the ground map
    m = GroundMap(np.append(stat_lat, lat_i), np.append(stat_lon, lon_i), border_size=20, \
        color_scheme='light')

    # Convert contour local coordinated to geo coordinates
    lat_cont = []
    lon_cont = []
    
    for x_cont, y_cont in zip(xx.ravel(), yy.ravel()):
        
        lat_c, lon_c, _ = local2LatLon(lat0, lon0, elev0, 1000*np.array([x_cont, y_cont, 0]))

        lat_cont.append(lat_c)
        lon_cont.append(lon_c)


    lat_cont = np.array(lat_cont).reshape(img_dim, img_dim)
    lon_cont = np.array(lon_cont).reshape(img_dim, img_dim)

    # Plot the time of arrival contours
    toa_conture = m.m.contourf(np.degrees(lon_cont), np.degrees(lat_cont), times_of_arrival, levels, zorder=3, \
        latlon=True, cmap='viridis')

    # Add a color bar which maps values to colors
    m.m.colorbar(toa_conture, label='Time of arrival (s)')


    # Plot stations
    m.scatter(stat_lat, stat_lon, c=stat_model_times_of_arrival, s=20, marker='o', edgecolor='0.5', \
        linewidths=1, cmap='viridis', vmin=0.00, vmax=toa_abs_max)

    # Plot station names
    for stat_entry in station_list:

        name, s_lat, s_lon = stat_entry[2:5]

        # Convert coordinates to map coordinates
        x, y = m.m(np.degrees(s_lon), np.degrees(s_lat))

        # Plot the text 
        txt = plt.text(x, y, name, horizontalalignment='left', verticalalignment='top', color='k')

        # Add border to text
        txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])

    # Plot intersection with the ground
    m.scatter(lat_i, lon_i, s=10, marker='x', c='g')

    # Plot the trajectory
    m.plot([lat_beg, lat_i], [lon_beg, lon_i], c='g')

    # Plot the wave release segment
    m.plot([wrph_lat, wrpl_lat], [wrph_lon, wrpl_lon], c='red', linewidth=2)


    # Plot wave release directions
    for i in range(len(x_stat)):
        
        wr_lat_i, wr_lon_i, wr_elev_i = local2LatLon(lat0, lon0, elev0, 1000*np.array([x_stat[i] + wr_vect_x[i], y_stat[i] + wr_vect_y[i], z_stat[i] + wr_vect_z[i]]))
        #wr_lat_i, wr_lon_i, wr_elev_i = local2LatLon(lat0, lon0, elev0, 1000*np.array([wrp_x[i], wrp_y[i], wrp_z[i]]))

        m.plot([stat_lat[i], wr_lat_i], [stat_lon[i], wr_lon_i], c='k', linewidth=1.0)


    plt.tight_layout()

    plt.savefig(os.path.join(setup.working_directory, setup.fireball_name) + '_map.png', dpi=300)

    plt.show()

    global sup
    global errors
    # Potential Supracenter locations
    sup = np.delete(sup, 0, 0)

    # Error in each potential Supracenter
    errors = np.delete(errors, 0)

    sup = np.delete(sup, -1, 0)

    if setup.run_mode != 'manual':
        errors = np.resize(errors, errors.size - 1)
        # Create search points file (used to replot previously found data)
        if setup.run_mode == 'search':
            
            print('Writing data file...')

            with open(os.path.join(setup.working_directory, setup.fireball_name) + '_plotted_points.csv', 'w') as f:
                f.write('lat, lon, x0, y0, t0, v, az, ze, err\n')

                for i in range(len(sup)):

                    alat, alon, _ = local2LatLon(lat0, lon0, elev0, np.array([sup[i, 0], sup[i, 1], 0]))
                    f.write('{:9.2f}, {:9.2f}, {:9.4f}, {:9.4f}, {:6.2f}, {:6.2f}, {:6.2f}, {:6.2f}, {:6.2f}\n'\
                        .format(np.degrees(alat), np.degrees(alon), sup[i, 0], sup[i, 1], sup[i, 2], sup[i, 3], sup[i, 4], sup[i, 5], errors[i]))
                alat, alon, _ = local2LatLon(lat0, lon0, elev0, np.array([x0, y0, 0]))

                zangle = np.radians(90 - np.degrees(zangle))
                azim = np.radians((np.degrees(azim))%360)

                f.write('{:9.2f}, {:9.2f}, {:9.4f}, {:9.4f}, {:6.2f}, {:6.2f}, {:6.4f}, {:6.4f}, {:6.2f}\n'\
                        .format(np.degrees(alat), np.degrees(alon), x0*1000, y0*1000, t0, v, azim, zangle, 0))

            print('Done writing data file...')

        # # Remove points above a certain threshold
        # for i in range(1):

        #     a = []
        #     std_error = np.std(errors)

        #     # Threshold
        #     lim = np.mean(errors) - 0*std_error

        #     for i in range(len(errors)):
        #         if errors[i] >= lim:
        #             a.append(i)

        #     errors = np.delete(errors, (a), axis=0)
        #     sup = np.delete(sup, (a), axis=0)


        # Create uncertainty errors
        # - the best points
        unc_errors = errors.copy()
        unc_sup = sup.copy()

        # std = np.std(unc_errors)

        # #Only keep the best errors for uncertainty
        # lim = np.mean(unc_errors) - 0.75*std
        # lim = np.mean(unc_errors) - np.sqrt(2*np.mean(unc_errors))
        
        # for i in range(len(unc_errors)):
        #     if unc_errors[i] >= lim:
        #         a.append(i)

        # unc_errors = np.delete(unc_errors, (a), axis=0)
        # unc_sup = np.delete(unc_sup, (a), axis=0)


        if len(unc_errors) == 0:
            unc_errors = errors
            unc_sup = sup

        else:
            try:
                # Set up axes
                fig = plt.figure(figsize=plt.figaspect(0.5))
                fig.set_size_inches(20.9, 11.7)
                ax5 = fig.add_subplot(111)

                bins = np.linspace(0, max(errors), max(errors)/2)

                if len(bins) <= 1:
                    bins = np.linspace(0, max(errors) + 1, max(errors)/2)

                # plot histrogram of errors showing both unfiltered (after removing large errors) and filtered (uncertainties)
                ax5.hist(    errors, bins, alpha=0.5, label='Unfiltered')
                ax5.hist(unc_errors, bins, alpha=0.5, label='Filtered')
                ax5.set_title('Filtered Error Data')
                ax5.set_xlabel('Error (s)')
                ax5.set_ylabel('Frequency')
                ax5.legend(loc='upper right')
                plt.show()
                plt.clf()
            except:
                print('ERROR: Cannot display histogram of errors!')

        # Set up axes
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(20.9, 11.7)


        ### FOUR-WINDOW PLOT ###
        #############################################################################################################

        # x, y positions
        ax1 = fig.add_subplot(222)

        # potential Supracenter points with color based on error
        sco = ax1.scatter(sup[:, 1]/1000, sup[:, 0]/1000, c=errors, cmap='inferno_r', s=4)
        #ax1.scatter(unc_sup[:, 1]/1000, unc_sup[:, 0]/1000, c='g', s=4)
        ax1.scatter(y0, x0, c = 'c', marker='x', s=50, linewidth=3)

        ax1.set_title("Position of Origin")
        ax1.set_ylabel('X (+ south) (km)')
        ax1.set_xlabel('Y (+ east) (km)')

        sc = ax1.scatter(y_stat, x_stat, c=toa_residuals, s=24, marker='^', cmap='viridis_r')
        ax1.set_ylim([setup.x_min, setup.x_max])
        ax1.set_xlim([setup.y_min, setup.y_max])

        # Plot station names
        for stat_entry in station_list:

            name, s_lat, s_lon = stat_entry[2:5]

            # Convert coordinates to map coordinates
            x, y, _ = latLon2Local(lat0, lon0, elev0, s_lat, s_lon, 0)

            # Plot the text 
            txt = ax1.text(y/1000, x/1000, name, horizontalalignment='left', verticalalignment='top', color='k', size=6)

            # Add border to text
            txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])

        # v, t plot
        ax2 = fig.add_subplot(221)
        ax2.scatter(sup[:, 2], sup[:, 3]/1000, c=errors, cmap='inferno_r', s=4)
        #ax2.scatter(unc_sup[:, 2], unc_sup[:, 3]/1000, c='g', s=4)
        ax2.scatter(t0, v/1000, c = 'c', marker='x', s=50, linewidth=3)

        ax2.set_title("Velocity vs Time")
        ax2.set_xlabel("Time (s)")
        ax2.set_ylabel("Velocity (km/s)")
        ax2.set_xlim([setup.t_min, setup.t_max])
        ax2.set_ylim([setup.v_min, setup.v_max])

        a = plt.colorbar(sc, ax=ax1)
        a.set_label("Station Residual (s)")

        # angle plot
        ra, dec = np.radians(sup[:, 4]), np.radians(90 - sup[:, 5])
        ax3 = fig.add_subplot(223)
        ax3 = CelestialPlot(ra, dec, projection='stere', bgcolor='w')

        ax3.scatter(ra, dec, c=errors, cmap='inferno_r', s=4)
        #ax3.scatter(np.radians(unc_sup[:, 4]), np.radians(unc_sup[:, 5]), c='g', s=4)

        ax3.scatter(azim, np.pi/2 - zangle, c = 'c', marker='x', s=50, linewidth=3)

        plt.title('zangle (from Horizontal) and Azimuth (+E due N) of Trajectory')

        # 2D plot
        ax4 = fig.add_subplot(224)
        ax4 = GroundMap(np.append(stat_lat, lat_i), np.append(stat_lon, lon_i), border_size=20, \
            color_scheme='light')
         # Add a color bar which maps values to colors
        ax4.m.colorbar(toa_conture, label='Time of arrival (s)')

        # Plot intersection with the ground
        ax4.scatter(lat_i, lon_i, s=10, marker='x', c='g')

        # Plot the trajectory
        ax4.plot([lat_beg, lat_i], [lon_beg, lon_i], c='g')

        # Plot the wave release segment
        ax4.plot([wrph_lat, wrpl_lat], [wrph_lon, wrpl_lon], c='red', linewidth=2)

        # Plot wave release directions
        for i in range(len(x_stat)):
            
            wr_lat_i, wr_lon_i, wr_elev_i = local2LatLon(lat0, lon0, elev0, 1000*np.array([x_stat[i] + wr_vect_x[i], y_stat[i] + wr_vect_y[i], z_stat[i] + wr_vect_z[i]]))
            #wr_lat_i, wr_lon_i, wr_elev_i = local2LatLon(lat0, lon0, elev0, 1000*np.array([wrp_x[i], wrp_y[i], wrp_z[i]]))

            ax4.plot([stat_lat[i], wr_lat_i], [stat_lon[i], wr_lon_i], c='k', linewidth=1.0)
        
        # Plot stations
        ax4.scatter(stat_lat, stat_lon, c=stat_model_times_of_arrival, s=20, marker='o', edgecolor='0.5', \
            linewidths=1, cmap='viridis', vmin=0.00, vmax=toa_abs_max)

        for stat_entry in station_list:

            name, s_lat, s_lon = stat_entry[2:5]

            # Convert coordinates to map coordinates
            x, y = ax4.m(np.degrees(s_lon), np.degrees(s_lat))

            # Plot the text 
            txt = plt.text(x, y, name, horizontalalignment='left', verticalalignment='top', color='k', size=6)

            # Add border to text
            txt.set_path_effects([PathEffects.withStroke(linewidth=1, foreground='w')])

        # Plot the time of arrival contours
        toa_conture = ax4.m.contourf(np.degrees(lon_cont), np.degrees(lat_cont), times_of_arrival, levels, zorder=3, \
            latlon=True, cmap='viridis')

        # Add a color bar which maps values to colors
        ax4.m.colorbar(toa_conture, label='Time of arrival (s)')

        lat_dot =     [0]*len(sup)
        lon_dot =     [0]*len(sup)
        # unc_lat_dot = [0]*len(unc_sup)
        # unc_lon_dot = [0]*len(unc_sup)
        for i in range(len(sup)):
            lat_dot[i], lon_dot[i], _ = local2LatLon(lat0, lon0, elev0, [sup[i, 0], sup[i, 1], 0])
        # for i in range(len(unc_sup)):
        #     unc_lat_dot[i], unc_lon_dot[i], _ = local2LatLon(lat0, lon0, elev0, [unc_sup[i, 0], unc_sup[i, 1], 0])

        lat_x, lon_x, _ = local2LatLon(lat0, lon0, elev0, [x0*1000, y0*1000, 0])

        ax4.scatter(lat_dot, lon_dot, c=errors, cmap='inferno_r', s=4)
        #ax4.scatter(unc_lat_dot, unc_lon_dot, c='g', s=4)

        ax4.scatter(lat_x, lon_x, c = 'c', marker='x', s=50, linewidth=3)

        fig.subplots_adjust(right=0.8)

        cb_ax = fig.add_axes([0.83, 0.1, 0.02, 0.8])
        fig.colorbar(sco, cax=cb_ax)

        plt.savefig(os.path.join(setup.working_directory, setup.fireball_name) + '_scatter.png', dpi=300)

        plt.title('Error of Solution (s)')

        plt.show()
    #####################
    # OUTPUT TEXT FILES #
    #####################

    wave_release_points_geo = [0]*len(station_list)

    for i in range(len(station_list)):
        wave_release_points_geo[i] = local2LatLon(lat0, lon0, elev0, [wave_release_points[i][0]*1000, wave_release_points[i][1]*1000, wave_release_points[i][2]])

    # Results .txt
    ##############
    with open(os.path.join(setup.working_directory, setup.fireball_name) + '_results.txt', 'w') as f:

        # Using PSO
        if setup.run_mode == 'search':
            f.write('Mode: Particle Swarm Search\n')
        elif setup.run_mode == 'replot':
            f.write('Mode: Replot Data\n')
        elif setup.run_mode == 'manual':
            f.write('Mode: Manual Search\n')
        else:
            f.write('Mode: ERROR\n')

        # Using restricted velocity
        if setup.v_fixed == '':
            f.write('Unrestricted Fireball Velocity\n')
        else:
            f.write('Fireball Velocity Restricted to {:5f}\n'.format(setup.v_fixed))

        # Using Winds
        if setup.enable_winds:
            f.write('Winds Enabled\n')
        else:
            f.write('Winds Disabled\n')

        f.write('x0: {:8.3f} km \n'.format(x0))
        f.write('y0: {:8.3f} km \n'.format(y0))
        f.write('t0: {:8.3f} s  \n'.format(t0))
        f.write('v: {:8.3f} m/s\n'.format(v ))
        f.write('Azimuth (+E of due N): {:6.2f}\n'.format(np.degrees(azim)))
        f.write('zangle (from horizon): {:6.2f}\n'.format(90 - np.degrees(zangle)))
        f.write('==========================================================================================================\n')
        f.write('  Station Name       Station Location          Residual       Travel Time         Wave Release Point      \n')
        f.write('                 (lat( N),lon( E),elev(m))       (s)              (s)          (lat( N),lon( E),elev(m))  \n')
        f.write('==========================================================================================================\n')
        for i in range(len(station_list)):
            f.write('      {:}        {:7.4f} {:7.4f} {:5.0f}        {:+7.3f}          {:6.2f}        {:7.3f} {:7.3f} {:7.3f}\n'\
                .format(station_list[i][2], np.degrees(stat_lat[i]), np.degrees(stat_lon[i]), stat_elev[i], res[i], stat_model_times_of_arrival[i], \
                    np.degrees(wave_release_points_geo[i][0]), np.degrees(wave_release_points_geo[i][1]), wave_release_points_geo[i][2]))              
        f.close()

    # Sounding Output .txt
    ######################

    # Compatible with darkflight
    # Produces weather profile at the low point of the trajectory

    if setup.weather_type != 'none':

        # Create weather profile at end point for use with darkflight
        # Use lowest point
        sounding, _ = supra.Supracenter.cyweatherInterp.getWeather([wrpl_x*1000, wrpl_y*1000, wrpl_z*1000], [wrpl_x*1000, wrpl_y*1000, 0], setup.weather_type, [x0, y0, 0], sounding)
        n, m = sounding.shape

        # Add model pressure (see the U. S. Standard Atmosphere Model of 1976)
        pressure_list = addPressure(sounding)

        # with open(setup.file_name + 'darkflight_sounding.txt', 'w') as f:

        #     for i in range(n):
        #         f.write(' {:8.2f} {:7.3f} {:7.3f} {:6.2f} {:6.2f}\n'.format(sounding[i, 0], sounding[i, 1], \
        #             sounding[i, 2], (sounding[i, 3]), float(pressure_list[i])))
        #     f.close()

    # If no weather profile is given, an isotropic model is given
    else:
        with open(os.path.join(setup.working_directory, setup.fireball_name) + 'darkflight_sounding.txt', 'w') as f:

            # fake height profile
            n = [0, 10000, 20000, 30000, 40000, 50000]

            # model pressures 
            p = [1000, 259, 49, 11, 3, 1]

            for i in range(len(n)):
                f.write(' {:8.2f} {:7.3f} {:7.3f} {:6.2f} {:6.2f}\n'.format(n[i], 0, 0, 0, p[i]))
            f.close()

    # Produce input file for darkflight
    ####################################

    with open(os.path.join(setup.working_directory, setup.fireball_name) + '_darkflight_input.ini', 'w') as f:
        pos = local2LatLon(lat0, lon0, elev0, [wrpl_x*1000, wrpl_y*1000, wrpl_z*1000])
        f.write('[General]\n')
        f.write('error=True\n')
        f.write('output={:}\n'.format(os.path.join(setup.working_directory, setup.fireball_name) + 'darkflight_trajectory_path'))
        f.write('atm={:}\n'.format(os.path.join(setup.working_directory, setup.fireball_name) + 'darkflight_sounding.txt'))
        f.write('\n')
        f.write('[Variables]\n')

        # Use final found values for darkflight, these values represent the end of the trajectory, not any point during the trajectory
        # and may not represent the actual darkflight variables
        f.write('lat = {:f}\nlon = {:f}\nz = {:f}\nv = {:f}\naz = {:f}\nze = {:f}\n'\
            .format(np.degrees(pos[0]), np.degrees(pos[1]), pos[2]/1000, v/1000, \
                (180 - np.degrees(azim))%360, (90 - np.degrees(zangle))))
        
        # Setup blank ini file
        # Blank values will be set as default
        # These can be further changed by the user
        f.write('mass_min = \n')
        f.write('mass_max = \n')
        f.write('tan = \n')
        f.write('sig = \n')
        f.write('den =\n')
        f.write('end = \n')
        f.write('shape = \n')
        f.write('dlat = \n')
        f.write('dlon = \n')
        f.write('dh = \n')
        f.write('dv = \n')
        f.write('daz = \n')
        f.write('dze = \n')
        f.write('max_iter = \n')
        f.write('drag = \n')
        f.write('\n')
        f.write('# This file was produced by SeismicTrajectory.py and should only be used as an estimation! Weather and variables may be inaccurate!')
        f.write('# The velocity produced in this file is based off of the trajectory solution, and is a very high estimate')
        f.close()

    ##########################################################################################################
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

    global sup 
    global errors

    sup = data[:, 2:8]
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

    except:
        print("ERROR: data_picks.csv files created previous to Jan. 8, 2019 are lacking a channel tag added. Redownloading the waveform files will likely fix this")

    # Date for weather
    ref_time = jd2Date(jd_ref)
    setup.ref_time = datetime.datetime(*(map(int, ref_time)))

    # Find search area in lat/lon (weather area)
    setup.search_area[0] = min(np.degrees(local2LatLon(lat0, lon0, elev0, [setup.x_min*1000, 0, 0])[0]), np.degrees(local2LatLon(lat0, lon0, elev0, [setup.x_max*1000, 0, 0])[0]))
    setup.search_area[1] = max(np.degrees(local2LatLon(lat0, lon0, elev0, [setup.x_min*1000, 0, 0])[0]), np.degrees(local2LatLon(lat0, lon0, elev0, [setup.x_max*1000, 0, 0])[0]))
    setup.search_area[2] = min(np.degrees(local2LatLon(lat0, lon0, elev0, [0, setup.y_min*1000, 0])[1]), np.degrees(local2LatLon(lat0, lon0, elev0, [0, setup.y_max*1000, 0])[1]))
    setup.search_area[3] = max(np.degrees(local2LatLon(lat0, lon0, elev0, [0, setup.y_min*1000, 0])[1]), np.degrees(local2LatLon(lat0, lon0, elev0, [0, setup.y_max*1000, 0])[1]))

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
    setup.lat_f, setup.lon_f, _ = latLon2Local(lat0, lon0, elev0, np.radians(setup.lat_f), np.radians(setup.lon_f), 0)   

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
        p0[3] *= 1000
        plotStationsAndTrajectory(station_list, p0, setup, sounding)

    else:
        print('Invalid mode! Use search, replot, or manual.')
        exit()

    #TEMPORARY EXIT FOR IF INPUTS ARE NOT AVAILABLE
    exit()

    ##########################################################################################################

    ### INPUTS ###
    ##########################################################################################################


    # DATA PATHS
    dir_paths = [
        "/local4/infrasound/Infrasound/Fireball/15-Sep-2007/seismic/UBh2"
        #"/local4/infrasound/Infrasound/Fireball/15-Sep-2011-SW_USA/is56"
        #"/local4/infrasound/Infrasound/Fireball/15-Sep-2011-SW_USA/is57"
        #"/local4/infrasound/Infrasound/Fireball/15-Sep-2007/seismic/LPAZ"
    ]

    site_files = [
        "UBh2.site"
        #"is56.site"
        #"is57.site"
        #"impact.site"
        ]

    wfdisc_files = [
        "UBh2.wfdisc"
        #"is56.wfdisc"
        #"is57.wfdisc"
        #"impact.wfdisc"
        ]

    # Average speed of sound in the atmosphere
    v_sound = 320 # m/s


    ##########################################################################################################

    seismic_data = []

    # Load seismic data from given files
    for dir_path, site_file, wfdisc_file in zip(dir_paths, site_files, wfdisc_files):

        # Load the seismic data from individual file
        file_data = loadCSSseismicData(dir_path, site_file, wfdisc_file)

        # Add all entries to the global list
        for entry in file_data:
            seismic_data.append(entry)


    # Merge all data from the same channels
    merged_data = mergeChannels(seismic_data)

    # Determine the earliest time from all beginning times
    ref_time = min([w.begin_time for _, merged_w, _, _ in merged_data for w in merged_w])

    # Setup the plotting
    f, axes = plt.subplots(nrows=len(merged_data), ncols=1, sharex=True)

    # Go through data from all channels
    for i, entry in enumerate(merged_data):

        # Select the current axis for plotting
        ax = axes[i]

        # Unpack the channels
        merged_sites, merged_wfdisc, merged_time, merged_waveform = entry

        # Go through all individual measurements
        for site, w, time_data, waveform_data in zip(merged_sites, merged_wfdisc, merged_time, merged_waveform):

            # Calculate the difference from the reference time
            t_diff = (w.begin_time - ref_time).total_seconds()

            # Offset the time data to be in accordance with the reference time
            time_data += t_diff

            # Plot the seismic data
            ax.plot(time_data/60, waveform_data, zorder=3, linewidth=0.5)

        # Add the station label
        ax.text(ax.get_xlim()[0], ax.get_ylim()[1], w.sta + " " + w.chan, va='top')

        ax.grid(color='0.9')


    plt.xlabel('Time (min)')

    plt.subplots_adjust(hspace=0)
    plt.show()

