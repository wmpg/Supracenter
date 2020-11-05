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

from supra.Supracenter.anglescan import anglescan
from supra.Utils.AngleConv import loc2Geo, angle2NDE, geo2Loc, angleBetweenVect
from supra.Utils.Classes import Position, Constants, Trajectory, Angle
from supra.Utils.pso import pso
from supra.Utils.Formatting import loadingBar
from supra.Supracenter.cyscan2 import cyscan
from wmpl.Formats.CSSseismic import loadCSSseismicData
from wmpl.Utils.TrajConversions import date2JD, jd2Date, raDec2ECI, geo2Cartesian, cartesian2Geo, raDec2AltAz, eci2RaDec, latLonAlt2ECEF, ecef2ENU, enu2ECEF, ecef2LatLonAlt
from wmpl.Utils.Math import vectMag, vectNorm, rotateVector, meanAngle
from wmpl.Utils.Plotting import Arrow3D, set3DEqualAxes
from wmpl.Utils.PlotMap import GroundMap
from wmpl.Utils.PlotCelestial import CelestialPlot

# def findPoints(setup):    
#     GRID_SPACE = 100
#     MIN_HEIGHT = 0
#     MAX_HEIGHT = 100000

#     u = setup.trajectory.vector.xyz
#     ground_point = setup.trajectory.pos_f.xyz
#     # find top boundary of line given maximum elevation of trajectory
#     if setup.trajectory.pos_i.elev != None:
#         scale = -setup.trajectory.pos_i.elev/u[2]

#     else:
#         scale = -100000/u[2]

#     # define line top boundary
#     top_point = ground_point - scale*u

#     ds = scale / (GRID_SPACE)

#     points = []

#     for i in range(GRID_SPACE + 1):
#         points.append(top_point + i*ds*u)

#     points = np.array(points)

#     offset = np.argmin(np.abs(points[:, 2] - MAX_HEIGHT))
#     bottom_offset = np.argmin(np.abs(points[:, 2] - MIN_HEIGHT))

#     points = np.array(points[offset:(bottom_offset+1)])

#     return points



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


def timeOfArrival(stat_coord, traj, bam, prefs, points, ref_loc=Position(0, 0, 0)):
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
    #x0, y0, t0, v, azim, zangle

    beta = math.sqrt(1 - (prefs.avg_sp_sound/traj.v_avg/1000)**2)

    # Difference from the reference point on the trajectory and the station
    b = stat_coord - traj.findGeo(0).xyz

    u = traj.getTrajVect()

    # Calculate the distance along the trajectory
    
    dt = abs(np.dot(b, -u))

    # Calculate the distance perpendicular to the trajectory
    dp = math.sqrt(abs(vectMag(b)**2 - dt**2))

    R = waveReleasePointWinds(stat_coord, bam, prefs, ref_loc, points, u)


    # if travel:
    #     # travel from trajectory only
    #     ti = R[3]*beta

    # else:
    # v = traj.getVelAtHeight(R[2])
        # if theo:
        #     # Calculate time of arrival
        #     ti = traj.t + dt/v + R[3]*beta
        # else:

            # Calculate time of arrival

    # ti = traj.t - dt/traj.v_avg + R[0]*beta
    ti = R[0]

    ti_pert = []
    for pert_R in R[1]:
        ti_pert.append(traj.t - dt/traj.v_avg + pert_R*beta)

    return ti, ti_pert

def waveReleasePointWindsContour(bam, ref_loc, points, div=37, mode='ballistic'):
    setup = bam.setup
    atmos = bam.atmos
    steps = 90
    alpha = np.linspace(0, 360*((steps-1)/steps), steps)
    alpha = np.radians(alpha)
    theta = setup.trajectory.azimuth.rad
    phi = setup.trajectory.zenith.rad
    tol = 25 #deg
    tol_fact = np.radians(np.linspace(-tol, tol, 10))
    v = [0]*(steps**2)

    results = []

    

    if mode == 'ballistic':
        n_steps = len(tol_fact)*len(points)*steps
        for tt, t in enumerate(tol_fact):
            
            # restriction on beta for initial launch angle to be 90 deg to trajectory
            beta = np.arctan(-1/np.tan(phi)/np.cos(theta-alpha)) + t
            # a = np.sin(phi)*np.cos(theta - alpha)
            # b = np.cos(phi)
            # c = np.cos(t)
            # d = np.sqrt(a**2 + b**2 - c**2)
            # beta = 2*np.arctan((a - d)/(b + c))
            for pp, p in enumerate(points):

                for ii, i in enumerate(range(steps)):

                    step = ii + pp*steps + tt*steps*len(points)
                    loadingBar("Contour Calculation", step, n_steps)

                    v[i] = np.array([np.sin(alpha[i])*np.sin(beta[i]),\
                                     np.cos(alpha[i])*np.sin(beta[i]),\
                                                    -np.cos(beta[i])]) 
                    s = p + p[2]/np.cos(beta[i])*v[i]
                    S = Position(0, 0, 0)
                    P = Position(0, 0, 0)
                    S.x, S.y, S.z = s[0], s[1], s[2]
                    P.x, P.y, P.z = p[0], p[1], p[2]
                    S.pos_geo(ref_loc)
                    P.pos_geo(ref_loc)

                    lats = [P.lat, S.lat]
                    lons = [P.lon, S.lon]
                    elev = [P.elev, S.elev]
                    # z_profile, _ = supra.Supracenter.cyweatherInterp.getWeather(p, s, setup.weather_type, \
                    #      ref_loc, copy.copy(sounding), convert=True)
                    z_profile, _ = atmos.getSounding(lats, lons, elev, spline=100)
                    res = anglescan(p, np.degrees(alpha[i]), np.degrees(beta[i]), z_profile, wind=True, debug=False)

                    # This is the limit in distance from the trajectory (hardcoded)
                    if np.sqrt(res[0]**2 + res[1]**2) <= 400000:
                        results.append(res)

    else:
        n_steps = len(tol_fact)*len(points)*steps

        beta = np.linspace(90 + 0.01, 180, steps)
        beta = np.radians(beta)
        p = points

        for ii, i in enumerate(range(steps)):
            for jj, j in enumerate(range(steps)):
                # print(np.degrees(beta[i]))
                step = jj + ii*steps
                loadingBar("Contour Calculation", step, n_steps)


                v[i*steps + j] = np.array([np.sin(alpha[i])*np.sin(beta[j]),\
                                 np.cos(alpha[i])*np.sin(beta[j]),\
                                                -np.cos(beta[j])]) 
                s = p + p[2]/np.cos(beta[j])*v[i*steps + j]

                S = Position(0, 0, 0)
                P = Position(0, 0, 0)
                S.x, S.y, S.z = s[0], s[1], s[2]
                P.x, P.y, P.z = p[0], p[1], p[2]
                S.pos_geo(ref_loc)
                P.pos_geo(ref_loc)

                lats = [P.lat, S.lat]
                lons = [P.lon, S.lon]
                elev = [P.elev, S.elev]

                z_profile, _ = atmos.getSounding(lats, lons, elev, spline=100)
                res = anglescan(p, np.degrees(alpha[i]), np.degrees(beta[j]), z_profile, wind=True, debug=False)
                
                if np.sqrt(res[0]**2 + res[1]**2) <= 200000:
                    results.append(res)

    return results

def angle2Geo(loc_coord, ref_loc):

    point = Position(0, 0, 0)
    point.x, point.y, point.z = loc_coord[0], loc_coord[1], loc_coord[2]
    point.pos_geo(ref_loc)

    return point

def waveReleasePointWinds(stat_coord, bam, prefs, ref_loc, points, u):
    #azim = (np.pi - azim)%(2*np.pi)
    # Break up the trajectory into points
    
    

    # Trajectory vector
    #u = np.array([-np.cos(azim)*np.sin(zangle), np.sin(azim)*np.sin(zangle), -np.cos(zangle)])

    # Cut down atmospheric profile to the correct heights, and interp


    
    a = (len(points))


    D = np.array(stat_coord)

    cyscan_res = []


    # Compute time of flight residuals for all stations
    for i in range(a):

        S = np.array(points[i])
        traj_point = angle2Geo(points[i], ref_loc)
        stat_point = angle2Geo(stat_coord, ref_loc)

        lats = [traj_point.lat, stat_point.lat]
        lons = [traj_point.lon, stat_point.lon]
        heights = [traj_point.elev, stat_point.elev]

        sounding, perturbations = bam.atmos.getSounding(lats, lons, heights)    

        A = cyscan(S, D, sounding, \
             wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)

        cyscan_res.append(A)
    
    T_nom = getTimes(np.array(cyscan_res), u, a)
    
    T_pert = []
    for p in range(len(perturbations)):

        cyscan_res = []
        for i in range(a):

            S = np.array(points[i])
            traj_point = angle2Geo(points[i], ref_loc)
            stat_point = angle2Geo(stat_coord, ref_loc)

            lats = [traj_point.lat, stat_point.lat]
            lons = [traj_point.lon, stat_point.lon]
            heights = [traj_point.elev, stat_point.elev]

            sounding, perturbations = bam.atmos.getSounding(lats, lons, heights)    

            A_p = cyscan(S, D, perturbations[p], \
             wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)

            cyscan_res.append(A_p)

        T_pert.append(getTimes(np.array(cyscan_res), u, a))

    return T_nom, T_pert

def getTimes(arr, u, a):
    ANGLE_TOL = 25 #degrees
    angle = [999]*a

    cyscan_res = arr
    T = cyscan_res[:, 0]
    az = cyscan_res[:, 1]
    tf = cyscan_res[:, 2]

    az = np.radians(az)
    tf = np.radians(180 - tf)

    mag_u = u/np.sqrt(u.dot(u))

    v = [0]*a

    for ii in range(a):
        v[ii] = np.array([np.sin(az[ii])*np.sin(tf[ii]), np.cos(az[ii])*np.sin(tf[ii]), -np.cos(tf[ii])])

        mag_v = v[ii]/np.sqrt(v[ii].dot(v[ii]))
        angle[ii] = np.absolute(90 - np.degrees(np.arccos(np.dot(mag_u, mag_v))))

    try:
        best_indx = np.nanargmin(angle)
    except:
        return np.array(np.nan)

    if angle[best_indx] > ANGLE_TOL:

        return np.array(np.nan)

    # S = np.array(points[best_indx])
    T = np.array(T[best_indx])

    return T
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

def planeConst(params, station_list, sounding, ref_pos, setup, rest_plane):
    x0, y0, t0, v, azim, zangle = params

    p_1 = np.array(x0, y0, 0)

    vect = getTrajectoryVector(azim, zangle)
    p_2 = -30000*vect

    if rest_plane.checkOnPlane(p_1, tol=100) and rest_plane.checkOnPlane(p_2, tol=100):
        return 1
    else:
        return -1

def trajSearch(params, station_list, sounding, ref_pos, setup, rest_plane):
    
    x0, y0, t0, v, azim, zangle = params

    pos_f = Position(x0, y0, 0)
    pos_f.x = 0
    pos_f.y = 0
    pos_f.z = 0
    pos_f.pos_geo(ref_pos)

    temp_traj = Trajectory(t0, v, zenith=Angle(zangle), azimuth=Angle(azim), pos_f=pos_f)

    points = temp_traj.findPoints(gridspace=100, min_p=np.max((17000, pos_f.elev)), max_p=50000)
    u = temp_traj.vector.xyz

    cost_value = 0

    for stn in station_list:
        t_theo = timeOfArrival(np.array([stn[3], stn[4], stn[5]]), temp_traj, setup, points, sounding=sounding, ref_loc=ref_pos, theo=True)

        t_obs = stn[6]

        cost_value += 2*((1 + (t_theo - t_obs)**2)**0.5 - 1)

        if np.isnan(t_theo): 
            if setup.debug:
                print(np.inf)
            return np.inf   

    if setup.debug:
        print(cost_value)
    return cost_value

# def timeResidualsAzimuth(params, stat_coord_list, arrival_times, setup, sounding, v_fixed=None, \
#         print_residuals=False, pool=[]):
#     """ Cost function for seismic fireball trajectory optimization. The function uses 

#     Arguments:
#         params: [list] Estimated parameters: x0, t0, t0, v, azim, elev.
#         stat_coord_list: [list of ndarrays] A list of station coordinates (x, y, z) in the reference coordinate system.
#         arrival_times: [list] A list of arrival times of the sound wave to the seismic station (in seconds 
#             from some reference time).
#         setup: [Object] Object containing all user-defined parameters
#         sounding: [ndarray] atmospheric profile of the search area

#     Keyword arguments:
#         azim_off: [float] Azimuth around which the given values are centred. If None (default), it is assumed 
#             that the azimuth is calculated +E of due North.
#         v_fixed: [float] Use a fixed velocity. Set to None to ignore

#     """
    
    
#     # Unpack estimated parameters
#     x0, y0, t0, v, azim, zangle = params

#     ref_pos = Position(setup.lat_centre, setup.lon_centre, 0)
#     cost_value = 0
    
#     pos_f = Position(x0, y0, 0)
#     pos_f.pos_geo(ref_pos)

#     temp_traj = Trajectory(t0, v, zenith=Angle(zangle), azimuth=Angle(azim), pos_f=pos_f)

#     points = setup.trajectory.findPoints(gridspace=100, min_p=12000, max_p=50000)
#     u = setup.trajectory.vector.xyz
    
#     # Go through all arrival times
#     for i, (t_obs, stat_coord) in enumerate(zip(arrival_times, stat_coord_list)):
#         ### Calculate the difference between the observed and the prediced arrival times ###
#         ######################################################################################################

#         # Calculate the time of arrival
#         ti = timeOfArrival(stat_coord, x0, y0, t0, v, np.radians(azim), np.radians(zangle), setup, points, u, sounding=sounding, ref_loc=ref_pos, theo=True)

#         # Smooth approximation of l1 (absolute value) loss
#         cost_value += 2*((1 + (t_obs - ti)**2)**0.5 - 1)

#         if np.isnan(ti): 
#             continue     
    
        ######################################################################################################
    # Save points for plotting later

    # Save points and errors for plotting
    # print(cost_value)
    # return cost_value



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

    t0 = time.time()
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
#multiprocessing.cpu_count()
        solutions.append([x, fopt])
        print('Computational  best estimation', fopt)

        print(solutions)

    ### TESTING WITHOUT PERTS
    if setup.perturb and False:

        #allTimes = [perturb_no, station_no, ball/frag, frag_no]
        try:
            perturb_times = allTimes.shape[0]

            p_arrival_times = allTimes[:, station_no, 0, 0] - float(ref_pick_time)

            x_perturb = [0]*perturb_times
            fopt_perturb = [0]*perturb_times

            for i in range(1, perturb_times):
                

                #remove nan
                #p_arrival_times[i] = [j for j in p_arrival_times[i] if j != j]

                # Use PSO for minimization, with perturbed arrival times
                x_perturb[i], fopt_perturb[i] = pso(timeResidualsAzimuth, lower_bounds, upper_bounds, args=(stat_coord_list, \
                    p_arrival_times[i], setup, sounding, v_fixed), maxiter=setup.maxiter, swarmsize=setup.swarmsize, \
                    phip=setup.phip, phig=setup.phig, debug=False, omega=setup.omega, processes=multiprocessing.cpu_count())
                
                print(x_perturb[i], fopt_perturb[i])
                print('Perturbation', i, 'best estimation', fopt_perturb[i])
        except AttributeError:
            x_perturb, fopt_perturb = [], []
            setup.perturb = False
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
    # for i in range(perturb_times):
    #     lat_fin, lon_fin, _ = loc2Geo(ref_pos.lat, ref_pos.lon, ref_pos.elev, [x_perturb[i][0], x_perturb[i][1], 0])
    #     print("Trajectory (Perturbation {:4d}) | Lat {:+10.4f} N Lon {:+10.4f} E t {:5.2f} s v {:7.4f} km/s Azimuth {:6.2f} deg fN Zenith {:5.2f} Error {:10.4f}"\
    #                     .format(i, lat_fin, lon_fin, x_perturb[i][2], x_perturb[i][3]/1000,\
    #                      x_perturb[i][4], x_perturb[i][5], fopt_perturb[i]))   
    #     results.append([lat_fin, lon_fin, x_perturb[i][2], x_perturb[i][3]/1000, x_perturb[i][4], x_perturb[i][5], fopt_perturb[i]])

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

    points = setup.trajectory.findPoints(gridspace=100)
    u = setup.trajectory.vector.xyz
    for stat_coord in stat_coord_list:

        # Calculate time of arrival
        ti = timeOfArrival(stat_coord, setup.trajectory, setup, points, sounding=sounding, ref_loc=Position(lat0, lon0, elev0))
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
