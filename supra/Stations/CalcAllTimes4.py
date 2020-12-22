import numpy as np
import os
import copy
from functools import partial
import multiprocessing

from supra.Utils.Classes import Position
from supra.Fireballs.SeismicTrajectory import timeOfArrival
from supra.Atmosphere.Parse import parseWeather
from supra.Supracenter.cyscan2 import cyscan
from supra.Utils.Formatting import *

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

class Times:
    ''' Object to store calculated times in
    '''
    def __init__(self):
        pass

    def __str__(self):
        return 'Times are Calculated'

def checkSkip(bam, prefs):
    """
    Returns TRUE if the time calculations can be skipped
    """

    # Check if the user has requested times to be recalculated
    if prefs.recalc_times:
        if prefs.debug:
            print('DEBUG: Not skipping CalcAllTimes - User override')
        return False

    for stn in bam.stn_list:
        
        # Check if any times
        if not hasattr(stn, 'times'):
            if prefs.debug:
                print('DEBUG: Not skipping CalcAllTimes - No previously calculated times')
            return False

        # Check Ballistic
        if not (len(stn.times.ballistic) >= 1 and prefs.ballistic_en):
            if prefs.debug:
                print('DEBUG: Not skipping CalcAllTimes - No nominal ballistic times')
            return False

        # Check Fragmentations
        if prefs.frag_en and len(stn.times.fragmentation) != len(bam.setup.fragmentation_point):
            if prefs.debug:
                print('DEBUG: Not skipping CalcAllTimes - No nominal fragmentation times')
            return False

        # Check Ballistic Perturb
        try:
            if prefs.pert_en and len(stn.times.ballistic[1]) != prefs.pert_num:
                if prefs.debug:
                    print('DEBUG: Not skipping CalcAllTimes - No perturbed ballistic arrivals')
                return False
        except IndexError:
            if prefs.debug:
                print('DEBUG: Not skipping CalcAllTimes - Cannot read perturbed ballistic arrivals')
            return False

        # Check Fragmentation Perturb
        try:
            if prefs.pert_en and len(stn.times.fragmentation[0][1]) != prefs.pert_num:
                if prefs.debug:
                    print('DEBUG: Not skipping CalcAllTimes - No perturbed fragmentation arrivals')
                return False
        except IndexError:
            if prefs.debug:
                print('DEBUG: Not skipping CalcAllTimes - Cannot read perturbed fragmentation arrivals')
            return False

    return True

def calcAllTimes(bam, prefs):
    ''' Calculates all arrivals to all stations
    '''

    #######################################
    # Check if times need to be calculated
    #######################################

    if checkSkip(bam, prefs):
        print('Skipped Calc Times')
        return bam.stn_list

    ####################
    # Times Calculation
    ####################

    ref_pos = Position(bam.setup.lat_centre, bam.setup.lon_centre, 0)
    no_of_frags = len(bam.setup.fragmentation_point)


    for ii, stn in enumerate(bam.stn_list):
        loadingBar('Calculating Station Times: ', ii+1, len(bam.stn_list))

        if not hasattr(stn, 'times'):
            stn.times = Times()

        stn.times.ballistic = []
        stn.times.fragmentation = []
        stn.metadata.position.pos_loc(ref_pos)

        ################
        # Fragmentation
        ################

        if prefs.frag_en:

            for i, frag in enumerate(bam.setup.fragmentation_point):

                offset = frag.time

                a = []
                
                supra = frag.position
                
                # convert to local coordinates based off of the ref_pos
                supra.pos_loc(ref_pos)

                lats = [supra.lat, stn.metadata.position.lat]
                lons = [supra.lon, stn.metadata.position.lon]
                heights = [supra.elev, stn.metadata.position.elev]

                sounding, perturbations = bam.atmos.getSounding(lats, lons, heights)

                # Travel time of the fragmentation wave
                f_time, frag_azimuth, frag_takeoff, frag_err = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), sounding, \
                    wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                    h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)           

                results = []

                for pert in perturbations:
                    temp = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), pert, \
                        wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                        h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)
                    temp[0] += offset
                    results.append(temp) 

                a.append([f_time + offset, frag_azimuth, frag_takeoff, frag_err])
                a.append(results)
                stn.times.fragmentation.append(a)

        ############
        # Ballistic
        ############

        if prefs.ballistic_en:

            u = np.array([bam.setup.trajectory.vector.x,
                          bam.setup.trajectory.vector.y,
                          bam.setup.trajectory.vector.z])

            angle_off = []
            X = []
            for i in range(len(bam.setup.fragmentation_point)):
                az = stn.times.fragmentation[i][0][1]
                tf = stn.times.fragmentation[i][0][2]

                az = np.radians(az)
                tf = np.radians(180 - tf)
                v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

                angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))
                X.append(bam.setup.fragmentation_point[i].position.elev)
            angle_off = np.array(angle_off)
            try:
                best_indx = np.nanargmin(abs(angle_off - 90))

            except ValueError:
                best_indx = None
                a.append(np.array([np.nan, np.nan, np.nan, np.nan]))
                for pert in perturbations:
                    a.append(np.array([np.nan, np.nan, np.nan, np.nan]))
                stn.times.ballistic.append(a)
                continue

            supra = bam.setup.fragmentation_point[best_indx].position
            supra.pos_loc(ref_pos)
            # Travel time of the fragmentation wave
            f_time, frag_azimuth, frag_takeoff, frag_err = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), sounding, \
                wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)           

            speed = bam.setup.trajectory.v
            distance = supra.pos_distance(bam.setup.trajectory.pos_f)
            timing = distance/speed

            results = []

            for pert in perturbations:
                e, b, c, d = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), pert, \
                    wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                    h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)
                e += timing
                results.append(np.array([e, b, c, d])) 

            a.append([f_time + timing, frag_azimuth, frag_takeoff, frag_err])
            a.append(results)
            stn.times.ballistic.append(a)

    return bam.stn_list