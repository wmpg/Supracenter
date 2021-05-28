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

def checkNomBall(bam, prefs):
    """ returns TRUE if the times have been calculated
    """
    for stn in bam.stn_list:
        # Check Ballistic
        if not (len(stn.times.ballistic) >= 1 and prefs.ballistic_en):

            return False

    return True

def checkNomFrag(bam, prefs):
    for stn in bam.stn_list:
        # Check Ballistic

        if not (prefs.frag_en and len(stn.times.fragmentation) == len(bam.setup.fragmentation_point)):

            return False

    return True

def checkPertBall(bam, prefs):
    for stn in bam.stn_list:
        # Check Ballistic

        try:
            if (prefs.pert_en and (len(stn.times.ballistic[0][1]) != prefs.pert_num)):

                return False
        except IndexError:
            return False

    return True

def checkPertFrag(bam, prefs):

    for stn in bam.stn_list:
        # Check Ballistic
        try:
            if (prefs.pert_en and len(stn.times.fragmentation[0][1]) != prefs.pert_num):

                return False
        except IndexError:
            return False

    return True

def checkSkip(bam, prefs):
    """
    Returns TRUE if the time calculations can be skipped
    """

    # Check if the user has requested times to be recalculated
    if prefs.recalc_times:
        if prefs.debug:
            print(printMessage("debug"), 'Not skipping CalcAllTimes - User override')
        return False

    for stn in bam.stn_list:
        # Check if any times
        if not hasattr(stn, 'times'):
            if prefs.debug:
                print(printMessage("debug"), 'Not skipping CalcAllTimes - No previously calculated times')
            return False

    # Check Ballistic
    if not checkNomBall(bam, prefs):
        if prefs.debug:
            print(printMessage("debug"), 'Not skipping CalcAllTimes - No nominal ballistic times')
        return False

    # Check Fragmentations
    if not checkNomFrag(bam, prefs):
        if prefs.debug:
            print(printMessage("debug"), 'Not skipping CalcAllTimes - No nominal fragmentation times')
        return False

    # Check Ballistic Perturb
    if not checkPertBall(bam, prefs):
        if prefs.debug:
            print(printMessage("debug"), 'Not skipping CalcAllTimes - No perturbed ballistic arrivals')
        return False

    # Check Fragmentation Perturb
    if not checkPertFrag(bam, prefs):
        if prefs.debug:
            print(printMessage("debug"), 'Not skipping CalcAllTimes - No perturbed fragmentation arrivals')
        return False


    return True

def printStatus(bam, prefs):

    print("")
    print("CalcAllTimes4.py Status:")
    print("################################")


    print("USER SKIP ", printTrue(prefs.recalc_times))
    print("PROGRAM SKIP ", printTrue(checkSkip(bam, prefs)))
    print("PERTURBATIONS ENABLED ", printTrue(prefs.pert_en))
    print("FRAGMENTATIONS ENABLED ", printTrue(prefs.frag_en))
    print("BALLISTIC ENABLED ", printTrue(prefs.ballistic_en))
    print("NOMINAL FRAGMENTATION CALCULATIONS ", printTrue(checkNomFrag(bam, prefs)))
    print("NOMINAL BALLISTIC CALCULATIONS " , printTrue(checkNomBall(bam, prefs)))
    print("PERTURBATIONS FRAGMENTATION CALCULATIONS ", printTrue(checkPertFrag(bam, prefs)))
    print("PERTURBATIONS BALLISTIC CALCULATIONS " , printTrue(checkPertBall(bam, prefs)))

    print("################################")
    print("")




def calcAllTimes(bam, prefs):
    ''' Calculates all arrivals to all stations
    '''

    #######################################
    # Check if times need to be calculated
    #######################################

    if checkSkip(bam, prefs):
        if prefs.debug:
            printStatus(bam, prefs)
        return bam.stn_list


    ####################
    # Times Calculation
    ####################

    ref_pos = Position(bam.setup.lat_centre, bam.setup.lon_centre, 0)

    if not hasattr(bam.setup, "fragmentation_point"):
        no_of_frags = 0
    else:
        no_of_frags = len(bam.setup.fragmentation_point)

    total_steps = 1 + (prefs.frag_en*no_of_frags*(1 + prefs.pert_en*prefs.pert_num) + prefs.ballistic_en*(1 + prefs.pert_en*prefs.pert_num))*len(bam.stn_list)

    step = 0
    
    

    for ii, stn in enumerate(bam.stn_list):
        

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

                step += 1
                loadingBar('Calculating Station Times: ', step, total_steps)

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

                if perturbations is not None:
                    for pert in perturbations:
                        step += 1
                        loadingBar('Calculating Station Times: ', step, total_steps)
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

            step += 1
            loadingBar('Calculating Station Times: ', step, total_steps)

            a = []
            # define line bottom boundary
            max_height = bam.setup.trajectory.pos_i.elev
            min_height = bam.setup.trajectory.pos_f.elev

            points = bam.setup.trajectory.trajInterp2(div=100, min_p=min_height, max_p=max_height)


            u = np.array([bam.setup.trajectory.vector.x,
                          bam.setup.trajectory.vector.y,
                          bam.setup.trajectory.vector.z])

            angle_off = []
            X = []
            for pt in points:
                
                S = Position(pt[0], pt[1], pt[2])


                lats = [S.lat, stn.metadata.position.lat]
                lons = [S.lon, stn.metadata.position.lon]
                heights = [S.elev, stn.metadata.position.elev]

                S.pos_loc(ref_pos)

                sounding, perturbations = bam.atmos.getSounding(lats, lons, heights)

                # Travel time of the fragmentation wave
                _, az, tf, _ = cyscan(np.array([S.x, S.y, S.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), sounding, \
                    wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                    h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)

                az = np.radians(az)
                tf = np.radians(180 - tf)
                v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

                angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))
                X.append(S.elev)
            angle_off = np.array(angle_off)
            try:
                best_indx = np.nanargmin(abs(angle_off - 90))

            except ValueError:
                best_indx = None
                a.append([np.nan, np.nan, np.nan, np.nan])
                results = []
                if perturbations is not None:
                    for pert in perturbations:
                        results.append(np.array([np.nan, np.nan, np.nan, np.nan]))
                a.append(results)
                stn.times.ballistic.append(a)
                continue

            supra = points[best_indx]
            ref_time = supra[3]
            supra = Position(supra[0], supra[1], supra[2])
            supra.pos_loc(ref_pos)


            lats = [supra.lat, stn.metadata.position.lat]
            lons = [supra.lon, stn.metadata.position.lon]
            heights = [supra.elev, stn.metadata.position.elev]

            sounding, perturbations = bam.atmos.getSounding(lats, lons, heights)
            # Travel time of the fragmentation wave
            f_time, frag_azimuth, frag_takeoff, frag_err = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), sounding, \
                wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)           

            speed = bam.setup.trajectory.v
            distance = supra.pos_distance(bam.setup.trajectory.pos_f)
            timing = distance/speed

            results = []
            if perturbations is not None:
                for pert in perturbations:
                    step += 1
                    loadingBar('Calculating Station Times: ', step, total_steps)
                    e, b, c, d = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), pert, \
                        wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                        h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)
                    e += timing
                    results.append(np.array([e, b, c, d])) 


            a.append([ref_time + f_time , frag_azimuth, frag_takeoff, frag_err])
            a.append(results)

            stn.times.ballistic.append(a)

    step += 1
    loadingBar('Calculating Station Times: ', step, total_steps)

    if prefs.debug:
        printStatus(bam, prefs)

    return bam.stn_list