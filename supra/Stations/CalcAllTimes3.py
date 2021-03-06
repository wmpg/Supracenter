import numpy as np
import os
import copy
from functools import partial
import multiprocessing

from supra.Utils.Classes import Position
from supra.Fireballs.SeismicTrajectory import timeOfArrival
from supra.Atmosphere.Parse import parseWeather
from supra.Supracenter.cyscan2 import cyscan
#from supra.Supracenter.faultscan import cyscan as faultscan

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

class Times:
    def __init__(self):
        pass

    def __str__(self):
        return 'Times are Calculated'

def checkSkip(bam, prefs):
    """
    Returns TRUE if the time calculations can be skipped
    """


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
        if prefs.pert_en and len(stn.times.ballistic[1]) != prefs.pert_num:
            if prefs.debug:
                print('DEBUG: Not skipping CalcAllTimes - No perturbed ballistic arrivals')
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

        if prefs.recalc_times:
            if prefs.debug:
                print('DEBUG: Not skipping CalcAllTimes - User override')
            return False

    return True

def calcAllTimes(bam, prefs):

    if checkSkip(bam, prefs):
        print('Skipped Calc Times')
        return bam.stn_list

    # define line bottom boundary
    if prefs.ballistic_en:
        try:
            points = bam.setup.trajectory.findPoints(gridspace=340, min_p=17000, max_p=50000)
        except AttributeError:
            points = []
    else:
        points = []


    # Ballistic Prediction
    ref_pos = Position(bam.setup.lat_centre, bam.setup.lon_centre, 0)
    no_of_frags = len(bam.setup.fragmentation_point)


    for ii, stn in enumerate(bam.stn_list):
        print("Station {:}/{:} ".format(ii+1, len(bam.stn_list)))

        if not hasattr(stn, 'times'):
            stn.times = Times()

        stn.times.ballistic = []
        stn.times.fragmentation = []
        stn.metadata.position.pos_loc(ref_pos)

        # Ballistic
        if prefs.ballistic_en:

            bam.setup.pos_f.pos_loc(ref_pos)


            # Time to travel from trajectory to station
            stn.times.ballistic = timeOfArrival(stn.metadata.position.xyz, bam.setup.trajectory, \
                                    bam, prefs, points, ref_loc=ref_pos)


        # Fragmentation
        if prefs.frag_en:

            for i, frag in enumerate(bam.setup.fragmentation_point):
                
                a = []
                
                supra = frag.position
                
                # convert to local coordinates based off of the ref_pos
                supra.pos_loc(ref_pos)

                #zProfile = zInterp(stn.position.z, supra.z, zProfile, div=37)
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
                    results.append(cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), pert, \
                        wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                        h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)) 

                a.append([f_time, frag_azimuth, frag_takeoff, frag_err])
                a.append(results)
                stn.times.fragmentation.append(a)

    return bam.stn_list