import numpy as np
import os
import copy
from functools import partial
import multiprocessing

from supra.Utils.Classes import Position, Trajectory, Angle
from supra.Fireballs.SeismicTrajectory import timeOfArrival
from supra.Atmosphere.Parse import parseWeather
from supra.Supracenter.cyscan2 import cyscan

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

    return True

def calcAllTimes(bam, prefs):
    ''' Calculates all arrivals to all stations
    '''
    for ii, stn in enumerate(bam.stn_list):
        print(stn.metadata.code, stn.metadata.position)
    exit()
    velocities = [11]
    zes = [5]
    az = 354.67
    t = 0
    pos_i = Position(48.1724466606, 13.0926245672, 50000)

    for vvv in range(len(velocities)):
        for z in range(len(zes)):
            points = []
            file_name = 'C:\\Users\\lmcfd\\Desktop\\Theoretical\\v{:}_ze{:}.csv'.format(str(velocities[vvv]), str(zes[z]))
            with open(file_name, 'w+') as f:
                new_traj = Trajectory(t, velocities[vvv]*1000, zenith=Angle(zes[z]), azimuth=Angle(az), pos_i=pos_i)
                points = new_traj.trajInterp2(div=150,\
                                      min_p=17000,\
                                      max_p=50000)

                az_temp = np.radians(az)
                ze_temp = np.radians(zes[z])

                u = np.array([np.sin(az_temp)*np.sin(ze_temp),
                              np.cos(az_temp)*np.sin(ze_temp),
                             -np.cos(ze_temp)])
                
                print(u)
                   
                # Generate points

                ####################
                # Times Calculation
                ####################

                ref_pos = new_traj.pos_f
                no_of_frags = len(points)

                f.write('Net, Code, Nom time, min time, max time, time range, mean, length, std, nom height, min height, max height, height range, mean, length, std\n')

                for ii, stn in enumerate(bam.stn_list):
                    print("Station {:}/{:} ".format(ii+1, len(bam.stn_list)))

                    stn.times = Times()

                    stn.times.ballistic = []
                    stn.times.fragmentation = []
                    stn.metadata.position.pos_loc(ref_pos)

                    ################
                    # Fragmentation
                    ################

                    if prefs.frag_en:
                        for i, frag in enumerate(points):
                            offset = frag[3]

                            a = []
                            supra = Position(frag[0], frag[1], frag[2])
                            
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

                                results.append(temp) 

                            a.append([f_time, frag_azimuth, frag_takeoff, frag_err])
                            a.append(results)
                            stn.times.fragmentation.append(a)

                    ############
                    # Ballistic
                    ############
                    if prefs.ballistic_en:
                        best_indx = 0
                        az = 0
                        tf = 0
                        v = []
                        h_es = []
                        t_es = []
                        angle_off = []
                        X = []
                        for i in range(len(points)):
                            az = stn.times.fragmentation[i][0][1]
                            tf = stn.times.fragmentation[i][0][2]

                            az = np.radians(az)
                            tf = np.radians(180 - tf)
                            v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])
                            angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))
                            X.append(points[i][2])
                        angle_off = np.array(angle_off)
                        try:
                            best_indx = np.nanargmin(abs(angle_off - 90))

                            if not abs(angle_off[best_indx] - 90) <= 25:
                                best_indx = None
                                raise ValueError


                        except ValueError:
                            best_indx = None

                        if best_indx is not None:
                            h_es.append(points[best_indx][2])
                            t_es.append(stn.times.fragmentation[best_indx][0][0])

                        for pp in range(len(perturbations)):
                            X = []
                            angle_off = []
                            for i in range(len(points)):
                                az = stn.times.fragmentation[i][1][pp][1]
                                tf = stn.times.fragmentation[i][1][pp][2]

                                az = np.radians(az)
                                tf = np.radians(180 - tf)
                                v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])
                                angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))
                                X.append(points[i][2])
                            angle_off = np.array(angle_off)
                            try:
                                best_indx = np.nanargmin(abs(angle_off - 90))

                                if not abs(angle_off[best_indx] - 90) <= 25:
                                    best_indx = None
                                    raise ValueError


                            except ValueError:
                                best_indx = None

                            if best_indx is not None:
                                h_es.append(points[best_indx][2])
                                t_es.append(stn.times.fragmentation[best_indx][1][pp][0])
                        # supra = Position(points[best_indx][0], points[best_indx][1], points[best_indx][2])
                        # supra.pos_loc(ref_pos)
                        # # Travel time of the fragmentation wave
                        # f_time, frag_azimuth, frag_takeoff, frag_err = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), sounding, \
                        #     wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                        #     h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)           

                        # speed = bam.setup.trajectory.v
                        # distance = supra.pos_distance(bam.setup.trajectory.pos_f)
                        # timing = distance/speed

                        # results = []
                        # print('ballistic', f_time)
                        # for pert in perturbations:
                        #     e, b, c, d = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z]), pert, \
                        #         wind=prefs.wind_en, n_theta=prefs.pso_theta, n_phi=prefs.pso_phi,
                        #         h_tol=prefs.pso_min_ang, v_tol=prefs.pso_min_dist)
                        #     e += timing
                        #     results.append(np.array([e, b, c, d])) 

                        # a.append([f_time + timing, frag_azimuth, frag_takeoff, frag_err])
                        # a.append(results)
                        # stn.times.ballistic.append(a)
                    # print csv

                    if len(t_es) > 0 and len(h_es) > 0:
                        f.write('{:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:},\n'.format(stn.metadata.network,\
                                stn.metadata.code, \
                                t_es[0], np.nanmin(t_es), np.nanmax(t_es), np.nanmax(t_es)-np.nanmin(t_es), np.nanmean(t_es), len(t_es), np.std(t_es),\
                                h_es[0], np.nanmin(h_es), np.nanmax(h_es), np.nanmax(h_es)-np.nanmin(h_es), np.nanmean(h_es), len(h_es), np.std(h_es)))
                    else:
                        f.write('{:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:},\n'.format(stn.metadata.network,\
                                stn.metadata.code, \
                                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan,\
                                np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)) 


                # import winsound
                # duration = 1000  # milliseconds
                # freq = 440  # Hz
                # winsound.Beep(freq, duration)
                # winsound.Beep(freq*2, duration)



    return None