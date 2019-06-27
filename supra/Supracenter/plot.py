""" Outputs Data through terminal, .txt files, and graphs"""

import os
import copy
import datetime
import time

import numpy as np
import matplotlib.pyplot as plt

from supra.Utils.AngleConv import loc2Geo, geo2Loc, angle2NDE
from supra.Utils.Classes import Constants
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.cyscan import cyscan
from supra.Supracenter.slowscan import slowscan

def pointsList(sounding, points):
    """ Helper function: Divides list sounding into 'length of points' divisions and expands each point into 
        every value of each sounding division. In short, expands points to have a length of sounding, while 
        keeping the same points values.

        Arguments:
            sounding: [ndarray] final interpolated weather profile
            points: [ndarray] list of points used in the interpolation

        Returns:
            locs: [ndarray] array of locations expanded in length to match the length of sounding

    """

    # Find length of arrays
    s = len(sounding)
    p = len(points)

    # Initialize arrays
    locs = np.array([0, 0])

    for i in range(s):

        # Add the correct value to the list
        locs = np.vstack((locs, points[int(p*i/s)]))

    # First row was all zeroes
    locs = np.delete(locs, 0, 0)

    return locs


def outputWeather(n_stations, x_opt, stns, setup, ref_pos, dataset, output_name, s_name, kotc, w):
    """ Function to calculate the results of the optimal Supracenter and export interpolated weather data
        from each station.

    Arguments:
        n_stations: [int] number of stations
        x_opt: [list] lat, lon, and height of the optimal supracenter
        xstn: [list] lat, lon, and height of all the stations
        setup: [Object] object storing user-defined setup
        consts: [Object] object storing physical constants
        ref_pos: [list] mean position of the stations to act as the origin for the local coordinate system
        dataset: [ndarray] atmospheric profile of search area
        output_name: [string] folder name to save output files in
        s_name: [list] name of all the stations
        kotc: [list] user-defined occurence time [HH, MM, SS]
        tobs: [list] station arrival times to each station
        w: [list] weights of each station


    Returns:
        time3D: [list] calculated travel time to each station
        az: [list] initial azimuth angle to each station
        tf: [list] initial takeoff angle to each station
        r: [list] residuals to each station            
    
    """
    consts = Constants()
    # Station Times
    tobs = stns[0:n_stations, 4]

    # Station Location
    xstn = stns[0:n_stations, 0:3]

    # Travel time to each station
    time3D = np.zeros(n_stations)

    # Initial azimuths angles of each station
    az = np.zeros(n_stations)

    # Initial takeoff angles of each station
    tf = np.zeros(n_stations)

    # difference in theoretical and simulated travel times
    sotc = np.zeros_like(tobs)

    #difference in theoretical and simulated travel times
    r = np.zeros_like(tobs)

    trace = []

    # Find parameters of optimal location
    print('Exporting weather profiles...')
    for j in range(n_stations):
        #print(np.array(x_opt), np.array(xstn[j, :]))
        # get interpolated weather profile
        sounding, points = getWeather(x_opt.xyz, np.array(xstn[j, :]), setup.weather_type, ref_pos, copy.copy(dataset), convert=True)

        # Rotate winds to match with coordinate system
        #sounding[:, 3] = np.radians(angle2NDE(np.degrees(sounding[:, 3])))

        # If weather interp was used
        if len(points) != 0:
            points = pointsList(sounding, points)

        # Export the interpolated weather profiles from the optimal Supracenter for each station
        with open(os.path.join(output_name, s_name[j] + '_sounding.txt'), 'w') as f:

            # With winds
            if setup.enable_winds == True:
                if setup.weather_type == 'custom' or setup.weather_type == 'none':
                    f.write('| Height (m) | Temp (K) | soundSpd (m/s) | wdSpd (m/s) | wdDir (deg fN) |\n')
                else:
                    f.write('| Latitude (deg N) | Longitude (deg E) | Height (m) | Temp (K) | soundSpd (m/s) | wdSpd (m/s) | wdDir (deg fN) |\n')
            
            # No winds
            else:
                if setup.weather_type == 'custom' or setup.weather_type == 'none':
                    f.write('| Height (m) | Temp (K) | soundSpd (m/s) |\n')
                else:
                    f.write('| Latitude (deg N) | Longitude (deg E) | Height (m) | Temp (K) | soundSpd (m/s) |\n')

            for ii in range(len(sounding)):

                if setup.enable_winds == True:
                    if setup.weather_type == 'custom' or setup.weather_type == 'none':
                        f.write('|  {:8.2f}  |  {:7.3f} |    {:6.4f}    |   {:7.3f}   |     {:6.2f}     |\n'\
                            .format(sounding[ii, 0], sounding[ii, 1]**2*consts.M_0/consts.GAMMA/consts.R, \
                                sounding[ii, 1], sounding[ii, 2], sounding[ii, 3]))
                    else:
                        f.write('|       {:7.4f}    |      {:7.4f}    |  {:8.2f}  |  {:7.3f} |    {:6.4f}    |   {:7.3f}   |     {:6.2f}     |\n'\
                            .format(points[ii][0], points[ii][1], sounding[ii, 0], sounding[ii, 1]**2*consts.M_0/consts.GAMMA/consts.R, \
                                sounding[ii, 1], sounding[ii, 2], sounding[ii, 3]))
                else:
                    if setup.weather_type == 'custom' or setup.weather_type == 'none':
                        f.write('|  {:8.2f}  |  {:7.3f} |    {:6.4f}    |\n'\
                            .format(sounding[ii, 0], sounding[ii, 1]**2*consts.M_0/consts.GAMMA/consts.R, sounding[ii, 1]))
                    else:
                        f.write('|       {:7.4f}    |      {:7.4f}    |  {:8.2f}  |  {:7.3f} |    {:6.4f}    |\n'\
                            .format(points[ii][0], points[ii][1], sounding[ii, 0], sounding[ii, 1]**2*consts.M_0/consts.GAMMA/consts.R, sounding[ii, 1]))

        # ray tracing function
        temp_trace, _, _ = slowscan(x_opt.xyz, np.array(xstn[j, :]), sounding, wind=setup.enable_winds, n_theta=setup.n_theta, n_phi=setup.n_phi, \
                                            precision=setup.angle_precision, tol=setup.angle_error_tol)
        time3D[j], _, _ = cyscan(x_opt.xyz, np.array(xstn[j, :]), sounding, wind=setup.enable_winds, n_theta=setup.n_theta, n_phi=setup.n_phi, \
                                           precision=setup.angle_precision, tol=setup.angle_error_tol)
        trace.append(temp_trace)
        # find residuals
        sotc[j] = tobs[j] - time3D[j]

    # for ii, element in enumerate(time3D):
    #     if np.isnan(element):
    #         w[ii] = 0

    # User defined occurrence time
    if kotc != None:
        motc = kotc
        index = []


    # elif setup.manual_fragmentation_search != '' and len(setup.manual_fragmentation_search) > 0:
    #     motc = setup.manual_fragmentation_search[3]

    # Unknown occurrence time
    else:
        index = np.isnan(sotc)
        sotc[index] = 0
        motc = np.dot(w, sotc)/sum(w)

    # Station residuals (from average)

    r = sotc - motc
    r[index] = np.nan

    return time3D, az, tf, r, motc, sotc, trace


def outputText(min_search, max_search, setup, results_arr, n_stations, s_name, xstn, tstn, w):
    
    """ Outputs the results of the search to the terminal and a text file 

    Arguments:
        min_search: [list] minimum lat, lon, and height of the search area
        max_search: [list] maximum lat, lon, and height of the search area
        single_point: [list] lat, lon, and height of a manual search point
        ref_time: [Object] datetime object of what the station arrival times are in reference to
        otc: [Object] datetime object of the time of occurrence
        kotc: [list] user-defined occurence time [HH, MM, SS]
        x_opt: [list] lat, lon, and height of the optimal supracenter
        f_opt: [float] value of error function
        n_stations: [int] number of stations
        setup: [Object] object storing user-defined setup
        s_name: [list] name of all the stations
        xstn: [list] lat, lon, and height of all the stations
        r: [list] residuals to each station
        w: [list] weights of each station
        az: [list] initial azimuth angle to each station
        tf: [list] initial takeoff angle to each station
        time3D: [list] calculated travel time to each station
        horz_dist: [list] land distance to each station
        output_name: [string] folder name to save output files in

    """
    single_point = setup.manual_fragmentation_search
    ref_time = setup.start_datetime
    output_name = setup.output_name
    run_times = len(results_arr)

    x_opt = [0]*run_times
    kotc = [0]*run_times
    otc = [0]*run_times
    f_opt = [0]*run_times
    r = [0]*run_times
    az = [0]*run_times
    tf = [0]*run_times
    time3D = [0]*run_times
    horz_dist = [0]*run_times
    
    for i in range(run_times):
        x_opt[i] = results_arr[i].x_opt
        kotc[i] = results_arr[i].kotc
        otc[i] = results_arr[i].otc
        f_opt[i] = results_arr[i].f_opt
        r[i] = results_arr[i].r
        az[i] = results_arr[i].az
        tf[i] = results_arr[i].tf
        time3D[i] = results_arr[i].time3D
        horz_dist[i] = results_arr[i].horz_dist

    x_opt = np.array(x_opt)

    ### Print Results ###
    print('......................................................................................................')
    print('The current search area extends from:')
    print('     {:7.4f} deg N, {:8.4f} deg E, {:9.2f} m Height\n'.format(*min_search))
    print('to:  {:7.4f} deg N, {:8.4f} deg E, {:9.2f} m Height\n'.format(*max_search))
    # if setup.search_type == 'pso':
    #     print('Search Type: Particle Swarm Optimization\n')
    # elif setup.search_type == 'gs':
    #     print('Search Type: Grid Search\n')
    if len(single_point) == 0:
        print('Best Fit for this search is at:\n')
    else:
        print('Given search point:\n')
    print('  {:7.4f} deg N, {:8.4f} deg E, {:9.2f} m Height\n'.format(x_opt[0, 0], x_opt[0, 1], x_opt[0, 2]))
    print('Reference Time {:02.0f}:{:02.0f}:{:02.0f}.{:06.0f}\n'.format(ref_time.hour, ref_time.minute, ref_time.second, ref_time.microsecond))
    if kotc == None:
        print('Time of Occurrence {:02.0f}:{:02.0f}:{:02.0f}.{:06.0f}\n'.format(otc[0].hour, otc[0].minute, otc[0].second, otc[0].microsecond))
    else:
        print('Known Occurrence Time {:02.0f}:{:02.0f}:{:02.0f}.{:06.0f}\n'.format(otc[0].hour, otc[0].minute, otc[0].second, otc[0].microsecond))   
    print('Value of Error function: {:8.3f}\n'.format(f_opt[0]))
    print(' ')
    # Difference between the occurence time and the reference time. Time after the reference time that the fragmentation occured
    time_diff = (otc[0].hour - ref_time.hour)*3600 + (otc[0].minute - ref_time.minute)*60 + otc[0].second-ref_time.second + (otc[0].microsecond - ref_time.microsecond)/100000
    if setup.perturb:
        print('###########################################################################################################')
        print('Perturbations')
        val, idx = min((val, idx) for (idx, val) in enumerate(f_opt))
        print('Best Perturbation: {:} Err: {:8.3f}'.format(idx, val))
        print('Time of Occurrence {:02.0f}:{:02.0f}:{:02.0f}.{:06.0f}\n'.format(otc[idx].hour, otc[idx].minute, otc[idx].second, otc[idx].microsecond))
        print('{:7.4f} (+ {:7.4f}) (- {:7.4f}) deg N'.format(x_opt[idx, 0], max(x_opt[:, 0])-x_opt[idx, 0], -min(x_opt[:, 0])+x_opt[idx, 0]))    
        print('{:7.4f} (+ {:7.4f}) (- {:7.4f}) deg E'.format(x_opt[idx, 1], max(x_opt[:, 1])-x_opt[idx, 1], -min(x_opt[:, 1])+x_opt[idx, 1]))    
        print('{:7.4f} (+ {:7.4f}) (- {:7.4f}) m'.format(x_opt[idx, 2], max(x_opt[:, 2])-x_opt[idx, 2], -min(x_opt[:, 2])+x_opt[idx, 2]))    
        print(' ')
        print('All Values')
        print('Latitude: Min: {:7.4f} deg N, Max: {:7.4f} deg N, Range: {:7.4} deg N'.format(min(x_opt[:, 0]), max(x_opt[:, 0]), max(x_opt[:, 0])-min(x_opt[:, 0])))
        print('Longitude: Min: {:7.4f} deg E, Max: {:7.4f} deg E, Range: {:7.4} deg E'.format(min(x_opt[:, 1]), max(x_opt[:, 1]), max(x_opt[:, 1])-min(x_opt[:, 1])))
        print('Elevation: Min: {:7.4f} m, Max: {:7.4f} m, Range: {:7.4} m'.format(min(x_opt[:, 2]), max(x_opt[:, 2]), max(x_opt[:, 2])-min(x_opt[:, 2])))

        print(' ')
    print('Current Residuals are as follows:')
    print('-------------------------------------------------------------------------------------------------------------------------------------------\n')
    print(' Station     Latitude (deg N)   Longitude (deg E)   Elevation(m)  Residual(s)   Weight   Azimuth(deg fN) Takeoff(deg fV)  Travel Time (s) Horz Dist (km) Arrival Time (s) \n')
    print('-------------------------------------------------------------------------------------------------------------------------------------------\n')
    for i in range(n_stations):
        r[0][i] = tstn[i] - time3D[0][i]-time_diff
        print('{:10}      {:9.4f}         {:9.4f}     {:9.2f}      {:+8.3f}      {:5.3f}     {:8.3f}      {:8.3f}        {:7.3f}        {:7.3f}       {:7.3f} \n'\
            .format(s_name[i], xstn[i, 0], xstn[i, 1], xstn[i, 2], -r[0][i], w[i], az[0][i], tf[0][i], time3D[0][i]+time_diff, horz_dist[0][i], tstn[i]))      

    # Write to file
    with open(os.path.join(output_name, '1output_data.txt'), 'w') as f:
        f.write('Results\n')
        f.write('\nTime of Execution:\n')
        f.write(datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S'))
        f.write('\n------------------------------------------------------------------------------------------------------\n')
        if setup.weather_type == 'none':
            f.write('Atmospheric Profile: Isotropic\n')
        elif setup.weather_type == 'custom':
            f.write('Atmospheric Profile: Custom\n')
        elif setup.weather_type == 'merra':
            f.write('Atmospheric Profile: MERRA-2\n')
        elif setup.weather_type == 'ecmwf':
            f.write('Atmospheric Profile: ECMWF\n')
        elif setup.weather_type == 'ukmo':
            f.write('Atmospheric Profile: UKMO\n')
        else: 
            f.write('Atmospheric Profile: ERROR\n')
        if setup.fit_type == 1:
            f.write('Fit Type: L1 Residual Fit\n')
        elif setup.fit_type == 2:
            f.write('Fit Type: L2 Residual Fit\n')
        else:
            f.write('Fit Type: Residual Fit of Order {:d}\n'.format(setup.fit_type))
            print('Warning: Only fit_type values of 1 and 2 are supported. ')
        if setup.enable_winds == False:
            f.write('Winds: Disabled\n')
        else:
            f.write('Winds: Enabled\n')
        # if setup.search_type == 'pso':
        #     f.write('Search Type: Particle Swarm Optimization\n')
        # elif setup.search_type == 'gs':
        #     f.write('Search Type: Grid Search\n')
        f.write('------------------------------------------------------------------------------------------------------\n')
        f.write('The current search area extends from:\n')
        f.write('     {:7.4f} deg N, {:8.4f} deg E, {:9.2f} m Height\n'.format(*min_search))
        f.write('to:  {:7.4f} deg N, {:8.4f} deg E, {:9.2f} m Height\n'.format(*max_search))
        if len(single_point) == 0:
            f.write('Best Fit for this search is at:\n')
        else:
            f.write('Given search point:\n')
        f.write('  {:7.4f} deg N, {:8.4f} deg E, {:9.2f} m Height\n'.format(x_opt[0, 0], x_opt[0, 1], x_opt[0, 2]))
        f.write('Reference Time:          {:02.0f}:{:02.0f}:{:02.0f}.{:06.0f}\n'.format(ref_time.hour, ref_time.minute, ref_time.second, ref_time.microsecond))
        if kotc == None:
            f.write('Time of Occurrence:      {:02.0f}:{:02.0f}:{:02.0f}.{:06.0f}\n'.format(otc[0].hour, otc[0].minute, otc[0].second, otc[0].microsecond))
        else:
            f.write('Known Occurrence Time:   {:02.0f}:{:02.0f}:{:02.0f}.{:06.0f}\n'.format(otc[0].hour, otc[0].minute, otc[0].second, otc[0].microsecond))  
        f.write('Value of Error function:     {:8.3f}\n'.format(f_opt[0]))
        f.write('\n')
        f.write('Current Residuals are as follows:\n')
        f.write('------------------------------------------------------------------------------------------------------------------------------------\n')
        f.write(' Station      Lat(deg N)    Lon(deg E)     Elevation(m)  Residual(s)   Weight   Azimuth(fN) Takeoff(fV)  Travel Time (s) H-Dist(km) Arrival Times(s)\n')
        f.write('------------------------------------------------------------------------------------------------------------------------------------\n')
        for i in range(n_stations):
            f.write('{:10}   {:9.4f}      {:9.4f}     {:9.2f}      {:+8.3f}      {:4.3f}     {:7.3f}      {:7.3f}        {:7.3f}    {:9.3f}      {:7.3f} \n'\
                .format(s_name[i], xstn[i, 0], xstn[i, 1], xstn[i, 2], -r[0][i], w[i], az[0][i], tf[0][i], time3D[0][i]+time_diff, horz_dist[0][i], tstn[i]))              
        f.close()

def scatterPlot(setup, results_arr, n_stations, xstn, s_name, dataset):

    """ Outputs a scatter plot of the search data and the optimal supracenter

    Arguments:
        single_point: [list] lat, lon, and height of a manual search point
        n_stations: [int] number of stations
        xstn: [list] lat, lon, and height of all the stations
        s_name: [list] name of all the stations
        r: [list] residuals to each station
        x_opt: [list] lat, lon, and height of the optimal supracenter
        reported_points: [list] name, lat, lon, and height of extra reference points to be plotted
        search: [list] min/max lat, lon, and height of the search boundary
        output_name: [string] folder name to save output files in
        ref_pos: [list] mean position of the stations to act as the origin for the local coordinate system
        sup: [ndarray] array of points searched with the search algorithm used for plotting
        errors: [ndarray] array of the errors in each of the points in sup for plotting

    Returns:
        min_search: [list] min lat, lon and height of the search area
        max_search: [list] max lat, lon and height of the search area
    """

    single_point = setup.manual_fragmentation_search
    r = results_arr[0].r
    x_opt = results_arr[0].x_opt
    reported_points = setup.reported_points
    search = setup.search
    output_name = setup.output_name
    ref_pos = [setup.ref_pos.lat, setup.ref_pos.lon, setup.ref_pos.elev]
    sup = results_arr[0].sup
    errors = results_arr[0].errors
    run_times = len(results_arr)

    ### Convert back to lat/lon
    # Search Region boundaries
    min_search = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], [search[0], search[2], search[4]])
    max_search = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], [search[1], search[3], search[5]])
    # Set up axes
    fig = plt.figure(figsize=plt.figaspect(0.5))
    fig.set_size_inches(20.9, 11.7)

    # automatic search
    if len(single_point) == 0:
        ax1 = fig.add_subplot(1, 2, 1, projection='3d')
        ax2 = fig.add_subplot(1, 2, 2, projection='3d')

    # manual search
    else:
        ax1 = fig.add_subplot(1, 1, 1, projection='3d')

    ### Labels
    ax1.set_title("Supracenter Locations")
    ax1.set_xlabel("Latitude (deg N)", linespacing=3.1)
    ax1.set_ylabel("Longitude (deg E)", linespacing=3.1)
    ax1.set_zlabel('Elevation (m)', linespacing=3.1)

    # plot station names and residuals
    for h in range(n_stations):

        # Convert station locations to geographic
        xstn[h, 0], xstn[h, 1], xstn[h, 2] = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], xstn[h, :])

        # Add station names
        ax1.text(xstn[h, 0], xstn[h, 1], xstn[h, 2],  '%s' % (s_name[h]), size=10, zorder=1, color='k')
        if len(single_point) == 0:
            ax2.text(xstn[h, 0], xstn[h, 1], xstn[h, 2],  '%s' % (s_name[h]), size=10, zorder=1, color='k')

    # Add stations with color based off of residual
    res = ax1.scatter(xstn[:, 0], xstn[:, 1], xstn[:, 2], c=abs(r), marker='^', cmap='viridis_r', depthshade=False)

    # Add point and label
    for j in range(run_times):
        if j == 0:
            ax1.scatter(x_opt[0], x_opt[1], x_opt[2], c = 'r', marker='*')
            ax1.text(x_opt[0], x_opt[1], x_opt[2], '%s' % ('Supracenter'), zorder=1, color='k')
        else:
            ax1.scatter(results_arr[j].x_opt[0], results_arr[j].x_opt[1], results_arr[j].x_opt[2], alpha=0.2, c = 'r', marker='*')

    # Get the limits of the plot

    x_min, y_min = search[0], search[2]
    x_max, y_max = search[1], search[3]

    img_dim = 30

    x_data = np.linspace(x_min, x_max, img_dim)
    y_data = np.linspace(y_min, y_max, img_dim)
    xx, yy = np.meshgrid(x_data, y_data)

    # Make an array of all plane coordinates
    plane_coordinates = np.c_[xx.ravel(), yy.ravel(), np.zeros_like(xx.ravel())]

    times_of_arrival = np.zeros_like(xx.ravel())

    x_opt = geo2Loc(ref_pos[0], ref_pos[1], ref_pos[2], x_opt[0], x_opt[1], x_opt[2])

    print('Creating contour plot...')
    # Calculate times of arrival for each point on the reference plane
    for i, plane_coords in enumerate(plane_coordinates):

        #plane_coords = geo2Loc(ref_pos[0], ref_pos[1], ref_pos[2], plane_coords[0], plane_coords[1], plane_coords[2])
        
        # Create interpolated atmospheric profile for use with cyscan
        sounding, points = getWeather(x_opt, plane_coords, setup.weather_type, ref_pos, copy.copy(dataset))

        # Use distance and atmospheric data to find path time
        ti, _, _ = cyscan(np.array(x_opt), np.array(plane_coords), sounding, \
                                wind=setup.enable_winds, n_theta=setup.n_theta, n_phi=setup.n_phi, \
                                precision=setup.angle_precision)
        if np.isnan(ti):
            #Make blank contour
            ti = -1
        times_of_arrival[i] = ti

    times_of_arrival = times_of_arrival.reshape(img_dim, img_dim)

    # Determine range and number of contour levels, so they are always centred around 0
    toa_abs_max = np.max([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
    toa_abs_min = np.min([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
    levels = np.linspace(toa_abs_min, toa_abs_max, 50)
    print(levels)

    # Plot colorcoded times of arrival on the surface
    toa_conture = ax1.contourf(xx, yy, times_of_arrival, levels, cmap='inferno', alpha=1.0)
    # Add a color bar which maps values to colors
    fig.colorbar(toa_conture, label='Time of arrival (s)', ax=ax1)


    # plot given reference points
    for q in range(len(reported_points[:])):
        ax1.scatter(reported_points[q][1], reported_points[q][2], reported_points[q][3], c='k', marker='X')
        ax1.text(reported_points[q][1], reported_points[q][2], reported_points[q][3], '%s' % (reported_points[q][0]), size=7, zorder=1, color='k')

    if len(single_point) == 0:
        # potential Supracenter points with color based on error
        for i in range(len(sup)):
            sup[i, 0], sup[i, 1], sup[i, 2] = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], sup[i, :])
        sc = ax2.scatter(sup[:, 0], sup[:, 1], sup[:, 2], c=errors, cmap='inferno_r', depthshade=False)
        a = plt.colorbar(sc, ax=ax2)
        a.set_label("Error in Supracenter (s)")
        ax2.scatter(xstn[:, 0], xstn[:, 1], xstn[:, 2], c=abs(r), marker='^', cmap='viridis_r', depthshade=False)
        ax2.set_title("Fits of Potential Supracenters")
        ax2.set_xlabel("Latitude (deg N)")
        ax2.set_ylabel("Longitude (deg E)")
        ax2.set_zlabel('Elevation (m)')

    if len(single_point) == 0 and setup.restrict_to_trajectory:
        # A = [0, 0, 0]
        # B = [0, 0, 0]
        # A[0], A[1], A[2] = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], [setup.lat_i, setup.lon_i, setup.elev_i])
        # B[0], B[1], B[2] = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], [setup.lat_f, setup.lon_f, setup.elev_f])
        ax1.plot([setup.lat_i, setup.lat_f], [setup.lon_i, setup.lon_f], [setup.elev_i*1000, setup.elev_f*1000], c='g')  
        ax2.plot([setup.lat_i, setup.lat_f], [setup.lon_i, setup.lon_f], [setup.elev_i*1000, setup.elev_f*1000], c='g')

    elif setup.restrict_to_trajectory:
        ax1.plot([setup.lat_i, setup.lat_f], [setup.lon_i, setup.lon_f], [setup.elev_i*1000, setup.elev_f*1000], c='g')  

    # colorbars
    b = plt.colorbar(res, ax=ax1)
    b.set_label("Station Residuals (s)")

    # automatic search
    if len(single_point) == 0:
        ax1.set_xlim3d(min_search[0], max_search[0])
        ax2.set_xlim3d(min_search[0], max_search[0])

        ax1.set_ylim3d(min_search[1], max_search[1])
        ax2.set_ylim3d(min_search[1], max_search[1])

        ax1.set_zlim3d(0, max_search[2])
        ax2.set_zlim3d(0, max_search[2])

    # manual search, don't pull up pso graph
    else:
        ax1.set_xlim3d(min_search[0], max_search[0])
        ax1.set_ylim3d(min_search[1], max_search[1])
        ax1.set_zlim3d(0, max_search[2])

    plt.tight_layout()
    plt.savefig(os.path.join(output_name, '1output_graph.png'))
    plt.show()
    
    plt.clf()
    plt.close()

    return min_search, max_search

def residPlot(results_arr, s_name, xstn, output_name, n_stations):
    """ outputs a 2D residual plot of the stations with the optimal supracenter

    Arguments:
        x_opt: [list] optimal supracenter position
        s_name: [list] list of station names
        xstn: [list] list of lat, lon, height positions of each station
        resid: [list] list of residuals to each station
        output_name: [string] folder to store the data in
        n_stations: [int] number of stations
    """
    x_opt = results_arr[0].x_opt
    resid = results_arr[0].r
    run_times = len(results_arr)
 # Residuals Plot
    fig = plt.figure()
    ax3 = fig.add_subplot(111)

    res = ax3.scatter(xstn[:, 0], xstn[:, 1], c=abs(resid), marker='^', cmap='viridis_r', s=21)
    
    for h in range(n_stations):

        # Add station names
        ax3.text(xstn[h, 0], xstn[h, 1],  '%s' % (s_name[h]), size=10, zorder=1, color='k')

    # Add point and label
    for i in range(run_times):
        if i == 0:
            ax3.scatter(x_opt[0], x_opt[1], c = 'r', marker='*', s=21)
            ax3.text(x_opt[0], x_opt[1],  '%s' % ('Supracenter'), size=10, zorder=1, color='k')
        else:
            ax3.scatter(results_arr[i].x_opt[0], results_arr[i].x_opt[1], alpha=0.2, c = 'r', marker='*', s=21)


    ax3.set_xlabel("Latitude (deg N)")
    ax3.set_ylabel("Longitude (deg E)")
    ax3.set_title("Station Residuals")

    c = plt.colorbar(res, ax=ax3)
    c.set_label("Station Residuals (s)")

    plt.savefig(os.path.join(output_name, '1output_residuals.png'))
    plt.show()