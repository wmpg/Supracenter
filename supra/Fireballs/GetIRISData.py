""" Functions for fetching USarray waveforms. """


from __future__ import print_function, division, absolute_import

import os
import sys
import datetime
import argparse
import copy
import urllib
import time

# Check version
if sys.version_info.major < 3:
    import urllib as urllibrary

else:
    import urllib.request as urllibrary

import numpy as np
import matplotlib.pyplot as plt
import scipy.signal
import obspy

from wmpl.Utils.Earth import greatCircleDistance
from wmpl.Utils.OSTools import mkdirP
from wmpl.Utils.PlotMap import GroundMap
from wmpl.Utils.Math import subsampleAverage

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Fireballs.SeismicTrajectory import timeOfArrival, waveReleasePoint, parseWeather, Constants
from supra.Fireballs.Program import configRead, configParse, position, station
from supra.Supracenter.cyscan import cyscan
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.SPPT import perturb
from supra.Supracenter.angleConv import geo2Loc, loc2Geo

DATA_FILE = 'data.txt'
C = ['r', 'g', 'm', 'k', 'y']
'''Reads input config files'''

try:
    # Python 2  
    import ConfigParser as configparser

except:
    # Python 3
    import configparser


import os
import sys

import datetime 
    

def getIRISStations(lat_centre, lon_centre, deg_radius, start_date, end_date, network='*', channel='BDF'):
    """ Retrieves seismic stations from the IRIS web service in a radius around given geographical position 
        which were active during the specified range of dates.

    Arguments:
        lat_centre: [float] Latitude of the query centre (+N, degrees).
        lon_centre: [float] Longitude of the query centre (+E, degrees).
        deg_radius: [float] Query radius (degrees).
        start_date: [str] First date when a station was recording (YYYY-MM-DD format).
        end_date: [str] Final date when a station was recording (YYYY-MM-DD format).
        
    Keyword arguments:
        network: [str] Filter retrieved stations by the given seismic network code. * by default, meaning
            stations for all networks will be retrieved.
        channel: [str] Seismograph channel. BDF by default.
        orfeus: [bool] If True, the European ORFEUS website will be queried instead of IRIS.

    Return:
        [list] A list of stations and their parameters.
    """

    # Use the European ORFEUS data access site
    try:
        # Construct ORFEUS URL
        iris_url = ("http://www.orfeus-eu.org/fdsnws/station/1/query?network={:s}&latitude={:.3f}&longitude={:.3f}" \
            "&maxradius={:.3f}&start={:s}&end={:s}&channel={:s}&format=text" \
            "&includerestricted=false&nodata=404").format(network, lat_centre, lon_centre, deg_radius, \
            start_date, end_date, channel)
    except:
        pass
        #print('Unable to get ORFEUS data!')

    try:
        # Construct IRIS URL
        iris_url2 = ("http://service.iris.edu/fdsnws/station/1/query?net={:s}&latitude={:.3f}&longitude={:.3f}" \
            "&maxradius={:.3f}&start={:s}&end={:s}&cha={:s}&nodata=404&format=text" \
            "&matchtimeseries=true").format(network, lat_centre, lon_centre, deg_radius, start_date, \
            end_date, channel)
    except:
        #print('Unable to get IRIS data!')
        pass

    # Initialize station arrays for both IRIS and ORFEUS
    stations_txt_1 = []
    stations_txt_2 = []

    station_list = []

    # Retrieve station list if there is data available
    # ORFEUS
    try:
        stations_txt_1 = urllibrary.urlopen(iris_url).read().decode('utf-8')
    except:
        #print('Unable to get ORFEUS data!')
        pass

    # IRIS
    try:
        stations_txt_2 = urllibrary.urlopen(iris_url2).read().decode('utf-8')
    except:
        #print('Unable to get IRIS data!')
        pass


    # Return an empty list if no stations were retrieved
    if not stations_txt_1 and not stations_txt_2:
        return station_list

    # Parse the ORFEUS stations
    if len(stations_txt_1) != 0:
        for entry in stations_txt_1.split('\n')[1:]:

            entry = entry.split('|')

            # Skip empty rows
            if len(entry) != 8:
                continue

            # Unpack the line
            network, station_code, lat, lon, elev, station_name, start_work, end_work = entry

            station_list.append([network, station_code, float(lat), float(lon), float(elev), station_name])

    # Parse the IRIS stations
    if len(stations_txt_2) != 0:
        for entry in stations_txt_2.split('\n')[1:]:

            entry = entry.split('|')

            # Skip empty rows
            if len(entry) != 8:
                continue
            
            # Unpack the line
            network, station_code, lat, lon, elev, station_name, start_work, end_work = entry

            station_list.append([network, station_code, float(lat), float(lon), float(elev), station_name])

    return station_list


def getIRISWaveformFiles(network, station_code, start_datetime, end_datetime, dir_path='.', channel='BDF'):
    """ Download weaveform files from the IRIS site. 
    
    Arguments:
        network: [str] Network code.
        station_code: [str] Station code.
        start_datetime: [datetime object] Datetime of the beginning of the data chunk to retrieve.
        end_datetime: [datetime object] Datetime of the end of the data chunk to retrieve. This cannot be
            more than 2 hours of data (this is not an IRIS limitation, but a limitation set here, so nobody
            clogs the IRIS servers with unreasonable requests).
        
    Keyword arguments:
        dir_path: [str] Full path to location where the miniSEED files will be saved.
        channel: [str] Seismograph channel. BDF by default.
        orfeus: [bool] If True, the European ORFEUS website will be queried instead of IRIS.

    Return:
        mseed_file_path: [str] Path to the downloaded miniSEED data file.
    """

    # Make strings from datetime objects
    sd = start_datetime
    start_time = "{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}".format(sd.year, sd.month, sd.day, \
        sd.hour, sd.minute, sd.second, sd.microsecond//1000)

    ed = end_datetime
    end_time = "{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}".format(ed.year, ed.month, ed.day, \
        ed.hour, ed.minute, ed.second, ed.microsecond//1000)

    #start_time = "2010-02-27T06:30:00.000"
    #end_time = "2010-02-27T10:30:00.000"

    # Check that no more then 2 hours of data is requested
    if (end_datetime - start_datetime).total_seconds() > 2*3600:
        print('No more than 2 hours of waveform data can be requested!')
        return None


    # Use the European ORFEUS data access site if specified
    try:
        # Construct ORFEUS URL
        iris_url = ("http://www.orfeus-eu.org/fdsnws/dataselect/1/query?network={:s}&station={:s}" \
                "&channel={:s}&start={:s}&end={:s}").format(network, station_code, channel, start_time, \
                end_time)
    except:
        print('Unable to access ORFEUS!')

    try:
        # Construct IRIS URL
        iris_url2 = ("http://service.iris.edu/fdsnws/dataselect/1/query?net={:s}&sta={:s}&cha={:s}" \
                "&start={:s}&end={:s}").format(network, station_code, channel, start_time, end_time)
    except:
        print('Unable to access IRIS!')


    # Construct a file name
    mseed_file = network + '_' + station_code + '_{:s}_'.format(channel) + start_time.replace(':', '.') \
        + '_' + end_time.replace(':', '.') + '.mseed'

    mseed_file_path = os.path.join(dir_path, mseed_file)

    # Get the miniSEED file

    # ORFEUS
    try:
        iris_file = urllibrary.urlopen(iris_url)
    except:
        print('Unable to access ORFEUS!')
        iris_file = False

    # IRIS
    try:
        iris_file2 = urllibrary.urlopen(iris_url2)
    except:
        print('Unable to access IRIS!')
        iris_file2 = False
    

    # Make sure both contents are not empty
    if iris_file or iris_file2:

        # Save the file if it has any length
        with open(mseed_file_path,'wb') as f:
            
            try:
                if iris_file:
                    f.write(iris_file.read())

                if iris_file2:
                    f.write(iris_file2.read())
                    
            except urllibrary.URLError:
                print('Connection error! Could not download the waveform!')

        if os.stat(mseed_file_path).st_size == 0:
                
            print("\rWarning {:}-{:} {:} is an empty file!                      ".format(network, station_code, channel))
            
            os.remove(mseed_file_path)
            return None
        else:
            return mseed_file

    else:
        return None



def writeStationAndWaveformsListFile(data_list, file_path):
    """ Writes the list of stations and miniSEED file names to disk. """

    with open(file_path, 'w') as f:
        
        for entry in data_list:
        
            #network, station_code, stat_lat, stat_lon, stat_elev, station_name, mseed_file = entry

            f.write('{:s}|{:s}|{:.6f}|{:.6f}|{:.6f}|{:s}|{:s}|{:s}\n'.format(*entry))



def readStationAndWaveformsListFile(file_path, rm_stat=[]):
    """ Reads the list of stations and miniSEED file names to disk. """

    if os.path.isfile(file_path):

        with open(file_path) as f:

            stn_list = []

            for line in f:

                if not line:
                    continue

                line = line.replace('\n', '')
                line = line.split('|')

                try:
                    network, station_code, stat_lat, stat_lon, stat_elev, station_name, channel, mseed_file = line
                except:
                    print("ERROR: Waveform files downloaded previous to Jan. 8, 2019 are lacking a channel tag added. Redownloading the waveform files will likely fix this")

                stat_pos = position(float(stat_lat), float(stat_lon), float(stat_elev))
                stn = station(network, station_code, stat_pos, channel, station_name, mseed_file)

                if station_code not in rm_stat:
                    stn_list.append(stn)
                else:
                    print('Excluding station: {:}'.format(network + '-' + station_code))

            return stn_list

    else:
        return []


def getAllWaveformFiles(lat_centre, lon_centre, deg_radius, start_datetime, end_datetime, network='*', \
    channel='BDF', dir_path='.'):
    """ Retrieves and saves waveforms as miniSEED files of all seismic stations from the IRIS web service in 
        a radius around given geographical position and for the given range of times.

    Arguments:
        lat_centre: [float] Latitude of the query centre (+N, degrees).
        lon_centre: [float] Longitude of the query centre (+E, degrees).
        deg_radius: [float] Query radius (degrees).
        start_datetime: [datetime object] Datetime of the beginning of the data chunk to retrieve.
        end_datetime: [datetime object] Datetime of the end of the data chunk to retrieve. This cannot be
            more than 2 hours of data (this is not an IRIS limitation, but a limitation set here, so nobody
            clogs the IRIS servers with unreasonable requests).
        
    Keyword arguments:
        network: [str] Filter retrieved stations by the given seismic network code. * by default, meaning
            stations for all networks will be retrieved.
        dir_path: [str] Full path to location where the miniSEED files will be saved.
        orfeus: [bool] If True, the European ORFEUS website will be queried instead of IRIS.

    """


    # Station activity date range
    sd = start_datetime
    start_date = "{:04d}-{:02d}-{:02d}".format(sd.year, sd.month, sd.day)
    ed = sd + datetime.timedelta(days=1)
    end_date = "{:04d}-{:02d}-{:02d}".format(ed.year, ed.month, ed.day)


    # Make the data directory
    mkdirP(dir_path)


    # Get a list of stations active on specified dates around the given location
    station_listBDF = getIRISStations(lat_centre, lon_centre, deg_radius, start_date, end_date, \
        network=network, channel='BDF')

    station_listBHZ = getIRISStations(lat_centre, lon_centre, deg_radius, start_date, end_date, \
        network=network, channel='BHZ')

    station_listHHZ = getIRISStations(lat_centre, lon_centre, deg_radius, start_date, end_date, \
        network=network, channel='HHZ')

    print('DOWNLOADED STATIONS:')
    print('BDF')
    if len(station_listBDF) == 0:
        print("No Stations in BDF")
    for i in range(len(station_listBDF)):
        print('{:2}-{:6} Lat: {:7.6} Lon: {:7.6} Elev: {:7.6} Loc: {:}'.format(*station_listBDF[i]))

    print('BHZ')
    if len(station_listBHZ) == 0:
        print("No Stations in BHZ")
    for i in range(len(station_listBHZ)):
        print('{:2}-{:6} Lat: {:7.6} Lon: {:7.6} Elev: {:7.6} Loc: {:}'.format(*station_listBHZ[i]))

    print('HHZ')
    if len(station_listHHZ) == 0:
        print("No Stations in HHZ")
    for i in range(len(station_listHHZ)):
        print('{:2}-{:6} Lat: {:7.6} Lon: {:7.6} Elev: {:7.6} Loc: {:}'.format(*station_listHHZ[i]))

    #Combine into one list of all stations
    station_list = station_listBDF + station_listBHZ + station_listHHZ

    # A list of station data and waveform files
    data_list = []

    no_stats = len(station_list)
    # Go through all stations, retrieve and save waveforms
    for ii, station_data in enumerate(station_list):

        # Unpack station info
        network, station_code, stat_lat, stat_lon, stat_elev, station_name = station_data

        sys.stdout.write("\rDownloading Station Data: {:5.2f} % {:}-{:}".format((ii+1)/(no_stats)*100, network, station_code))
        sys.stdout.flush()
        time.sleep(0.001)
        
        if station_data in station_listBDF:
            # Retreive the waveform of the given station
            mseed_file = getIRISWaveformFiles(network, station_code, start_datetime, end_datetime, \
                channel='BDF', dir_path=dir_path)
            station_data.append('BDF')
        elif station_data in station_listHHZ:
            mseed_file = getIRISWaveformFiles(network, station_code, start_datetime, end_datetime, \
                channel='HHZ', dir_path=dir_path)
            station_data.append('HHZ')
        elif station_data in station_listBHZ:
            mseed_file = getIRISWaveformFiles(network, station_code, start_datetime, end_datetime, \
                channel='BHZ', dir_path=dir_path)
            station_data.append('BHZ')
        else:
            print("WARNING: Cannot find station {:}".format(station_code))

        if mseed_file != None:
            station_data.append(mseed_file)
            data_list.append(station_data)

    # Save the list of station parameters and data files to disk
    writeStationAndWaveformsListFile(data_list, os.path.join(dir_path, DATA_FILE))
    print('')
    print('Data file: ', DATA_FILE, 'written!')



def butterworthBandpassFilter(lowcut, highcut, fs, order=5):
    """ Butterworth bandpass filter.

    Argument:
        lowcut: [float] Lower bandpass frequency (Hz).
        highcut: [float] Upper bandpass frequency (Hz).
        fs: [float] Sampling rate (Hz).

    Keyword arguments:
        order: [int] Butterworth filter order.

    Return:
        (b, a): [tuple] Butterworth filter.

    """

    # Calculate the Nyquist frequency
    nyq = 0.5*fs

    low = lowcut/nyq
    high = highcut/nyq

    # Init the filter
    b, a = scipy.signal.butter(order, [low, high], btype='band')

    return b, a



def convolutionDifferenceFilter(waveform_data):
    """ Apply the convolution filter on data as suggested in Kalenda et al. (2014). """

    # Apply the filter 
    filtered_data = np.convolve(waveform_data, [-0.5, 1.0, -0.5], mode='same')

    # Detrend data
    filtered_data = filtered_data - np.mean(filtered_data)

    return filtered_data


def plotStationMap(dir_path, data_list, lat_centre, lon_centre, setup, sounding, ax=None, isc_data=None):
    """ Plots the map of siesmic stations from loaded data file. """

    fig = plt.figure(figsize=plt.figaspect(0.5))
    fig.set_size_inches(20.9, 11.7)
    if ax is None:
        ax = plt.gca()

    # Find unique networks
    # networks = [entry[0] for entry in data_list]
    # stat = [entry[1] for entry in data_list]

    # net_isc = []

    # lats=[]
    # lons=[]
    # Extra stations 
    if isc_data is not None:

        all_stns = data_list + isc_data

        # Remove duplicates
        # k = sorted(isc_data)
        # isc_data = [k[i] for i in range(len(k)) if i == 0 or k[i] != k[i-1]]

        # for line in isc_data:
            
        #     # Only use stations within 5 degrees of lat and lon
        #     if abs(line[2] - lat_centre) < 5 and abs(line[3] - lon_centre) < 5:

        #         lats.append(np.radians(line[2]))
        #         lons.append(np.radians(line[3]))
        #         net_isc.append(line[5])

    # # Extract the list of station locations
    # lat_list = [np.radians(entry[2]) for entry in data_list]
    # lon_list = [np.radians(entry[3]) for entry in data_list]

    if len(all_stns) == 0:
        print("ERROR: No stations to plot!")
        exit()

    lats = []
    lons = []
    for i in range(len(all_stns)):
        lats.append(all_stns[i].position.lat_r)
        lons.append(all_stns[i].position.lon_r)

    # Plot stations and extra stations
    m = GroundMap(lats, lons, ax=ax, color_scheme='light')

    # Plot different networks with different colours
    for stn in all_stns:

        # # Extract the list of station locations
        # lat_net_list = [np.radians(entry[2]) for entry in data_list]
        # lon_net_list = [np.radians(entry[3]) for entry in data_list]

        m.scatter(stn.position.lat_r, stn.position.lon_r, s=2, label=stn.network)

        # for i in range(len(lat_net_list)):
        x, y = m.m(stn.position.lon, stn.position.lat)

        plt.text(x, y, stn.network + '-' + stn.code, horizontalalignment='left', verticalalignment='top', color='k', fontsize=8)
                
        # if stat[i] in setup.rm_stat:
        #     pass
        #     # print('Excluding station: {:}'.format(networks[i] + '-' + stat[i]))
        # else:
        #     if stat[i] in setup.high_f:
        #         m.scatter(lat_net_list[i], lon_net_list[i], s=25, c='g')
        #     elif stat[i] in setup.high_b:
        #         m.scatter(lat_net_list[i], lon_net_list[i], s=25, c='b')

    # # if len(lats) != 0:
    # for i in range(len(net_isc)):

    #     x, y = m.m(np.degrees(lons[i]), np.degrees(lats[i]))

    #     plt.text(x, y, net_isc[i], horizontalalignment='left', verticalalignment='top', color='k', fontsize=8)

    lx, ly = m.m(lon_centre, lat_centre)

    # # All extra stations added
    # if isc_data is not None:
        
    #     for i in range(len(net_isc)):

    #         # Convert coordinates to map coordinates
    #         x, y = m.m(np.degrees(lons[i]), np.degrees(lats[i]))

    #         # Plot extra stations
    #         m.scatter(lats[i], lons[i], marker='^', c='k', s=1, )

    #         # Plot the text 
    #         #plt.text(x, y, net_isc[i], horizontalalignment='left', verticalalignment='top', color='k', fontsize=8)
        
    #         data_list.append(isc_data[i])

    # Plot source location    
    m.scatter([np.radians(lat_centre)], [np.radians(lon_centre)], marker='*', c='yellow', edgecolor='k', \
        linewidth=0.1, label='Source')


    # Plot the trajectory or fragmentation point if given
    if setup.show_fragmentation_waveform or setup.show_ballistic_waveform:

        if setup.show_fragmentation_waveform:

            for i, line in enumerate(setup.fragmentation_point):

                # Fragmentation plot
                m.scatter([np.radians(float(line[0]))], [np.radians(float(line[1]))], c=C[(i+1)%4], marker='x')

        # Extract coordinates of the reference station
        ref_pos = position(lat_centre, lon_centre, 0)

        # # Calculate the coordinates of the trajectory intersection with the ground
        # lat_i, lon_i, elev_i = local2LatLon(float(np.radians(lat0)), float(np.radians(lon0)), float(0), \
        #     np.array([float(setup.lat_f), float(setup.lon_f), 0]))

        # Calculate the coordinate of the beginning of the trajectory
        # lat_beg, lon_beg = np.radians(float(np.degrees(setup.lat_i)) - np.cos(np.radians(setup.azim))), \
        #                    np.radians(float(np.degrees(setup.lon_i)) - np.sin(np.radians(setup.azim)))

        if setup.show_ballistic_waveform:

            # Plot intersection with the ground
            m.scatter(setup.traj_f.lat_r, setup.traj_f.lon_r, s=10, marker='x', c='b')

            # Plot the trajectory
            m.plot([setup.traj_i.lat_r, setup.traj_f.lat_r], [setup.traj_i.lon_r, setup.traj_f.lon_r], c='b')

            # Get the limits of the plot
            # (approximately a box around the deg_radius)
            x_min = setup.traj_f.lon - 100000*setup.deg_radius
            x_max = setup.traj_f.lon + 100000*setup.deg_radius
            y_min = setup.traj_f.lat - 100000*setup.deg_radius
            y_max = setup.traj_f.lat + 100000*setup.deg_radius

            # Grid size of the contour plot
            img_dim = setup.contour_res
            x_data = np.linspace(x_min, x_max, img_dim)
            y_data = np.linspace(y_min, y_max, img_dim)
            xx, yy = np.meshgrid(x_data, y_data)

            # # Make an array of all plane coordinates
            plane_coordinates = np.c_[xx.ravel(), yy.ravel(), np.zeros_like(xx.ravel())]

            times_of_arrival = np.zeros_like(xx.ravel())

            # print('Creating contour plot...')
            # # Calculate times of arrival for each point on the reference plane
            # Calculate times of arrival for each point on the reference plane
            az = np.radians(setup.azim)
            ze = np.radians(setup.zangle)

            # vector of the fireball trajectory
            traj_vect = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

            #traj_vect = np.array([np.cos(az)*np.cos(ze), np.sin(az)*np.cos(ze), -np.sin(ze)])

            for i, plane_coords in enumerate(plane_coordinates):
                
                # Print out percent done
                if (i + 1) % 10 == 0:
                    sys.stdout.write("\rDrawing Contour: {:.2f} %".format(100*(i + 1)/img_dim**2))
                    sys.stdout.flush()
                    time.sleep(0.001)

                setup.traj_f.pos_loc(ref_pos)

                # Point on the trajectory where the plane coordinate arrival came from
                p = waveReleasePoint(plane_coords, setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, np.radians(setup.azim), \
                                            np.radians(setup.zangle), setup.v_sound)

                # Coordinate transformation (rotate 90 deg CCW)
                #p[0], p[1] = -p[1], p[0]

                # vector between the wave release point and the plane coordinate
                d_vect = plane_coords - p

                # Since the arrivals are always perpendicular to the fireball trajectory, only take arrivals where the dot product
                # of the vectors are small. This may not hold true for weather?
                #print(np.dot(d_vect, traj_vect))
                #if np.dot(d_vect, traj_vect) < setup.dot_tol:

                ti = timeOfArrival(plane_coords, setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, np.radians(setup.azim), \
                                                np.radians(setup.zangle), setup, sounding=sounding, ref_loc=[ref_pos.lat_r, ref_pos.lon_r, ref_pos.elev], travel=True, fast=True)

                # escape value for when there is no arrival
                #else:
                #    ti = np.nan


                times_of_arrival[i] = ti + setup.t0

            # if there is no arrival, set to the maximum value on the contour
            max_time = np.nanmax(times_of_arrival)
            for i in range(len(times_of_arrival)):
                if np.isnan(times_of_arrival[i]):
                    times_of_arrival[i] = max_time

            times_of_arrival = times_of_arrival.reshape(img_dim, img_dim)

            # Determine range and number of contour levels, so they are always centred around 0
            toa_abs_max = np.max([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
            #  toa_abs_min = np.min([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
            levels = np.linspace(0, toa_abs_max, 25)

            ### Convert contour local coordinated to geo coordinates
            lat_cont = []
            lon_cont = []

            for x_cont, y_cont in zip(xx.ravel(), yy.ravel()):
                
                lat_c, lon_c, _ = loc2Geo(ref_pos.lat, ref_pos.lon, ref_pos.elev, np.array([x_cont, y_cont, 0]))

                lat_cont.append(lat_c)
                lon_cont.append(lon_c)

            lat_cont = np.array(lat_cont).reshape(img_dim, img_dim)
            lon_cont = np.array(lon_cont).reshape(img_dim, img_dim)

            try:
                # Plot the time of arrival contours
                toa_conture = m.m.contourf(lon_cont, lat_cont, times_of_arrival, levels, zorder=3, \
                    latlon=True, cmap='viridis_r', alpha=0.3)
                m.m.colorbar(toa_conture, label='Time of arrival (s)')
            except:
                print("WARNING: Unable to plot the contour for the ballistic trajectory.")

            # # Plot colorcoded times of arrival on the surface
            # toa_conture = self.m.m.contour(xx, yy, times_of_arrival, levels, cmap='inferno',\
            #      alpha=1.0, zorder=3, latlon=False)
            # # Add a color bar which maps values to colors
            


    ax.set_title('Source location: {:.6f}, {:.6f}'.format(lat_centre, lon_centre))
    
    #plt.savefig('/home/luke/Desktop/_map.png', dpi=300)
    #plt.savefig(os.path.join(setup.output_folder, 'map.png'), dpi=300)


def plotAllWaveforms(dir_path, stn_list, setup, sounding, ax=None, waveform_window=None,\
    difference_filter_all=False):
    """ Bandpass filter and plot all waveforms from the given data list. 

    Keyword arguments:
        waveform_window: [int] If given, the waveforms will be cut around the modelled time of arrival line
            with +/- waveform_window/2 seconds. None by default, which means the whole waveform will be 
            plotted.
        difference_filter_all: [bool] If True, the Kalenda et al. (2014) difference filter will be applied
                on the data plotted in the overview plot of all waveforms.
    """

    # Initialize variables
    v_sound = setup.v_sound
    t0 = setup.t0
    lat_centre = setup.lat_centre
    lon_centre = setup.lon_centre

    if ax is None:
        ax = plt.gca()

    max_wave_value = 0
    min_wave_value = np.inf
    min_time = np.inf
    max_time = 0

    # # Add extra stations from the config file
    # if setup.stations is not None:
    #     for line in setup.stations:

    #         # Prevent adding duplicates
    #         if line in data_list:
    #             continue

    #         data_list.append(line)

    lats = []
    lons = []
    for i in range(len(stn_list)):
        lats.append(stn_list[i].position.lat)
        lons.append(stn_list[i].position.lon)

    # Azimuth from source point to station (degrees +N of due E)
    az = np.arctan2(lat_centre - np.array(lats), lon_centre - np.array(lons))
    # az - normalized values for color-coding, 
    # az_n - original values for text

    az_n = copy.copy(az)
    #az_n = (90 - np.degrees(az_n)%360)%360

    # normalize azimuths
    for i in range(len(az)):
        az[i] += abs(min(az))
        az[i] /= (max(az) + abs(min(az)))

    # Go though all stations and waveforms
    bad_stats = []

    for idx, stn in enumerate(stn_list):


        sys.stdout.write('\rPlotting: {:} {:}              '.format(stn.network, stn.code))
        sys.stdout.flush()
        time.sleep(0.001)

        mseed_file_path = os.path.join(dir_path, stn.file_name)

        try:
            
            # Read the miniSEED file
            if os.path.isfile(mseed_file_path):

                mseed = obspy.read(mseed_file_path)

            else:
                bad_stats.append(idx)
                print('File {:s} does not exist!'.format(mseed_file_path))
                continue


        except TypeError as e:
            bad_stats.append(idx)
            print('Opening file {:} failed with error: {:}'.format(mseed_file_path, e))
            continue

        # Find channel with BHZ, HHZ, or BDF

        for i in range(len(mseed)):
            if mseed[i].stats.channel == 'BDF':
                stn.channel = 'BDF'
                stream = i

        for i in range(len(mseed)):
            if mseed[i].stats.channel == 'BHZ':
                stn.channel = 'BHZ'
                stream = i

        for i in range(len(mseed)):
            if mseed[i].stats.channel == 'HHZ':
                stn.channel = 'HHZ'
                stream = i

        for i in range(len(mseed)):
            if mseed[i].stats.channel == 'EHZ':
                stn.channel = 'EHZ'
                stream = i


        # Unpack miniSEED data
        delta = mseed[stream].stats.delta
        waveform_data = mseed[stream].data

        # Extract time
        start_datetime = mseed[stream].stats.starttime.datetime
        end_datetime = mseed[stream].stats.endtime.datetime

        stn.offset = (start_datetime - setup.start_datetime).total_seconds()

        # Skip stations with no data
        if len(waveform_data) == 0:
            continue

        # Apply the Kalenda et al. (2014) difference filter instead of Butterworth
        if difference_filter_all:

            waveform_data = convolutionDifferenceFilter(waveform_data)

        else:

            ### BANDPASS FILTERING ###

            # Init the butterworth bandpass filter
            butter_b, butter_a = butterworthBandpassFilter(0.8, 5.0, 1.0/delta, order=6)

            # Filter the data
            waveform_data = scipy.signal.filtfilt(butter_b, butter_a, waveform_data)

            # Average and subsample the array for quicker plotting (reduces 40Hz to 10Hz)
            waveform_data = subsampleAverage(waveform_data, 4)
            delta *= 4

            ##########################


        # Calculate the distance from the source point to this station (kilometers)
        station_dist = greatCircleDistance(np.radians(lat_centre), np.radians(lon_centre), stn.position.lat_r, stn.position.lon_r)

        # Construct time array, 0 is at start_datetime
        time_data = np.arange(0, (end_datetime - start_datetime).total_seconds(), delta)

        # Cut the waveform data length to match the time data
        waveform_data = waveform_data[:len(time_data)]
        time_data = time_data[:len(waveform_data)] + stn.offset
        
        # Detrend the waveform and normalize to fixed width
        waveform_data = waveform_data - np.mean(waveform_data)

        #waveform_data = waveform_data/np.percentile(waveform_data, 99)*2
        waveform_data = waveform_data/np.max(waveform_data)*10

        # Add the distance to the waveform
        waveform_data += station_dist


        # Cut the waveforms around the time of arrival, if the window for cutting was given.
        if waveform_window is not None:

            # Time of arrival
            toa = station_dist/(v_sound/1000) + t0

            # Cut the waveform around the time of arrival
            crop_indices = (time_data >= toa - waveform_window/2) & (time_data <= toa + waveform_window/2)
            time_data = time_data[crop_indices]
            waveform_data = waveform_data[crop_indices]
            

            # Skip plotting if array empty
            if len(time_data) == 0:
                continue

        # Replace all NaNs with 0s
        waveform_data = np.nan_to_num(waveform_data, 0)
        
        max_time = np.max([max_time, np.max(time_data)])
        min_time = np.min([min_time, np.min(time_data)])

        # Keep track of minimum and maximum waveform values (used for plotting)
        max_wave_value = np.max([max_wave_value, np.max(waveform_data)])
        min_wave_value = np.min([min_wave_value, np.min(waveform_data)])

        if setup.colortoggle:
            c = plt.cm.plasma(az[idx])
        else:
            c = None
            
        #if data_list[idx][1].strip() not in setup.rm_stat: 
            
        # Plot the waveform on the the time vs. distance graph
        ax.plot(waveform_data, time_data, c=c, alpha=0.7, linewidth=0.2, zorder=2)

        if stn.code in setup.rm_stat:
            print('Excluding station: {:}'.format(stn.network + '-' + stn.code))
        else:
            # Print the name of the station
            # Fragmentation
            if stn.code in setup.high_f:
                ax.text(np.mean(waveform_data), np.max(time_data), "{:} - {:} \n Az: {:5.1f}".format(stn.network, stn.code, az_n[idx]), \
                    rotation=270, va='bottom', ha='center', size=7, zorder=2, color="g")
            # Ballistic
            elif stn.code in setup.high_b:
                ax.text(np.mean(waveform_data), np.max(time_data), "{:} - {:} \n Az: {:5.1f}".format(stn.network, stn.code, az_n[idx]), \
                    rotation=270, va='bottom', ha='center', size=7, zorder=2, color="b")
        
            else:
                ax.text(np.mean(waveform_data), np.max(time_data), "{:} - {:} \n Az: {:5.1f}".format(stn.network, stn.code, az_n[idx]), \
                    rotation=270, va='bottom', ha='center', size=7, zorder=2, color="w")

    toa_line_time = np.linspace(0, max_time, 10)

    # Plot the constant sound speed line (assumption is that the release happened at t = 0)
    ax.plot((toa_line_time)*v_sound/1000, (toa_line_time + t0), color='r', alpha=0.25, linewidth=1, \
        zorder=2, label="$V_s = " + "{:d}".format(int(v_sound)) + r" \rm{ms^{-1}}$")

    # Reference location for the local coordinate system
    ref_pos = position(lat_centre, lon_centre, 0)
    # Ballistic Prediction
    b_time = [0]*len(stn_list)
    b_dist = [0]*len(stn_list)
    rb_dist = [0]*len(stn_list)

    good_stats = (x for x in range(len(stn_list)) if x not in bad_stats)
    
    print('')

    if setup.perturb_times <= 0 and setup.perturb:
        print("ERROR: perturb_times must be greater than 0")

    # for ptb_n in range(setup.perturb_times):
    #     # Manual search for ballistic wave
    #     if ptb_n > 0:
    #         print("STATUS: Perturbation: {:}".format(ptb_n))
    #         sounding_p = perturb(sounding, setup.perturb_method)
    #     else:
            # sounding_p = sounding
    sounding_p = sounding

    if setup.show_ballistic_waveform:
        # Input coordinate type. True - coordinates are given as lat/lon. False - coordinates are given in local 
        # coordinates in reference to the source center

        # Convert to local coordinates
        setup.traj_f.pos_loc(ref_pos)


        for stn in good_stats:

            if stn_list[stn].code.strip() not in setup.rm_stat: 
                # Station location in local coordinates

                stn_list[stn].position.pos_loc(ref_pos)

                # Time to travel from trajectory to station
                b_time[i] = timeOfArrival([stn_list[stn].position.x, stn_list[stn].position.y, stn_list[stn].position.z], setup.traj_f.x/1000, setup.traj_f.y/1000, setup.t0, 1000*setup.v, \
                                            np.radians(setup.azim), np.radians(setup.zangle), setup, sounding=sounding_p, fast=True)# - setup.t + setup.t0

                # Point on trajectory where wave is released
                bx, by, bz = waveReleasePoint([stn_list[stn].position.x, stn_list[stn].position.y, stn_list[stn].position.z], setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, \
                                            np.radians(setup.azim), np.radians(setup.zangle), setup.v_sound)

                # Distance from source center to station
                b_dist[i] = ((stn_list[stn].position.x)**2 + (stn_list[stn].position.y)**2)**0.5

                # Distance from ballistic wave to station
                rb_dist[i] = ((stn_list[stn].position.x - bx)**2 + (stn_list[stn].position.y - by)**2 + (stn_list[stn].position.z - bz)**2)**0.5

                # Convert to km
                b_dist[i] /= 1000
                rb_dist[i] /= 1000

            else:
                b_dist[i], b_time[i], rb_dist[i] = np.nan, np.nan, np.nan
            
    # Plot Ballistic Prediction
    #if ptb_n == 0:
    ax.scatter(b_dist, b_time, c='b', marker='_', s=100, label='Ballistic', zorder=3)
    #else:
    #    ax.scatter(b_dist, b_time, c='b', marker='_', s=100, alpha=0.3, zorder=3)

    # Fragmentation Prediction
    f_time = [0]*len(stn_list)
    f_dist = [0]*len(stn_list)
    rf_dist = [0]*len(stn_list)

    # Manual search for fragmentation waves
    if setup.show_fragmentation_waveform:

        if len(setup.fragmentation_point) == 0:
            print("ERROR: Cannot plot fragmentation if there is no fragmentation point. Set show_fragmentation_waveform = False if not using.")
            exit()

        for j, line in enumerate(setup.fragmentation_point):
            # Supracenter location in local coordinates
            supra = position(float(line[0]), float(line[1]), float(line[2]))

            supra.pos_loc(ref_pos)

            for i, stn in enumerate(stn_list):

                if stn.code.strip() not in setup.rm_stat:
                    if stn in bad_stats:
                        f_dist[i], f_time[i], rf_dist[i] = np.nan, np.nan, np.nan
                    # Station location in local coordinates

                    stn.position.pos_loc(ref_pos)

                    ###### DIFFERENT WEATHERS HERE ######
                    if setup.weather_type == 'none':
                        zProfile = np.array([[0, setup.v_sound, 0, 0], [10000, setup.v_sound, 0, 0]])

                    else:   
                    # Cut down atmospheric profile to the correct heights, and interp
                        zProfile, _ = getWeather(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), setup.weather_type, \
                            [ref_pos.lat, ref_pos.lon, ref_pos.elev], sounding_p, convert=True)

                    # Time to travel from Supracenter to station 
                    f_time[i], _, _ = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True)
                    
                    # Add reference time
                    f_time[i] += float(line[3])
                    
                    # Distance from source center to station
                    f_dist[i] = ((stn.position.x)**2 + (stn.position.y)**2)**0.5

                    # Distance from Supracenter to station
                    rf_dist[i] = ((stn.position.x - supra.x)**2 + (stn.position.y - supra.y)**2 + (stn.position.z - supra.z)**2)**0.5

                    # Convert to km
                    f_dist[i] /= 1000
                    rf_dist[i] /= 1000

                else:
                    f_dist[i], f_time[i], rf_dist[i] = np.nan, np.nan, np.nan
                    # Plot Fragmentation Prediction
            # if ptb_n == 0:
            ax.scatter(f_dist, f_time, c=C[(j+1)%4], marker='_', s=100, label='Fragmentation {:}'.format(j+1), zorder=3)
            # else:
                # ax.scatter(f_dist, f_time, c=C[(j+1)%4], marker='_', s=100, alpha=0.3, zorder=3)
            #ax.scatter(f_dist[i], f_time[i], c=C[(j+1)%4], marker='_', s=100, label='Fragmentation', zorder=3)


        
    ax.set_xlabel('Distance (km)')
    ax.set_ylabel('Time (s)')

    #ax.set_ylim(min_time - 200, max_time + 500)
    ax.set_xlim(0, max_wave_value)

    ax.grid(color='#ADD8E6', linestyle='dashed', linewidth=0.5, alpha=0.7)

    # Export station distance file
    with open(os.path.join(dir_path, 'output.txt'), 'w') as f:

        f.write('Station Lat(deg N) Lon(deg E) Elev(m) Az(+E dN) Ball_d(km) Ball_t(s) Frag_d(km) Frag_t(s)\n')
        for i, stn in enumerate(stn_list):
            f.write('{:8}, {:8.4f}, {:8.4f}, {:7.2f}, {:8.3f},  {:7.2f},  {:7.2f},  {:7.2f},  {:7.2f}\n'\
             .format(str(stn.network) + '-' + str(stn.code), stn.position.lat, stn.position.lon, \
                    stn.position.elev, az_n[i], rb_dist[i], b_time[i], rf_dist[i], f_time[i]))



if __name__ == "__main__":


    ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description="""
                ~~GetIRISData~~ 
    Plot seismic/infrasound data at distances
    away from a given source to show station
    arrival times of a signal

    Denis Vida, Luke McFadden
    """,
        formatter_class=argparse.RawTextHelpFormatter)

    arg_parser.add_argument('input_file', type=str, help='Path to Supracenter input file.')

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    #################

    setup = configRead(cml_args.input_file)
    configParse(setup, 'get_data')

    # 90 - zenith is used in the program
    #setup.zangle = 90 - setup.zangle

    # Create fireball folder
    if not os.path.exists(setup.working_directory):
        os.makedirs(setup.working_directory)

    #Build seismic data path
    dir_path = os.path.join(setup.working_directory, setup.fireball_name)

    ##########################################################################################################

    if setup.get_data:
        ### Download all waveform files which are within the given geographical and temporal range ###
        ##########################################################################################################
        getAllWaveformFiles(setup.lat_centre, setup.lon_centre, setup.deg_radius, setup.start_datetime, \
            setup.end_datetime, network='*', channel='all', dir_path=dir_path)
        ##########################################################################################################

    data_file_path = os.path.join(dir_path, DATA_FILE)

    if os.path.isfile(data_file_path):
        
        stn_list = readStationAndWaveformsListFile(data_file_path)

    else:
        print('Station and waveform data file not found! Download the waveform files first!')
        sys.exit()

    # Init the constants
    consts = Constants()

    setup.search_area = [0, 0, 0, 0]
    setup.search_area[0] = setup.lat_centre - setup.deg_radius 
    setup.search_area[1] = setup.lat_centre + setup.deg_radius
    setup.search_area[2] = setup.lon_centre - setup.deg_radius
    setup.search_area[3] = setup.lon_centre + setup.deg_radius

    sounding = parseWeather(setup, consts)

    ### Plot station map ###

    ##########################################################################################################
    if setup.plot_all_stations:
        isc_stn_list = setup.stations
    else:
        isc_stn_list = None

    plotStationMap(dir_path, stn_list, setup.lat_centre, setup.lon_centre, setup, sounding, \
        isc_data=isc_stn_list)

    #plt.legend(loc='upper left')
    plt.savefig(os.path.join(dir_path, "all_channels_{:s}_stations.png".format( \
        str(setup.start_datetime).replace(':', '.'))), dpi=300)

    plt.show()

    ##########################################################################################################

    ### Filter and plot all downloaded waveforms ###
    ##########################################################################################################
    stn_list = stn_list + isc_stn_list
    plotAllWaveforms(dir_path, stn_list, setup, sounding, difference_filter_all=setup.difference_filter_all)

    plt.title('Source location: {:.6f}, {:.6f}, Reference time: {:s} UTC, channel: {:s}'.format(setup.lat_centre, \
        setup.lon_centre, str(setup.start_datetime), '*'), fontsize=7)
    plt.legend(loc='lower right')

    plt.savefig(os.path.join(dir_path, "{:s}_waveforms.png".format(str(setup.start_datetime).replace(':', '.'))), dpi=300)


    plt.show()