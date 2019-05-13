""" Tool for marking times of arrival of seismic/air waves in IRIS data. """

from __future__ import print_function, division, absolute_import

import os
import sys
import datetime
import argparse

import obspy
import numpy as np
import scipy.signal
import time
import copy

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.widgets import Button
from matplotlib.widgets import Slider

# # Force backend
# plt.switch_backend("TkAgg")

DATA_FILE = 'data.txt'
OUTPUT_CSV = 'data_picks.csv'

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Fireballs.GetIRISData import readStationAndWaveformsListFile, butterworthBandpassFilter, \
    plotAllWaveforms, convolutionDifferenceFilter
from supra.Fireballs.SeismicTrajectory import timeOfArrival, waveReleasePoint, parseWeather, Constants
from supra.Fireballs.Program import configRead, configParse, position
from supra.Supracenter.angleConv import geo2Loc, loc2Geo
from supra.Supracenter.cyscan import cyscan
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.SPPT import perturb
from supra.Supracenter.weatherGraph import graphWeather
from wmpl.Utils.Earth import greatCircleDistance
from wmpl.Utils.PlotMap import GroundMap
from wmpl.Utils.TrajConversions import datetime2JD

global arrTimes 
global sounding


def trajCalc(setup):
    """ Creates trajectory between point A and the ground (B) based off of the initial position and the angle of travel

    Arguments:
    setup: [Object] ini file parameters

    Returns: 
    A [list] lat/lon/elev of the tail of the trajectory
    B [list] lat/lon/elev of the head of the trajectory
    """
    print("Trajectory Calculation...")

    ref = position(setup.lat_centre, setup.lon_centre, 0)

    # Tail of the trajectory
    A = geo2Loc(ref.lat, ref.lon, ref.elev, setup.lat_i, setup.lon_i, setup.elev_i*1000)

    # convert angles to radians
    ze = np.radians(setup.zangle)
    az = np.radians(setup.azim)

    # Create trajectory vector
    traj = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

    # How far along the trajectory until it reaches the ground
    n = -A[2]/traj[2]

    # B is the intersection between the trajectory vector and the ground
    B = A + n*traj

    # Convert back to geo coordinates
    B = np.array(loc2Geo(ref.lat, ref.lon, ref.elev, B))
    A = np.array(loc2Geo(ref.lat, ref.lon, ref.elev, A))

    print("Created Trajectory between A and B:")
    print("     A = {:10.4f}N {:10.4f}E {:10.2f}m".format(A[0], A[1], A[2]))
    print("     B = {:10.4f}N {:10.4f}E {:10.2f}m".format(B[0], B[1], B[2]))

    return A, B


def calcAllTimes(stn_list, setup, sounding):
        """ 
        Method to calculate all arrival times and place them into an array, so that they are not recalculated every
        time a new waveform is opened
                
        Arguments:
        data_list: [ndarray] array containing all station names and locations
        setup: [object] ini file parameters
        sounding: [ndarray] atmospheric profile

        Returns:
        allTimes: [ndarray] a ndarray that stores the arrival times for all stations over every perturbation
        [perturbation, station, 0 - ballistic/ 1 - fragmentation, frag number (0 for ballistic)]
        """

        #All perturbation happens here
        allTimes = [0]*setup.perturb_times

        # Ballistic Prediction
        ref_pos = position(setup.lat_centre, setup.lon_centre, 0)

        no_of_frags = len(setup.fragmentation_point)

        # array of frags and ballistic arrivals have to be the same size. So, minimum can be 1
        if no_of_frags == 0:
            no_of_frags = 1

        # Initialize variables
        b_time = 0

        consts = Constants()


        # For temporal perturbations, fetch the soudning data for the hour before and after the event
        if setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(setup, consts, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(setup, consts, time=-1)

        else:
            sounding_u = []
            sounding_l = []

        d_time = 2*(setup.perturb_times*len(stn_list)*no_of_frags)
        count = 0

        #number of perturbations
        for ptb_n in range(setup.perturb_times):

            if ptb_n > 0:
                
                if setup.debug:
                    print("STATUS: Perturbation {:}".format(ptb_n))

                # generate a perturbed sounding profile
                sounding_p = perturb(setup, sounding, setup.perturb_method, \
                    sounding_u=sounding_u, sounding_l=sounding_l, \
                    spread_file=setup.perturbation_spread_file, lat=setup.lat_centre, lon=setup.lon_centre)
            else:

                # if not using perturbations on this current step, then return the original sounding profile
                sounding_p = sounding

            # Initialize station times array
            stnTimes = [0]*len(stn_list)

            #number of stations
            for n, stn in enumerate(stn_list):

                # For ballistic arrivals
                if setup.show_ballistic_waveform:
                    bTimes = [0]*no_of_frags
                    for i in range(no_of_frags):
                        count += 1
                        sys.stdout.write("\rCalculating all times: {:5.2f} % ".format(count/d_time * 100))
                        sys.stdout.flush()
                        time.sleep(0.001)
                        #need filler values to make this a numpy array with fragmentation
                        if i == 0:
                            
                            stn.position.pos_loc(ref_pos)
                            setup.traj_f.pos_loc(ref_pos)

                            # Time to travel from trajectory to station
                            b_time = timeOfArrival([stn.position.x, stn.position.y, stn.position.z], setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, \
                                                        np.radians(setup.azim), np.radians(setup.zangle), setup, sounding=sounding_p, travel=False, fast=False, ref_loc=[ref_pos.lat, ref_pos.lon, ref_pos.elev])# + setup.t 


                            bTimes[i] = b_time
                        else:
                            bTimes[i] = np.nan
                else:
                    bTimes = [np.nan]*no_of_frags

                # Fragmentation Prediction
                f_time = np.array([0]*no_of_frags)

                # If manual fragmentation search is on
                if setup.show_fragmentation_waveform:
                    fTimes = [0]*no_of_frags
                    for i, line in enumerate(setup.fragmentation_point):
                        count += 1
                        sys.stdout.write("\rCalculating all times: {:5.2f} % ".format(count/d_time * 100))

                        # location of supracenter
                        supra = position(float(line[0]), float(line[1]), float(line[2]))
                        
                        # convert to local coordinates based off of the ref_pos
                        supra.pos_loc(ref_pos)

                        # convert station coordinates to local coordinates based on the ref_pos
                        stn.position.pos_loc(ref_pos)

                        # Cut down atmospheric profile to the correct heights, and interp
                        zProfile, _ = getWeather(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), setup.weather_type, \
                                [ref_pos.lat, ref_pos.lon, ref_pos.elev], copy.copy(sounding_p), convert=True)

                        # Travel time of the fragmentation wave
                        f_time, _, _ = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True, \
                            n_theta=setup.n_theta, n_phi=setup.n_phi, precision=setup.angle_precision, tol=setup.angle_error_tol)

                        fTimes[i] = f_time + line[3]

                else:

                    # Repack all arrays into allTimes array
                    fTimes = [np.nan]*no_of_frags

                stnTimes[n] = ([np.array(bTimes), np.array(fTimes)])

            allTimes[ptb_n] = np.array(stnTimes)

        allTimes = np.array(allTimes)

        # Save as .npy file to be reused in SeismicTrajectory and Supracenter
        np.save(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times'), allTimes)
        print("All Times File saved as {:}".format(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times.npy')))

        return allTimes


class WaveformPicker(object):
    def __init__(self, dir_path, setup, sounding, data_list, waveform_window=600, \
        difference_filter_all=False, stn_list=[]):
        """

        Arguments:
            data_list: [list]

        Keyword arguments:
            waveform_window: [int] Number of seconds for the wavefrom window.
            difference_filter_all: [bool] If True, the Kalenda et al. (2014) difference filter will be applied
                on the data plotted in the overview plot of all waveforms.
        """
        self.setup = setup
        self.dir_path = dir_path

        self.v_sound = setup.v_sound
        self.t0 = setup.t0

        self.stn_list = data_list

        # Filter out all stations for which the mseed file does not exist
        filtered_stn_list = []

        names = []
        lats = []
        lons = []

        for stn in self.stn_list:
            
            mseed_file = stn.file_name
            mseed_file_path = os.path.join(self.dir_path, mseed_file)

            if os.path.isfile(mseed_file_path):
                filtered_stn_list.append(stn)

            else:
                print('mseed file does not exist:', mseed_file_path)

        self.stn_list = filtered_stn_list


        self.lat_centre = setup.lat_centre
        self.lon_centre = setup.lon_centre

        self.waveform_window = waveform_window


        self.current_station = 0
        self.current_wavefrom_raw = None
        self.current_wavefrom_delta = None
        self.current_waveform_processed = None

        # List of picks
        self.pick_list = []

        self.pick_group = 0

        # Define a list of colors for groups
        self.pick_group_colors = ['r', 'g', 'm', 'k', 'y']

        # Current station map handle
        self.current_station_scat = None

        # Station waveform marker handle
        self.current_station_all_markers = None

        # Picks on all waveform plot handle
        self.all_waves_picks_handle = None

        # Handle for pick text
        self.pick_text_handle = None
        self.pick_markers_handles = []

        # handle for pick marker on the wavefrom
        self.pick_wavefrom_handle = None


        # Default bandpass values
        self.bandpass_low_default = 2.0
        self.bandpass_high_default = 8.0

        # Flag indicating whether CTRL is pressed or not
        self.ctrl_pressed = False


        ### Sort stations by distance from source ###

        # Calculate distances of station from source
        self.source_dists = []


        for stn in self.stn_list:

            stat_name, stat_lat, stat_lon = stn.code, stn.position.lat, stn.position.lon

            names.append(stat_name)
            lats.append(stat_lat)
            lons.append(stat_lon)

            # Calculate the distance in kilometers
            dist = greatCircleDistance(np.radians(setup.lat_centre), np.radians(setup.lon_centre), \
                np.radians(stat_lat), np.radians(stat_lon))

            self.source_dists.append(dist)

        # Get sorted arguments
        dist_sorted_args = np.argsort(self.source_dists)

        # Sort the stations by distance
        self.stn_list = [self.stn_list[i] for i in dist_sorted_args]
        self.source_dists = [self.source_dists[i] for i in dist_sorted_args]
        #############################################

        # Init the plot framework
        self.initPlot(setup, sounding)

        # Extract the list of station locations
        self.lat_list = [stn.position.lat_r for stn in self.stn_list]
        self.lon_list = [stn.position.lon_r for stn in self.stn_list]

        # Init ground map
        self.m = GroundMap(self.lat_list, self.lon_list, ax=self.ax_map, color_scheme='light')

        for stn in self.stn_list:

            # Plot stations
            if stn.code in setup.high_f:
                self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='g', s=2)
            elif stn.code in setup.high_b:
                self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='b', s=2)
            else:
                self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='k', s=2)

        # Manual Supracenter search
        if setup.show_fragmentation_waveform:
            
            # Fragmentation plot
            for i, line in enumerate(setup.fragmentation_point):
                self.m.scatter([np.radians(float(line[0]))], [np.radians(float(line[1]))], c=self.pick_group_colors[(i+1)%4], marker='x')

        # Extract coordinates of the reference station
        ref_pos = position(setup.lat_centre, setup.lon_centre, 0)

        # Plot source location
        self.m.scatter([np.radians(setup.lat_centre)], [np.radians(setup.lon_centre)], marker='*', c='yellow')

        # Manual trajectory search
        if setup.show_ballistic_waveform:

            # Plot the trajectory with the bottom point known
            self.m.plot([setup.traj_i.lat_r, setup.traj_f.lat_r], [setup.traj_i.lon_r, setup.traj_f.lon_r], c='b')
            # Plot intersection with the ground
            self.m.scatter(setup.traj_f.lat_r, setup.traj_f.lon_r, s=10, marker='x', c='b')

            ### CONTOUR ###

            # Get the limits of the plot
            x_min = setup.lat_f - 100000*setup.deg_radius
            x_max = setup.lat_f + 100000*setup.deg_radius
            y_min = setup.lon_f - 100000*setup.deg_radius
            y_max = setup.lon_f + 100000*setup.deg_radius

            img_dim = int(setup.contour_res)
            x_data = np.linspace(x_min, x_max, img_dim)
            y_data = np.linspace(y_min, y_max, img_dim)
            xx, yy = np.meshgrid(x_data, y_data)


            # # Make an array of all plane coordinates
            plane_coordinates = np.c_[xx.ravel(), yy.ravel(), np.zeros_like(xx.ravel())]

            times_of_arrival = np.zeros_like(xx.ravel())

            az = np.radians(setup.azim)
            ze = np.radians(setup.zangle)

            # vector of the trajectory of the fireball
            traj_vect = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

            for i, plane_coords in enumerate(plane_coordinates):

                # Print out percentage complete
                if (i + 1) % 10 == 0:
                    sys.stdout.write("\rDrawing Contour: {:.2f} %".format(100*(i + 1)/img_dim**2))
                    sys.stdout.flush()
                    time.sleep(0.001)

                setup.traj_f.pos_loc(ref_pos)
                # Point on the trajectory where the sound wave that will hit the plane_coord originated from

                p = waveReleasePoint(plane_coords, setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, az, \
                                          ze, setup.v_sound)

                # # vector between the wave release point and the plane coordinate
                d_vect = plane_coords - p

                # Since the arrivals are always perpendicular to the fireball trajectory, only take arrivals where the dot product
                # of the vectors are small.

                if abs(np.dot(d_vect/1000, traj_vect)) < setup.dot_tol:

                    # time of arrival from the trajectory
                    ti = timeOfArrival(plane_coords, setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, \
                                       az, ze, setup, sounding=sounding, ref_loc=[ref_pos.lat, ref_pos.lon, 0], travel=True, fast=True)# - setup.t + setup.t0

                # escape value for if sound never reaches the plane_coord
                else:
                   ti = np.nan

                times_of_arrival[i] = ti + setup.t0

            print('')

            # if sound never reaches the plane_coord, set to maximum value of the contour
            max_time = np.nanmax(times_of_arrival)
            for i in range(len(times_of_arrival)):
                if np.isnan(times_of_arrival[i]):
                    times_of_arrival[i] = max_time

            times_of_arrival = times_of_arrival.reshape(img_dim, img_dim)

            # Determine range and number of contour levels, so they are always centred around 0
            toa_abs_max = np.max([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
            #  toa_abs_min = np.min([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
            levels = np.linspace(0, toa_abs_max, 25)

            # Convert contour local coordinated to geo coordinates
            lat_cont = []
            lon_cont = []

            for x_cont, y_cont in zip(xx.ravel(), yy.ravel()):
                
                lat_c, lon_c, _ = loc2Geo(ref_pos.lat, ref_pos.lon, ref_pos.elev, np.array([x_cont, y_cont, 0]))

                lat_cont.append(lat_c)
                lon_cont.append(lon_c)

            lat_cont = np.array(lat_cont).reshape(img_dim, img_dim)
            lon_cont = np.array(lon_cont).reshape(img_dim, img_dim)

            # Plot the time of arrival contours
            toa_conture = self.m.m.contourf(lon_cont, lat_cont, times_of_arrival, levels, zorder=3, \
                latlon=True, cmap='viridis_r', alpha=0.5)

            # # Add a color bar which maps values to colors
            self.m.m.colorbar(toa_conture, label='Time of arrival (s)')

        if setup.arrival_times_file != '':
            try:
                self.arrTimes = np.load(setup.arrival_times_file)
                print("Reading in arrival times file...")
            except:
                print("WARNING: Unable to load allTimes_file {:} . Please check that file exists".format(setup.arrival_times_file))
                self.arrTimes = calcAllTimes(self.stn_list, setup, sounding)
        else:  
            # Calculate all arrival times
            self.arrTimes = calcAllTimes(self.stn_list, setup, sounding)
        
        self.updatePlot(setup)

    

    def initPlot(self, setup, sounding):
        """ Initializes the plot framework. """


        ### Init the basic grid ###

        gs = gridspec.GridSpec(4, 3, width_ratios=[1, 1, 1], height_ratios=[10, 5, 1, 1])


        # All waveform axis
        self.ax_all_waves = plt.subplot(gs[0, 0:2])

        # Map axis
        self.ax_map = plt.subplot(gs[0, 2])

        # Waveform axis
        self.ax_wave = plt.subplot(gs[1, :])

        # Register a mouse press event on the waveform axis
        plt.gca().figure.canvas.mpl_connect('button_press_event', self.onWaveMousePress)

        # Init what happes when a keyboard key is pressed or released
        plt.gca().figure.canvas.mpl_connect('key_press_event', self.onKeyPress)
        plt.gca().figure.canvas.mpl_connect('key_release_event', self.onKeyRelease)

        # Register window resize
        plt.gca().figure.canvas.mpl_connect('resize_event', self.onResize)

        # Previous button axis
        self.ax_prev_btn = plt.subplot(gs[2, 0])

        # Next button axis
        self.ax_next_btn = plt.subplot(gs[3, 0])

        # Bandpass options
        bandpass_gridspec = gridspec.GridSpecFromSubplotSpec(5, 1, subplot_spec=gs[2:4, 1])
        self.ax_bandpass_low = plt.subplot(bandpass_gridspec[0])
        self.ax_bandpass_high = plt.subplot(bandpass_gridspec[1])
        self.ax_bandpass_button = plt.subplot(bandpass_gridspec[2])

        # Spectrogram button
        self.ax_specgram_btn = plt.subplot(bandpass_gridspec[3])

        # Convolution filter button
        self.ax_convolution_button = plt.subplot(bandpass_gridspec[4])
        

        # Pick list
        picks_gridspec = gridspec.GridSpecFromSubplotSpec(4, 1, subplot_spec=gs[2:, 2])
        self.ax_picks = plt.subplot(picks_gridspec[:3])

        # Disable ticks and axis lines
        self.ax_picks.set_axis_off()

        # Export CSV button
        self.ax_export_csv_btn = plt.subplot(picks_gridspec[3])


        ############################


        # Create 'prev' button
        self.prev_btn = Button(self.ax_prev_btn, 'Previous')
        self.prev_btn.on_clicked(self.decrementStation)

        # Create 'next' button
        self.next_btn = Button(self.ax_next_btn, 'Next')
        self.next_btn.on_clicked(self.incrementStation)


        # Bandpass sliders and button
        self.bandpass_low_slider = Slider(self.ax_bandpass_low, 'Low:', 0.1, 5, \
            valinit=self.bandpass_low_default)
        self.bandpass_high_slider = Slider(self.ax_bandpass_high, 'High:', 3, 40, \
            valinit=self.bandpass_high_default, slidermin=self.bandpass_low_slider)
        
        self.bandpass_button = Button(self.ax_bandpass_button, 'Bandpass filter')
        self.bandpass_button.on_clicked(self.filterBandpass)


        # Spectrogram button
        self.specgram_btn = Button(self.ax_specgram_btn, 'Spectrogram of raw data')
        self.specgram_btn.on_clicked(self.showSpectrogram)

        # Convolution filter button
        self.convolution_button = Button(self.ax_convolution_button, '$S_i - (S_{i-1} + S_{i+1})/2$ filter')
        self.convolution_button.on_clicked(self.filterConvolution)


        # Export CSV button
        self.export_csv_btn = Button(self.ax_export_csv_btn, 'Export CSV: ' + OUTPUT_CSV)
        self.export_csv_btn.on_clicked(self.exportCSV)
    
        # Plot all waveforms
        plotAllWaveforms(self.dir_path, list(self.stn_list), setup, sounding, ax=self.ax_all_waves, \
            waveform_window=self.waveform_window, difference_filter_all=setup.difference_filter_all)

    def onKeyPress(self, event):
        
        if event.key == 'control':
            self.ctrl_pressed = True


        elif event.key == '+':
            
            # Increment the pick group
            self.pick_group += 1

            self.updatePlot()


        elif event.key == '-':
            # Decrement the pick group

            if self.pick_group > 0:
                self.pick_group -= 1

            self.updatePlot()


    def onKeyRelease(self, event):

        if event.key == 'control':
            self.ctrl_pressed = False



    def onResize(self, event):

        # Perform tight layout when window is resized
        plt.tight_layout()




    def addPick(self, pick_group, station_no, pick_time):
        """ Adds the pick to the list of picks. """

        self.pick_list.append([pick_group, station_no, pick_time])

        self.updatePickList()


    def removePick(self, station_no, pick_time_remove):
        """ Removes the pick from the list of picks with the closest time. """


        if len(self.pick_list):

            closest_pick_indx = None
            min_time_diff = np.inf

            # Go though all stations and find the pick with closest to the given time
            for i, entry in enumerate(self.pick_list):

                pick_grp, stat_no, pick_time = entry

                # Check if current station
                if stat_no == self.current_station:

                    time_diff = abs(pick_time - pick_time_remove)

                    # Store the minimum time difference
                    if time_diff < min_time_diff:

                        min_time_diff = time_diff
                        closest_pick_indx = i



            if closest_pick_indx is not None:
                
                # Remove pick on the given station closest the given time
                self.pick_list.pop(closest_pick_indx)


                self.updatePickList()




    def updatePickList(self):
        """ Updates the list of picks on the screen for the given station. """

        stations_with_picks = []
        all_pick_times = []

        for entry in self.pick_list:

            pick_grp, station_no, pick_time = entry

            stations_with_picks.append(station_no)
            all_pick_times.append(pick_time)


        # Remove old picks on all wavefrom plot
        if self.all_waves_picks_handle is not None:
            self.all_waves_picks_handle.remove()


        # Get distances of of pick stations
        dists_with_picks = [self.source_dists[stat_no] for stat_no in stations_with_picks]

        # Mark picks on the all waveform plot
        self.all_waves_picks_handle = self.ax_all_waves.scatter(dists_with_picks, all_pick_times, \
            marker='*', s=50, c='r')

        self.updatePlot(draw_waveform=False)


    def updatePickTextAndWaveMarker(self):
        """ Updates the list of picks on the screen. """


        current_station_groups = []
        current_station_picks = []

        for entry in self.pick_list:

            pick_grp, station_no, pick_time = entry

            # Take picks taken on the current station
            if station_no == self.current_station:
                current_station_groups.append(pick_grp)
                current_station_picks.append(pick_time)

        # Remove old pick text
        if self.pick_text_handle is not None:
            self.pick_text_handle.remove()

        # Generate the pick string
        pick_txt_str  = 'Change group: +/-\n'
        pick_txt_str += 'Add/remove pick: CTRL + left/right click \n'
        pick_txt_str += '\n'
        pick_txt_str += 'Current group: {:5d}\n\n'.format(self.pick_group)
        pick_txt_str += 'Picks: Group, Time\n'
        pick_txt_str += '------\n'
        pick_txt_str += "\n".join(["{:5d},   {:.2f}".format(gr, pt) for gr, pt in zip(current_station_groups, \
            current_station_picks)])

        # Print picks on screen
        self.pick_text_handle = self.ax_picks.text(0, 1, pick_txt_str, va='top', fontsize=7)


        # Remove old pick markers
        for handle in self.pick_markers_handles:
            try:
                handle.remove()
            except:
                pass


        self.pick_markers_handles = []


        if len(current_station_picks) > 0:

            # Get a list of colors per groups
            color_list = [self.pick_group_colors[grp%len(self.pick_group_colors)] \
                for grp in current_station_groups]


            # Get the Y coordinate of the pick in the waveform plot
            if self.current_waveform_processed is not None:
                pick_y_list = []
                for pick_time in current_station_picks:
                    
                    pick_indx = np.abs(self.current_waveform_time - pick_time).argmin()
                    pick_y = self.current_waveform_processed[pick_indx]
                    pick_y_list.append(pick_y)

            else:
                pick_y_list = [0]*len(current_station_picks)

            # Set pick marker on the current wavefrom
            scat_handle = self.ax_wave.scatter(current_station_picks, pick_y_list, marker='*', \
                c=color_list, s=50)

            self.pick_markers_handles.append(scat_handle)

            # Plot group numbers above picks
            #self.pick_wavefrom_text_handles = []
            for c, grp, pt in zip(color_list, current_station_groups, current_station_picks):
                txt_handle = self.ax_wave.text(pt, 0, str(grp), color=c, ha='center', va='bottom')

                self.pick_markers_handles.append(txt_handle)


    def onWaveMousePress(self, event):

        # Check if the mouse was pressed within the waveform axis
        if event.inaxes == self.ax_wave:

            # Check if CTRL is pressed
            if self.ctrl_pressed:

                pick_time = event.xdata

                # Check if left button was pressed
                if event.button == 1:

                    # Extract network and station code
                    net, station_code = self.stn_list[self.current_station].network, self.stn_list[self.current_station].code

                    print('Adding pick on station {:s} at {:.2f}'.format(net + ": " + station_code, \
                        pick_time))

                    self.addPick(self.pick_group, self.current_station, pick_time)


                # Check if right button was pressed
                elif event.button == 3:
                    print('Removing pick...')

                    self.removePick(self.current_station, pick_time)





    def incrementStation(self, event=None):
        """ Increments the current station index. """

        self.current_station += 1

        if self.current_station >= len(self.stn_list):
            self.current_station = 0

        while self.checkExists() == False:
            self.current_station += 1
            if self.current_station >= len(self.stn_list):
                self.current_station = 0

        self.updatePlot()


    def decrementStation(self, event=None):
        """ Decrements the current station index. """

        self.current_station -= 1

        if self.current_station < 0:
            self.current_station = len(self.stn_list) - 1

        while self.checkExists() == False:
            self.current_station -= 1
            if self.current_station < 0:
                self.current_station = len(self.stn_list) - 1

        self.updatePlot()



    def markCurrentStation(self):
        """ Mark the position of the current station on the map. """

        # Extract current station
        stn = self.stn_list[self.current_station]

        if self.current_station_scat is None:

            # Mark the current station on the map
            self.current_station_scat = self.m.scatter([stn.position.lat_r], [stn.position.lon_r], s=20, \
                edgecolors='r', facecolors='none')

        else:

            # Calculate map coordinates
            stat_x, stat_y = self.m.m(stn.position.lon, stn.position.lat)

            # Set the new position
            self.current_station_scat.set_offsets([stat_x, stat_y])


    def checkExists(self):
        """
        Checks if the current waveform is readable
        """

        # Extract current station
        stn = self.stn_list[self.current_station]

        # Get the miniSEED file path
        mseed_file_path = os.path.join(self.dir_path, stn.file_name)

        try:
            
            if os.path.isfile(mseed_file_path):
                pass
            else:
                print('File {:s} does not exist!'.format(mseed_file_path))
                return False

        except TypeError as e:
            
            print('Opening file {:s} failed with error: {:s}'.format(mseed_file_path, str(e)))
            return False

        try:
            obspy.read(mseed_file_path)

        except TypeError:
            print('mseed file could not be read:', mseed_file_path)
            return False

        return True

    def drawWaveform(self, waveform_data=None):
        """ Draws the current waveform from the current station in the wavefrom window. Custom wavefrom 
            can be given an drawn, which is used when bandpass filtering is performed. 

        """
        setup = self.setup

        # Clear waveform axis
        self.ax_wave.cla()

        # Extract current station
        stn = self.stn_list[self.current_station]

        # Get the miniSEED file path
        mseed_file_path = os.path.join(self.dir_path, stn.file_name)

        # Try reading the mseed file, if it doesn't work, skip to the next frame
        try:
            mseed = obspy.read(mseed_file_path)

        except TypeError:
            print('mseed file could not be read:', mseed_file_path)
            #self.incrementStation()
            return None

        # Unpact miniSEED data
        delta = mseed[0].stats.delta
        start_time = mseed[0].stats.starttime
        end_time = mseed[0].stats.endtime
        

        # Check if the waveform data is already given or not
        if waveform_data is None:
            waveform_data = mseed[0].data

            # Store raw data for bookkeeping on first open
            self.current_wavefrom_raw = waveform_data


        # Convert the beginning and the end time to datetime objects
        start_datetime = start_time.datetime
        end_datetime = end_time.datetime

        self.current_wavefrom_delta = delta
        self.current_waveform_time = np.arange(0, (end_datetime - start_datetime).total_seconds() + delta, \
            delta)


        # ### BANDPASS FILTERING ###

        # # Init the butterworth bandpass filter
        # butter_b, butter_a = butterworthBandpassFilter(self.bandpass_low_default, \
        #     self.bandpass_high_default, 1.0/delta, order=6)

        # # Filter the data
        # waveform_data = scipy.signal.filtfilt(butter_b, butter_a, waveform_data)

        # ##########################


        # Construct time array, 0 is at start_datetime
        time_data = np.copy(self.current_waveform_time)

        # Cut the waveform data length to match the time data
        waveform_data = waveform_data[:len(time_data)]
        time_data = time_data[:len(waveform_data)]


        # Trim the ends to avoid issues with differenced data
        time_data = time_data[1:-1]
        waveform_data = waveform_data[1:-1]

        # Store currently plotted waveform
        self.current_waveform_processed = waveform_data

        # Calculate the time of arrival assuming constant propagation with the given speed of sound
        t_arrival = self.source_dists[self.current_station]/(self.v_sound/1000) + self.t0
        #t_arrival = 0

        # Calculate the limits of the plot to be within the given window limit
        time_win_min = t_arrival - self.waveform_window/2
        time_win_max = t_arrival + self.waveform_window/2

        # Plot the wavefrom
        self.ax_wave.plot(time_data, waveform_data, color='k', linewidth=0.2, zorder=3)

        # Initialize variables
        b_time = 0

        print('####################')
        print("Current Station: {:}".format(stn.name))
        print("Channel: {:}".format(stn.channel))
        # If manual ballistic search is on
        if setup.show_ballistic_waveform:

            # Plot Ballistic Prediction
            b_time = self.arrTimes[0, self.current_station, 0, 0]
            
            # check if nan
            if b_time == b_time:
                self.ax_wave.plot([b_time]*2, [np.min(waveform_data), np.max(waveform_data)], c='b', label='Ballistic', zorder=3)
                
                print("Ballistic Arrival: {:.3f} s".format(b_time))
            else:
                print("No Ballistic Arrival")

            for i in range(setup.perturb_times):
                if i >= 1:
                    try:
                        self.ax_wave.plot([self.arrTimes[i, self.current_station, 0, 0]]*2, \
                         [np.min(waveform_data), np.max(waveform_data)], alpha=0.3, c='b', zorder=3)
                    except:
                        pass
        # Fragmentation Prediction

        # If manual fragmentation search is on
        if setup.show_fragmentation_waveform:

            for i, line in enumerate(setup.fragmentation_point):

                f_time = self.arrTimes[0, self.current_station, 1, i]
            #     # check if nan
                if f_time == f_time:
                    # Plot Fragmentation Prediction
                    self.ax_wave.plot([f_time]*2, [np.min(waveform_data), np.max(waveform_data)], c=self.pick_group_colors[(i+1)%4], label='Fragmentation', zorder=3)
                    
                    if len(setup.fragmentation_point) > 1:
                        #self.ax_wave.text(f_time, np.min(waveform_data), 'Frag{:}'.format(i+1))
                        self.ax_wave.text(f_time, np.min(waveform_data) + int(i)/(len(setup.fragmentation_point))*(np.max(waveform_data) - np.min(waveform_data)), '{:.1f} km'.format(line[2]/1000))
                
                    print('Fragmentation {:} Arrival: {:.3f} s'.format(i+1, f_time))

                else:
                    print('No Fragmentation {:} Arrival'.format(i+1))

                for j in range(setup.perturb_times):
                    if j >= 1:
                        try:
                            self.ax_wave.plot([self.arrTimes[j, self.current_station, 1, i]]*2, [np.min(waveform_data),\
                                 np.max(waveform_data)], alpha=0.3,\
                                 c=self.pick_group_colors[(i+1)%4], zorder=3)
                        except:
                            pass

        # Set the time limits to be within the given window
        self.ax_wave.set_xlim(time_win_min, time_win_max)

        self.ax_wave.grid(color='#ADD8E6', linestyle='dashed', linewidth=0.5, alpha=0.5)

        #self.ax_wave.legend()

        # Add text with station label
        if stn.code in setup.high_f:
            self.ax_wave.text(time_win_min, np.max(waveform_data), stn.network + ": " + stn.code \
                + "(" + stn.channel + ")" +", {:d} km".format(int(self.source_dists[self.current_station])) , va='top', ha='left', color='g')

        elif stn.code in setup.high_b:
            self.ax_wave.text(time_win_min, np.max(waveform_data), stn.network + ": " + stn.code \
                + "(" + stn.channel + ")" + ", {:d} km".format(int(self.source_dists[self.current_station])) , va='top', ha='left', color='b')

        else:
            self.ax_wave.text(time_win_min, np.max(waveform_data), stn.network + ": " + stn.code \
                + "(" + stn.channel + ")" + ", {:d} km".format(int(self.source_dists[self.current_station])) , va='top', ha='left', color='k')


    def markStationWaveform(self):
        """ Mark the currently shown waveform in the plot of all waveform. """

        if self.current_station_all_markers is not None:
            for marker in self.current_station_all_markers:
                marker.remove()


        # Calculate the position
        dist = self.source_dists[self.current_station]

        # Calcualte the time of arrival
        t_arrival = self.source_dists[self.current_station]/(self.v_sound/1000) + self.t0

        # Plot the marker
        marker1 = self.ax_all_waves.scatter(dist, t_arrival - 500, marker='^', s=100, linewidths=3, c='k', 
            alpha=1, zorder=3)

        marker2 = self.ax_all_waves.scatter(dist, t_arrival + 500, marker='v', s=100, linewidths=3, c='k', 
            alpha=1, zorder=3)


        self.current_station_all_markers = [marker1, marker2]


    def showSpectrogram(self, event):
        """ Show the spectrogram of the waveform in the current window. """


        # Get time limits of the shown waveform
        x_min, x_max = self.ax_wave.get_xlim()

        # Extract the time and waveform
        crop_window = (self.current_waveform_time >= x_min) & (self.current_waveform_time <= x_max)
        wave_arr = self.current_wavefrom_raw[crop_window]


        ### Show the spectrogram ###
        
        fig = plt.figure()
        ax_spec = fig.add_subplot(111)

        ax_spec.specgram(wave_arr, Fs=1.0/self.current_wavefrom_delta, cmap=plt.cm.inferno, \
            xextent=(x_min, x_max))

        ax_spec.set_xlabel('Time (s)')
        ax_spec.set_ylabel('Frequency (Hz)')

        fig.show()

        ###


    def filterBandpass(self, event):
        """ Run bandpass filtering using values set on sliders. """

        # Get bandpass filter values
        bandpass_low = self.bandpass_low_slider.val
        bandpass_high = self.bandpass_high_slider.val


        # Limit the high frequency to be lower than the Nyquist frequency
        max_freq = (1.0/self.current_wavefrom_delta)/2

        if bandpass_high > max_freq:
            bandpass_high = max_freq - 0.1

            self.bandpass_high_slider.set_val(bandpass_high)
        

        # Init the butterworth bandpass filter
        butter_b, butter_a = butterworthBandpassFilter(bandpass_low, bandpass_high, \
            1.0/self.current_wavefrom_delta, order=6)

        # Filter the data
        waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(self.current_wavefrom_raw))


        # Plot the updated waveform
        self.drawWaveform(waveform_data)


    def filterConvolution(self, event):
        """ Apply the convolution filter on data as suggested in Kalenda et al. (2014). """

        waveform_data = convolutionDifferenceFilter(self.current_wavefrom_raw)

        self.drawWaveform(waveform_data)


    def updatePlot(self, draw_waveform=True):
        """ Update the plot after changes. """

        # Mark the position of the current station on the map
        self.markCurrentStation()

        # Plot the wavefrom from the current station
        if draw_waveform:
            self.drawWaveform()

        # Set an arrow pointing to the current station on the waveform
        self.markStationWaveform()

        # Update the pick list text and plot marker on the waveform
        self.updatePickTextAndWaveMarker()

        # Reset bandpass filter values to default
        self.bandpass_low_slider.set_val(self.bandpass_low_default)
        self.bandpass_high_slider.set_val(self.bandpass_high_default)

        plt.draw()
        plt.pause(0.001)


    def exportCSV(self, event):
        """ Save picks to a CSV file. """

        # Open the output CSV
        with open(os.path.join(self.dir_path, OUTPUT_CSV), 'w') as f:

            # Write the header
            f.write('Pick group, Network, Code, Lat, Lon, Elev, Pick JD, Pick time, station_number \n')

            # Go through all picks
            for entry in self.pick_list:

                # Unpack pick data
                pick_group, station_no, pick_time = entry

                # # Extract current station
                stn = self.stn_list[station_no]

                # # Unpack the station entry
                # net, station_code, stat_lat, stat_lon, stat_elev, station_name, channel, mseed_file = stat_entry

                # Get the miniSEED file path
                mseed_file_path = os.path.join(self.dir_path, stn.file_name)

                try:
                    
                    if os.path.isfile(mseed_file_path):
                        
                        # Read the miniSEED file
                        mseed = obspy.read(mseed_file_path)

                    else:
                        print('File {:s} does not exist!'.format(mseed_file_path))
                        continue

                except TypeError as e:

                    print('Opening file {:s} failed with error: {:s}'.format(mseed_file_path, e))
                    continue

                # Find datetime of the beginning of the file
                start_datetime = mseed[0].stats.starttime.datetime

                # Calculate Julian date of the pick time
                pick_jd = datetime2JD(start_datetime + datetime.timedelta(seconds=pick_time))


                # Write the CSV entry
                f.write("{:d}, {:s}, {:s}, {:.6f}, {:.6f}, {:.2f}, {:.8f}, {:}, {:}\n".format(pick_group, stn.network, \
                    stn.code, stn.position.lat, stn.position.lon, stn.position.elev, pick_jd, pick_time, station_no))


        print('CSV written to:', OUTPUT_CSV)


if __name__ == "__main__":  

        ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description="""
            ~~MakeIRISPicks~~ 
    Tool for marking times of arrival of 
    seismic/infrasound data in IRIS data. 

    Denis Vida
    """,
        formatter_class=argparse.RawTextHelpFormatter)

    arg_parser.add_argument('input_file', type=str, help='Path to Supracenter input file.')

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    #################

    setup = configRead(cml_args.input_file)
    configParse(setup, 'picks')

    ##########################################################################################################
    if not os.path.exists(setup.working_directory):
        os.makedirs(setup.working_directory)

        #Build seismic data path
    dir_path = os.path.join(setup.working_directory, setup.fireball_name)

    # Load the station and waveform files list
    data_file_path = os.path.join(dir_path, DATA_FILE)
    if os.path.isfile(data_file_path):
        
        stn_list = readStationAndWaveformsListFile(data_file_path, rm_stat=setup.rm_stat)

    else:
        print('Station and waveform data file not found! Download the waveform files first!')
        sys.exit()

    if setup.stations is not None:
        stn_list = stn_list + setup.stations

    # Init the constants
    consts = Constants()
    setup.search_area = [0, 0, 0, 0]
    setup.search_area[0] = setup.lat_centre - setup.deg_radius 
    setup.search_area[1] = setup.lat_centre + setup.deg_radius
    setup.search_area[2] = setup.lon_centre - setup.deg_radius
    setup.search_area[3] = setup.lon_centre + setup.deg_radius

    sounding = parseWeather(setup, consts)

    if setup.weather_type != 'none':
        S = position(setup.lat_centre, setup.lon_centre, 45000)
        D = position(setup.lat_centre, setup.lon_centre,     0)
        try:
            graphWeather(setup.weather_name, setup.weather_type, [S.lat, S.lon, S.elev], [D.lat, D.lon, D.elev], 'spread_r', 100, spread_file=setup.spread_file, setup=setup)
        except:
            pass
    arrTimes = np.array([])

    if len(stn_list) == 0:
        print('ERROR: No stations!')
        exit()

    # Init the wave\from picker
    WaveformPicker(dir_path, setup, sounding, stn_list, difference_filter_all=setup.difference_filter_all)

    plt.tight_layout()

    #plt.waitforbuttonpress(timeout=-1)

    plt.show()