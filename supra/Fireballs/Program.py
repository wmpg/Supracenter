import numpy as np
import argparse
import os
import time
import datetime
import matplotlib.pyplot as plt

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

try:
    # Python 2  
    import ConfigParser as configparser

except:
    # Python 3
    import configparser

from supra.Supracenter.angleConv import angle2NDE, loc2Geo, geo2Loc
#from supra.Fireballs.MakeIRISPicks import trajCalc
# from supra.Fireballs.GetIRISData import getAllWaveformFiles, readStationAndWaveformsListFile, plotStationMap, plotAllWaveforms
# from supra.Fireballs.MakeIRISPicks import WaveformPicker
#from supra.Fireballs.SeismicTrajectory import latLon2Local, local2LatLon
# from supra.Supracenter.mainmenu import main

DATA_FILE = 'data.txt'

class position:
    """ A latitude/longitude/elevation object describing a position in 3D space
    Give lat [-90, 90], lon [-180, 180] in degrees
    Give elev in m 
    """

    def __init__(self, lat, lon, elev):

        self.lat  = lat
        self.lon  = lon
        self.elev = elev
        self.lat_r = np.radians(lat)
        self.lon_r = np.radians(lon) 
    
    def pos_loc(self, ref_pos):
        """
        Converts a position object to local coordinates based off of a reference position

        Arguments:
            ref_pos: [position Obj] the position of a reference location 
            used to convert another position object to local coordinates
        """
        self.x, self.y, self.z = geo2Loc(ref_pos.lat, ref_pos.lon, ref_pos.elev, self.lat, self.lon, self.elev)

class station:

    """
    A station object containing information and position of the station
    """

    def __init__(self, network, code, position, channel, name, file_name):
        """
        Arguments:
            network: [String] the network of the station
            code: [String] the code of the station
            position: [position Obj] a position object of the station containing the lat/lon/elev of the station
            channel: [String] the channel of the station (BHZ, HHZ, BDF, etc.)
            name: [String] A string describing the station's name or location (for display purposes only)
            file_name: [String] Location of the .mseed file inside the working directory
        """

        self.network = network
        self.code = code
        self.position = position
        self.channel = channel
        self.file_name = file_name
        self.name = name

def configRead(ini_file):

    class Config:

        def __init__(self):
            pass

    setup = Config()
    config = configparser.ConfigParser()

    if not os.path.exists(ini_file):
        print("FILE ERROR: The input file {:} does not exist!".format(ini_file))
        exit()

    config.read(ini_file)

    ### [General] ###


    try:
        setup.debug = config.get('General', 'debug')
    except:
        setup.debug = 'False'

    setup.debug = (setup.debug.lower() == 'true')

    try:
        setup.fireball_name = config.get('General', 'fireball_name')
        if setup.debug:
            print('DEBUG: fireball_name detected {:}'.format(setup.fireball_name))
    except:
        print("INI ERROR: [General] (fireball_name) Expected fireball_name. Please give the fireball a name, which will be the name of the folder containing its data.")
        exit()

    try:
        setup.run_mode = config.get('General', 'run_mode')
    except:
        #if statements
        print("INI WARNING: [General] (run_mode) Expected run_mode. Choose one of the following: search, replot, manual")

    try:
        setup.difference_filter_all = config.get('General', 'difference_filter_all')
    except:
        print('INI WARNING: [General] (difference_filter_all) Expected difference_filter_all, setting to default value: False')
        setup.difference_filter_all = 'False'

    setup.difference_filter_all = (setup.difference_filter_all.lower() == 'true')

    try:
        setup.get_data = config.get('General', 'get_data')
    except:
        print('INI WARNING: [General] (get_data) Expected get_data, setting to defaul value: False')    
        setup.get_data = 'False'

    setup.get_data = (setup.get_data.lower() == 'true')


    ### [Files] ###

    try:
        setup.working_directory = config.get('Files', 'working_directory')
    except:
        print('INI ERROR: [Files] (working_directory) Expected working_directory. Please give the directory of where a new fireball folder will be stored.')
        exit()

    try:
        setup.arrival_times_file = config.get('Files', 'arrival_times_file')
        if setup.debug:
            print("DEBUG: arrival_times_file detected: {:}".format(setup.arrival_times_file))
    except:
        setup.arrival_times_file = ''
        if setup.debug:
            print("DEBUG: arrival_times_file not detected.")

    try:
        setup.sounding_file = config.get('Files', 'sounding_file')
    except:
        setup.sounding_file = ''
        if setup.debug:
            print("DEBUG: sounding_file not detected.")

    try:
        setup.perturbation_spread_file = config.get('Files', 'perturbation_spread_file')
    except:
        setup.perturbation_spread_file = ''
        if setup.debug:
            print("DEBUG: perturbation_spread_file not detected.")

    try:
        setup.station_picks_file = config.get('Files', 'station_picks_file')
    except:
        setup.station_picks_file = ''
        if setup.debug:
            print("DEBUG: station_picks_file not detected.")

    try:
        setup.replot_points_file = config.get('Files', 'replot_points_file')
    except:
        setup.replot_points_file = ''
        if setup.debug:
            print("DEBUG: replot_points_file not detected.")

    ### [Parameters] ###
    
    try:
        setup.lat_centre = float(config.get('Parameters', 'lat_centre'))
    except:
        setup.lat_centre = ''
        if setup.debug:
            print("DEBUG: lat_centre not detected.")

    try:
        setup.lon_centre = float(config.get('Parameters', 'lon_centre'))
    except:
        setup.lon_centre = ''
        if setup.debug:
            print("DEBUG: lon_centre not detected.")

    try:
        setup.deg_radius = float(config.get('Parameters', 'deg_radius'))
    except:
        setup.deg_radius = ''
        if setup.debug:
            print("DEBUG: deg_radius not detected.")

    try:
        setup.start_datetime = datetime.datetime.strptime(config.get('Parameters', 'start_datetime'), "%Y-%m-%d %H:%M:%S.%f")
    except:
        setup.start_datetime = ''
        if setup.debug:
            print("DEBUG: start_datetime not detected.")

    try:
        setup.end_datetime = datetime.datetime.strptime(config.get('Parameters', 'end_datetime'), "%Y-%m-%d %H:%M:%S.%f")
    except:
        setup.end_datetime = ''
        if setup.debug:
            print("DEBUG: end_datetime not detected.")

    try:
        setup.v_sound = float(config.get('Parameters', 'v_sound'))
    except:
        setup.v_sound = ''
        if setup.debug:
            print("DEBUG: v_sound not detected.")

    ### [Ballistic] ###

    try:
        setup.t0 = float(config.get('Ballistic', 't0'))
    except:
        setup.t0 = ''
        if setup.debug:
            print("DEBUG: t0 not detected.")

    try:
        setup.v = float(config.get('Ballistic', 'v'))
    except:
        setup.v = ''
        if setup.debug:
            print("DEBUG: v0 not detected.")

    try:
        setup.azim = float(config.get('Ballistic', 'azim'))
    except:
        setup.azim = ''
        if setup.debug:
            print("DEBUG: azim not detected.")

    try:
        setup.zangle = float(config.get('Ballistic', 'zangle'))
    except:
        setup.zangle = ''
        if setup.debug:
            print("DEBUG: zangle not detected.")

    try:
        setup.lat_i = float(config.get('Ballistic', 'lat_i'))
    except:
        setup.lat_i = ''
        if setup.debug:
            print("DEBUG: lat_i not detected.")

    try:
        setup.lon_i = float(config.get('Ballistic', 'lon_i'))
    except:
        setup.lon_i = ''
        if setup.debug:
            print("DEBUG: lon_i not detected.")

    try:
        setup.elev_i = float(config.get('Ballistic', 'elev_i'))
    except:
        setup.elev_i = ''
        if setup.debug:
            print("DEBUG: elev_i not detected.")

    try:
        setup.lat_f = float(config.get('Ballistic', 'lat_f'))
    except:
        setup.lat_f = ''
        if setup.debug:
            print("DEBUG: lat_f not detected.")

    try:
        setup.lon_f = float(config.get('Ballistic', 'lon_f'))
    except:
        setup.lon_f = ''
        if setup.debug:
            print("DEBUG: lon_f not detected.")

    try:
        setup.elev_f = float(config.get('Ballistic', 'elev_f'))
    except:
        setup.elev_f = ''
        if setup.debug:
            print("DEBUG: elev_f not detected.")

    try:
        setup.show_ballistic_waveform = config.get('Ballistic', 'show_ballistic_waveform')
    except:
        setup.show_ballistic_waveform = 'false'

    setup.show_ballistic_waveform = (setup.show_ballistic_waveform.lower() == 'true')

    ### [Fragmentation] ###

    try:
        setup.show_fragmentation_waveform = config.get('Fragmentation','show_fragmentation_waveform')
    except:
        setup.show_fragmentation_waveform = 'false'

    try:
        l =                  config.get('Fragmentation', 'fragmentation_point')
        setup.fragmentation_point = eval(l)
    except:
        print("INI WARNING: [Fragmentation] (fragmentation_point) Could not parse fragmentation_point. Needed: [[ lat, lon, elev, time ]]")
        setup.fragmentation_point = []
        setup.show_fragmentation_waveform = 'false'

    setup.show_fragmentation_waveform = (setup.show_fragmentation_waveform.lower() == 'true')

    try:
        l = config.get('Fragmentation', 'manual_fragmentation_search')
        setup.manual_fragmentation_search = eval(l)
    except:
        setup.manual_fragmentation_search = ''

    ### [Restrictions] ### 
   
    try:
        setup.v_fixed = float(config.get('Restrictions', 'v_fixed'))
    except:
        setup.v_fixed = ''
        if setup.debug:
            print("DEBUG: v_fixed not detected.")

    try:
        setup.azimuth_min = float(config.get('Restrictions', 'azimuth_min'))
    except:
        setup.azimuth_min = ''
        if setup.debug:
            print("DEBUG: azimuth_min not detected.")

    try:
        setup.azimuth_max = float(config.get('Restrictions', 'azimuth_max'))
    except:
        setup.azimuth_max = ''
        if setup.debug:
            print("DEBUG: azimuth_max not detected.")

    try:
        setup.zangle_min = float(config.get('Restrictions', 'zangle_min'))
    except:
        setup.zangle_min = ''
        if setup.debug:
            print("DEBUG: zangle_min not detected.")

    try:
        setup.zangle_max = float(config.get('Restrictions', 'zangle_max'))
    except:
        setup.zangle_max = ''
        if setup.debug:
            print("DEBUG: zangle_max not detected.")

    try:
        setup.x_min = float(config.get('Restrictions', 'x_min'))
    except:
        setup.x_min = ''
        if setup.debug:
            print("DEBUG: x_min not detected.")

    try:
        setup.x_max = float(config.get('Restrictions', 'x_max'))
    except:
        setup.x_max = ''
        if setup.debug:
            print("DEBUG: x_max not detected.")

    try:
        setup.y_min = float(config.get('Restrictions', 'y_min'))
    except:
        setup.y_min = ''
        if setup.debug:
            print("DEBUG: y_min not detected.")

    try:
        setup.y_max = float(config.get('Restrictions', 'y_max'))
    except:
        setup.y_max = ''
        if setup.debug:
            print("DEBUG: y_max not detected.")

    try:
        setup.t_min = float(config.get('Restrictions', 't_min'))
    except:
        setup.t_min = ''
        if setup.debug:
            print("DEBUG: t_min not detected.")

    try:
        setup.t_max = float(config.get('Restrictions', 't_max'))
    except:
        setup.t_max = ''
        if setup.debug:
            print("DEBUG: t_max not detected.")

    try:
        setup.v_min = float(config.get('Restrictions', 'v_min'))
    except:
        setup.v_min = ''
        if setup.debug:
            print("DEBUG: v_min not detected.")

    try:
        setup.v_max = float(config.get('Restrictions', 'v_max'))
    except:
        setup.v_max = ''
        if setup.debug:
            print("DEBUG: v_max not detected.")

    try:
        setup.max_error = float(config.get('Restrictions', 'max_error'))
    except:
        setup.max_error = ''
        if setup.debug:
            print("DEBUG: max_error not detected.")

    try:
        setup.restricted_time = datetime.strptime(config.get('Restrictions', 'restricted_time'), "%Y-%m-%d %H:%M:%S.%f")
    except:
        setup.restricted_time = ''
        if setup.debug:
            print("DEBUG: restricted_time not detected.")

    try:
        setup.min_time = float(config.get('Restrictions', 'min_time'))
    except:
        setup.min_time = ''
        if setup.debug:
            print("DEBUG: min_time not detected.")

    try:
        setup.max_time = float(config.get('Restrictions', 'max_time'))
    except:
        setup.max_time = ''
        if setup.debug:
            print("DEBUG: max_time not detected.")

    try:
        setup.traj_tol = float(config.get('Restrictions', 'traj_tol'))
    except:
        setup.traj_tol = ''
        if setup.debug:
            print("DEBUG: traj_tol not detected.")

    try:
        setup.weight_distance_min = float(config.get('Restrictions', 'weight_distance_min'))
    except:
        setup.weight_distance_min = ''
        if setup.debug:
            print("DEBUG: weight_distance_min not detected.")

    try:
        setup.weight_distance_max = float(config.get('Restrictions', 'weight_distance_max'))
    except:
        setup.weight_distance_max = ''
        if setup.debug:
            print("DEBUG: weight_distance_max not detected.")

    try:
        setup.restrict_to_trajectory = config.get('Restrictions', 'restrict_to_trajectory')
    except:
        setup.restrict_to_trajectory = 'false'

    setup.restrict_to_trajectory = (setup.restrict_to_trajectory.lower() == 'true')    

    try:
        setup.search_area = list((config.get('Restrictions', 'search_area')).replace(' ', '').split(','))
    except:
        setup.search_area = ''
        if setup.debug:
            print("DEBUG: search_area not detected.")

    try:
        setup.search_height = list((config.get('Restrictions', 'search_height')).replace(' ', '').split(','))
    except:
        setup.search_height = ''
        if setup.debug:
            print("DEBUG: search_height not detected.")
    
    try:
        setup.search_area = [float(i) for i in setup.search_area]
    except:
        setup.search_area = ''  
        if setup.debug:
            print("DEBUG: search_area not detected.")

    try:
        setup.search_height =   [float(i) for i in setup.search_height]
    except:
        setup.search_height = ''
        if setup.debug:
            print("DEBUG: search_height not detected.")

    ### [Atmosphere] ###

    try:
        setup.enable_winds = config.get('Atmosphere', 'enable_winds')
    except:
        setup.enable_winds = 'false'

    setup.enable_winds = (setup.enable_winds.lower() == 'true')

    try:
        setup.weather_type = config.get('Atmosphere', 'weather_type')
        setup.weather_type = setup.weather_type.lower()
    except:
        setup.weather_type = 'none'

    try:
        setup.grid_size = float(config.get('Atmosphere', 'grid_size'))
    except:
        setup.grid_size = ''
        if setup.debug:
            print("DEBUG: grid_size not detected.")

    ### [Perturbations] ### 

    try:
        setup.perturb = config.get('Perturbations', 'perturb')
    except:
        setup.perturb = 'false'

    setup.perturb = (setup.perturb.lower() == 'true')

    try:
        setup.perturb_method = config.get('Perturbations', 'perturb_method')
        setup.perturb_method = setup.perturb_method.lower()
    except:
        setup.perturb_method = ''
        if setup.debug:
            print("DEBUG: perturb_method not detected.")

    try:
        setup.perturb_times = int(config.get('Perturbations', 'perturb_times'))
    except:
        setup.perturb_times = ''
        if setup.debug:
            print("DEBUG: perturb_times not detected.")

    try:
        setup.observe_frag_no = int(config.get('Perturbations', 'observe_frag_no'))
    except:
        setup.observe_frag_no = ''
        if setup.debug:
            print("DEBUG: observe_frag_no not detected.")

    ### [Speed] ### 

    try:
        setup.fast_ballistic = config.get('Speed', 'fast_ballistic')
    except:
        setup.fast_ballistic = 'true'

    setup.fast_ballistic = (setup.fast_ballistic.lower() == 'true')

    try:
        setup.fit_type = int(config.get('Speed', 'fit_type'))
    except:
        setup.fit_type = ''
        if setup.debug:
            print("DEBUG: fit_type not detected.")

    try:
        setup.n_theta  = int(config.get('Speed', 'n_theta'))
    except:
        setup.n_theta  = ''
        if setup.debug:
            print("DEBUG: n_theta not detected.")

    try:
        setup.n_phi = int(config.get('Speed', 'n_phi'))
    except:
        setup.n_phi = ''
        if setup.debug:
            print("DEBUG: n_phi not detected.")

    try:
        setup.angle_precision = float(config.get('Speed', 'angle_precision'))
    except:
        setup.angle_precision = ''
        if setup.debug:
            print("DEBUG: angle_precision not detected.")

    try:
        setup.angle_error_tol = float(config.get('Speed', 'angle_error_tol'))
    except:
        setup.angle_error_tol = ''
        if setup.debug:
            print("DEBUG: angle_error_tol not detected.")

    ### [PSO] ###

    try:
        setup.maxiter = int(config.get('PSO', 'maxiter'))
    except:
        setup.maxiter = ''
        if setup.debug:
            print("DEBUG: max_error not detected.")

    try:
        setup.swarmsize = int(config.get('PSO', 'swarmsize'))
    except:
        setup.swarmsize = ''
        if setup.debug:
            print("DEBUG: swarmsize not detected.")

    try:
        setup.run_times = int(config.get('PSO', 'run_times'))
    except:
        setup.run_times = ''
        if setup.debug:
            print("DEBUG: run_times not detected.")

    try:
        setup.minstep = float(config.get('PSO', 'minstep'))
    except:
        setup.minstep = ''
        if setup.debug:
            print("DEBUG: minstep not detected.")

    try:
        setup.minfunc = float(config.get('PSO', 'minfunc'))
    except:
        setup.minfunc = ''
        if setup.debug:
            print("DEBUG: minfunc not detected.")

    try:
        setup.phip = float(config.get('PSO', 'phip'))
    except:
        setup.phip = ''
        if setup.debug:
            print("DEBUG: phip not detected.")

    try:
        setup.phig = float(config.get('PSO', 'phig'))
    except:
        setup.phig = ''
        if setup.debug:
            print("DEBUG: phig not detected.")

    try:
        setup.omega = float(config.get('PSO', 'omega'))
    except:
        setup.omega = ''
        if setup.debug:
            print("DEBUG: omega not detected.")

    try:
        setup.pso_debug = config.get('PSO', 'pso_debug')
    except:
        setup.pso_debug = ''
        if setup.debug:
            print("DEBUG: pso_debug not detected.")

    ### [Graphing] ###

    try:
        setup.plot_all_stations = config.get('Graphing', 'plot_all_stations')
    except:
        setup.plot_all_stations = 'true'

    setup.plot_all_stations = (setup.plot_all_stations.lower() == 'true')

    try:
        setup.colortoggle = config.get('Graphing', 'colortoggle')
    except:
        setup.colortoggle = 'true'

    setup.colortoggle = (setup.colortoggle.lower() == 'true')

    try:
        setup.dot_tol = int(config.get('Graphing', 'dot_tol'))
    except:
        setup.dot_tol = ''
        if setup.debug:
            print("DEBUG: dot_tol not detected.")

    try:
        setup.contour_res = int(config.get('Graphing', 'contour_res'))
    except:
        setup.contour_res = ''
        if setup.debug:
            print("DEBUG: contour_res not detected.")

    try:
        l = config.get('Graphing', 'high_f')
        setup.high_f = eval(l)
    except:
        setup.high_f = ''
        if setup.debug:
            print("DEBUG: high_f not detected.")

    try:
        l = config.get('Graphing', 'high_b')
        setup.high_b = eval(l)
    except:
        setup.high_b = ''
        if setup.debug:
            print("DEBUG: high_b not detected.")

    try:
        l = config.get('Graphing', 'rm_stat')
        setup.rm_stat = eval(l)
    except:
        setup.rm_stat = ''
        if setup.debug:
            print("DEBUG: rm_stat not detected.")

    try:
        setup.img_dim = int(config.get('Graphing', 'img_dim'))
    except:
        setup.img_dim = ''
        if setup.debug:
            print("DEBUG: img_dim not detected.")

    try:
        l = config.get('Graphing', 'reported_points')
        setup.reported_points = eval(l)
    except:
        setup.reported_points = ''
        if setup.debug:
            print("DEBUG: reported_points not detected.")

    ### [Extra Stations] ###
    try:
        l = config.get('Extra Stations', 'stations')
        station_list = eval(l)
        setup.stations = []

        for line in station_list:

            pos = position(line[2], line[3], line[4])
            stn = station(line[0], line[1], pos, line[5], line[6], line[7])
            setup.stations.append(stn)

    except:
        setup.stations = ''
        if setup.debug:
            print("DEBUG: stations not detected.")

    return setup

def configParse(setup, ini_type):
    '''
    Parses a config object and outputs a list of variables missing
    '''

    if ini_type == "picks" or ini_type == "get_data":

        print("INI TYPE: 'picks'")

        error_message = """
            CORE REQUIRED VARIABLES
            (The program will not run without these)
            fireball_name
            ini_type
            working_directory
            lat_centre
            lon_centre
            start_datetime
            end_datetime
            """

        if setup.lat_centre == '':
            print(error_message)
            print("Unable to parse lat_centre")
            exit()

        if setup.lon_centre == '':
            print(error_message)
            print("Unable to parse lon_centre")
            exit()

        if setup.start_datetime == '':
            print(error_message)
            print("Unable to parse start_datetime")
            exit()

        if setup.end_datetime == '':
            print(error_message)
            print("Unable to parse end_datetime")
            exit()


    elif ini_type == "supracenter":

        print("INI TYPE: 'supracenter'")

        error_message = """
            CORE REQUIRED VARIABLES
            (The program will not run without these)
            fireball_name
            ini_type
            working_directory
            station_picks_file
            start_datetime
            search_area
            """

        if setup.station_picks_file == '':
            print(error_message)
            print("Unable to parse station_picks_file")
            exit()

        if setup.start_datetime == '':
            print(error_message)
            print("Unable to parse start_datetime")
            exit()

        if setup.search_area == '':
            print(error_message)
            print("Unable to parse search_area")
            exit()

    elif ini_type == 'trajectory':

        print("INI TYPE: 'trajectory'")

        error_message = """
            CORE REQUIRED VARIABLES
            (The program will not run without these)
            fireball_name
            ini_type
            working_directory
            station_picks_file
            start_datetime
            lat_centre
            lon_centre
            """

        if setup.station_picks_file == '':
            print(error_message)
            print("Unable to parse station_picks_file")
            exit()

        if setup.start_datetime == '':
            print(error_message)
            print("Unable to parse start_datetime")
            exit()

        if setup.lat_centre == '':
            print(error_message)
            print("Unable to parse lat_centre")
            exit()

        if setup.lon_centre == '':
            print(error_message)
            print("Unable to parse lon_centre")
            exit()

    else:
        print("INI ERROR: [General] (ini_type) Unable to parse ini_type. Choose either: picks, supracenter, trajectory")
        exit()

    if setup.deg_radius == '':
        print("SETTING TO DEFAULT: Variable undetected, deg_radius -> 2")
        setup.deg_radius = 2

    if setup.v_sound == '':
        print("SETTING TO DEFAULT: Variable undetected, v_sound -> 310")
        setup.v_sound = 310

    #check if enough trajectory data is given
    if (ini_type == 'picks' and setup.show_ballistic_waveform == True) or (ini_type == 'supracenter' and setup.restrict_to_trajectory == True):
        pass
        #Either two points or one point and angles

    if setup.perturb_times == '':
        setup.perturb_times = 0
    setup.perturb_times += 1

    if setup.weight_distance_min == '' or setup.weight_distance_max == '':
        setup.weight_distance_max = 0
        setup.weight_distance_min = 0
        print("SETTING TO DEFAULT: Variable undetected, weight_distance_min -> 0")
        print("SETTING TO DEFAULT: Variable undetected, weight_distance_max -> 0")

    if setup.min_time == '':
        setup.min_time = -3600
        print("SETTING TO DEFAULT: Variable undetected, min_time -> -3600")

    if setup.max_time == '':
        setup.max_time = 3600
        print("SETTING TO DEFAULT: Variable undetected, max_time ->  3600")

    if setup.minfunc == '':
        setup.minfunc = 1e-8
        print("SETTING TO DEFAULT: Variable undetected, minfunc -> 1e-8")

    if setup.minstep == '':
        setup.minstep = 1e-8
        print("SETTING TO DEFAULT: Variable undetected, minstep -> 1e-8")
    
    if setup.phip == '':
        setup.phip = 0.5
        print("SETTING TO DEFAULT: Variable undetected, phip -> 0.5")
    
    if setup.phig == '':
        setup.phig = 0.5
        print("SETTING TO DEFAULT: Variable undetected, phig -> 0.5")
    
    if setup.omega == '':
        setup.omega = 0.5
        print("SETTING TO DEFAULT: Variable undetected, omega -> 0.5")

    if setup.azimuth_min == '':
        setup.azimuth_min = 0
        print("SETTING TO DEFAULT: Variable undetected, azimuth_min -> 0")

    if setup.azimuth_max == '':
        setup.azimuth_max = 360
        print("SETTING TO DEFAULT: Variable undetected, azimuth_max -> 360")

    if setup.zangle_min == '':
        setup.zangle_min = 0
        print("SETTING TO DEFAULT: Variable undetected, zenith_min -> 0")

    if setup.zangle_max == '':
        setup.zangle_max = 90
        print("SETTING TO DEFAULT: Variable undetected, zenith_max -> 90")

    if setup.x_min == '':
        setup.x_min = -200
        print("SETTING TO DEFAULT: Variable undetected, x_min -> -200")

    if setup.x_max == '':
        setup.x_max = 200
        print("SETTING TO DEFAULT: Variable undetected, x_max -> 200")

    if setup.y_min == '':
        setup.y_min = -200
        print("SETTING TO DEFAULT: Variable undetected, y_min -> -200")

    if setup.y_max == '':
        setup.y_max = 200
        print("SETTING TO DEFAULT: Variable undetected, y_max -> 200")

    if setup.t_min == '':
        setup.t_min = -200
        print("SETTING TO DEFAULT: Variable undetected, t_min -> -200")

    if setup.t_max == '':
        setup.t_max = 200
        print("SETTING TO DEFAULT: Variable undetected, t_max -> 200")

    if setup.v_min == '':
        setup.v_min = 11
        print("SETTING TO DEFAULT: Variable undetected, v_min -> 11")

    if setup.v_max == '':
        setup.v_max = 70
        print("SETTING TO DEFAULT: Variable undetected, v_max -> 70")

    if setup.max_error == '':
        setup.max_error = 1000
        print("SETTING TO DEFAULT: Variable undetected, max_error -> 1000")

    if setup.search_area == '':
        setup.search_area = [0, 0, 0, 0]
        print("WARNING: search_area undetected")

    # set up file names
    setup.arrival_times_file = os.path.join(setup.working_directory, setup.fireball_name, setup.arrival_times_file)
    setup.sounding_file = os.path.join(setup.working_directory, setup.fireball_name, setup.sounding_file)
    setup.perturbation_spread_file = os.path.join(setup.working_directory, setup.fireball_name, setup.perturbation_spread_file)
    setup.station_picks_file = os.path.join(setup.working_directory, setup.fireball_name, setup.station_picks_file)
    setup.replot_points_file = os.path.join(setup.working_directory, setup.fireball_name, setup.replot_points_file)

    print("Trajectory Calculation...")
    # find missing trajectory points

    if setup.elev_f == '':
        setup.elev_f = 0

    # Case 0: everything is known
    if setup.lat_i != '' and setup.lon_i != '' and setup.elev_i != '' and setup.lat_f != '' and setup.lon_f != ''  and setup.azim != '' and setup.zangle != '':
        print("Trajectory is known")    

    # Case 1: two points are known
    elif setup.lat_i != '' and setup.lon_i != '' and setup.elev_i != '' and setup.lat_f != '' and setup.lon_f != ''  and setup.azim == '' and setup.zangle == '':
        print("Two points are known")

        a = geo2Loc(setup.lat_f, setup.lon_f, setup.elev_f, setup.lat_i, setup.lon_i, setup.elev_i*1000)
        b = np.array([0, 0, 0])

        d = np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)

        setup.zangle = np.degrees(np.arctan2(d, abs(b[2] - a[2])))

        setup.azim = np.degrees(np.arctan2(b[0] - a[0], b[1] - a[1]))%360

        print("Calculated Trajectory")
        print("Given: Point A: {:8.2f} N {:8.2f} E {:8.2f} m".format(setup.lat_i, setup.lon_i, setup.elev_i))
        print("       Point B: {:8.2f} N {:8.2f} E {:8.2f} m".format(setup.lat_f, setup.lon_f, setup.elev_f))
        print("Calculated: Angles: Zenith {:8.2f} deg , Azimuth {:8.2f} deg".format(setup.azim, setup.zangle))

    # Case 2: top point and angle are known
    elif setup.lat_i != '' and setup.lon_i != '' and setup.elev_i != '' and setup.lat_f == '' and setup.lon_f == ''  and setup.azim != '' and setup.zangle != '':
        print("Top point and angles are known")

        # Tail of the trajectory
        A = geo2Loc(setup.lat_centre, setup.lon_centre, 0, setup.lat_i, setup.lon_i, setup.elev_i*1000)

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
        B = np.array(loc2Geo(setup.lat_centre, setup.lon_centre, 0, B))
        A = np.array(loc2Geo(setup.lat_centre, setup.lon_centre, 0, A))

        print("Created Trajectory between A and B:")
        print("     A = {:10.4f}N {:10.4f}E {:10.2f}m".format(A[0], A[1], A[2]))
        print("     B = {:10.4f}N {:10.4f}E {:10.2f}m".format(B[0], B[1], B[2]))

        setup.lat_f = B[0]
        setup.lon_f = B[1]
        setup.elev_f = B[2]


    # Case 3: bottom point and angle are known
    elif setup.lat_i == '' and setup.lon_i == '' and setup.elev_i == '' and setup.lat_f != '' and setup.lon_f != ''  and setup.azim != '' and setup.zangle != '':
        print("Bottom point and angles are known")

        B = np.array([0, 0, 0])

        # convert angles to radians
        ze = np.radians(setup.zangle)
        az = np.radians(setup.azim)

        # Create trajectory vector
        traj = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

        # backwards propegate the trajectory until it reaches 100000 m up
        n = 100000/traj[2]

        # B is the intersection between the trajectory vector and the ground
        A = n*traj

        # Convert back to geo coordinates
        B = np.array(loc2Geo(setup.lat_centre, setup.lon_centre, 0, B))
        A = np.array(loc2Geo(setup.lat_centre, setup.lon_centre, 0, A))

        print("Created Trajectory between A and B:")
        print("     A = {:10.4f}N {:10.4f}E {:10.2f}m".format(A[0], A[1], A[2]))
        print("     B = {:10.4f}N {:10.4f}E {:10.2f}m".format(B[0], B[1], B[2]))

        setup.lat_i = A[0]
        setup.lon_i = A[1]
        setup.elev_i = A[2]

    # Case 4: nothing is known
    else:
        print("Trajectory unknown")

    #turn coordinates into position objects
    setup.traj_i = position(setup.lat_i, setup.lon_i, setup.elev_i)
    setup.traj_f = position(setup.lat_f, setup.lon_f, setup.elev_f)

def iniBuilder():
    pass

