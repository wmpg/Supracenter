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

from supra.Supracenter.angleConv import angle2NDE, loc2Geo, geo2Loc, strToBool
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
        self.offset = 0

def configWrite(ini_file, setup):
    config = configparser.RawConfigParser()

    config.add_section('General')
    config.set('General', 'fireball_name', str(setup.fireball_name))
    config.set('General', 'difference_filter_all', str(setup.difference_filter_all))
    config.set('General', 'get_data', str(setup.get_data))
    config.set('General', 'run_mode', str(setup.run_mode))
    config.set('General', 'debug', str(setup.debug))

    config.add_section('Files')
    config.set('Files', 'working_directory', str(setup.working_directory))
    config.set('Files', 'arrival_times_file', str(setup.arrival_times_file))
    config.set('Files', 'sounding_file', str(setup.sounding_file))
    config.set('Files', 'perturbation_spread_file', str(setup.perturbation_spread_file))
    config.set('Files', 'station_picks_file', str(setup.station_picks_file))
    config.set('Files', 'replot_points_file', str(setup.replot_points_file))

    config.add_section('Parameters')
    config.set('Parameters', 'lat_centre', str(setup.lat_centre))
    config.set('Parameters', 'lon_centre', str(setup.lon_centre))
    config.set('Parameters', 'deg_radius', str(setup.deg_radius))
    config.set('Parameters', 'start_datetime', str(setup.start_datetime))
    config.set('Parameters', 'end_datetime', str(setup.end_datetime))                    
    config.set('Parameters', 'v_sound', str(setup.v_sound))    

    config.add_section('Ballistic')
    config.set('Ballistic', 't0', str(setup.t0))
    config.set('Ballistic', 'v', str(setup.v))
    config.set('Ballistic', 'azim', str(setup.azim))
    config.set('Ballistic', 'zangle', str(setup.zangle))
    config.set('Ballistic', 'lat_i', str(setup.lat_i))
    config.set('Ballistic', 'lon_i', str(setup.lon_i))
    config.set('Ballistic', 'elev_i', str(setup.elev_i))
    config.set('Ballistic', 'lat_f', str(setup.lat_f))
    config.set('Ballistic', 'lon_f', str(setup.lon_f))
    config.set('Ballistic', 'elev_f', str(setup.elev_f))
    config.set('Ballistic', 'show_ballistic_waveform', str(setup.show_ballistic_waveform))

    config.add_section('Fragmentation')
    config.set('Fragmentation', 'fragmentation_point', str(setup.fragmentation_point))
    config.set('Fragmentation', 'show_fragmentation_waveform', str(setup.show_fragmentation_waveform))
    config.set('Fragmentation', 'manual_fragmentation_search', str(setup.manual_fragmentation_search))

    config.add_section('Restrictions')
    config.set('Restrictions', 'v_fixed', str(setup.v_fixed))
    config.set('Restrictions', 'azimuth_min', str(setup.azimuth_min))
    config.set('Restrictions', 'azimuth_max', str(setup.azimuth_max))
    config.set('Restrictions', 'zangle_min', str(setup.zangle_min))
    config.set('Restrictions', 'zangle_max', str(setup.zangle_max))
    config.set('Restrictions', 'x_min', str(setup.x_min))
    config.set('Restrictions', 'x_max', str(setup.x_max))
    config.set('Restrictions', 'y_min', str(setup.y_min))
    config.set('Restrictions', 'y_max', str(setup.y_max))
    config.set('Restrictions', 't_min', str(setup.t_min))
    config.set('Restrictions', 't_max', str(setup.t_max))
    config.set('Restrictions', 'v_min', str(setup.v_min))
    config.set('Restrictions', 'v_max', str(setup.v_max))
    config.set('Restrictions', 'max_error', str(setup.max_error))
    config.set('Restrictions', 'enable_restricted_time', str(setup.enable_restricted_time))
    config.set('Restrictions', 'restricted_time', str(setup.restricted_time))  
    config.set('Restrictions', 'min_time', str(setup.min_time))
    config.set('Restrictions', 'max_time', str(setup.max_time))
    config.set('Restrictions', 'traj_tol', str(setup.traj_tol))
    config.set('Restrictions', 'restrict_to_trajectory', str(setup.restrict_to_trajectory))
    config.set('Restrictions', 'weight_distance_min', str(setup.weight_distance_min))
    config.set('Restrictions', 'weight_distance_max', str(setup.weight_distance_max))
    config.set('Restrictions', 'search_area', str(setup.search_area))
    config.set('Restrictions', 'search_height', str(setup.search_height))

    config.add_section('Atmosphere')
    config.set('Atmosphere', 'enable_winds', str(setup.enable_winds))
    config.set('Atmosphere', 'weather_type', str(setup.weather_type))
    config.set('Atmosphere', 'grid_size', str(setup.grid_size))

    config.add_section('Perturbations')
    config.set('Perturbations', 'perturb', str(setup.perturb))
    config.set('Perturbations', 'perturb_method', str(setup.perturb_method))
    config.set('Perturbations', 'perturb_times', str(setup.perturb_times))
    config.set('Perturbations', 'observe_frag_no', str(setup.observe_frag_no))

    config.add_section('Speed')
    config.set('Speed', 'fast_ballistic', str(setup.fast_ballistic))
    config.set('Speed', 'fit_type', str(setup.fit_type))
    config.set('Speed', 'n_theta', str(setup.n_theta))
    config.set('Speed', 'n_phi', str(setup.n_phi))
    config.set('Speed', 'angle_precision', str(setup.angle_precision))
    config.set('Speed', 'angle_error_tol', str(setup.angle_error_tol))

    config.add_section('PSO')
    config.set('PSO', 'maxiter', str(setup.maxiter))
    config.set('PSO', 'swarmsize', str(setup.swarmsize))
    config.set('PSO', 'run_times', str(setup.run_times))
    config.set('PSO', 'minfunc', str(setup.minfunc))
    config.set('PSO', 'minstep', str(setup.minstep))
    config.set('PSO', 'phip', str(setup.phip))
    config.set('PSO', 'phig', str(setup.phig))
    config.set('PSO', 'omega', str(setup.omega))
    config.set('PSO', 'pso_debug', str(setup.pso_debug))

    config.add_section('Graphing')
    config.set('Graphing', 'plot_all_stations', str(setup.plot_all_stations))
    config.set('Graphing', 'colortoggle', str(setup.colortoggle))
    config.set('Graphing', 'dot_tol', str(setup.dot_tol))
    config.set('Graphing', 'contour_res', str(setup.contour_res))
    config.set('Graphing', 'high_f', str(setup.high_f))
    config.set('Graphing', 'high_b', str(setup.high_b))
    config.set('Graphing', 'rm_stat', str(setup.rm_stat))
    config.set('Graphing', 'img_dim', str(setup.img_dim))
    config.set('Graphing', 'reported_points', str(setup.reported_points))
    
    config.add_section('Extra Stations')
    config.set('Extra Stations', 'stations', str(setup.stations))

    with open(ini_file, 'w') as configfile:
        config.write(configfile)

def tryFloat(statement):
    
    try:
        return float(statement)
    except:
        return ''

def tryInt(statement):

    try:
        return int(statement)
    except:
        return ''

def tryBool(statement):

    statement = statement.lower()

    return strToBool(statement)

def tryStr(statement):

    return statement

def tryDateTime(statement):

    try:
        try:
            return datetime.datetime.strptime(statement, "%Y-%m-%d %H:%M:%S.%f")
        except:
            return datetime.datetime.strptime(statement, "%Y-%m-%d %H:%M:%S")
    except:
        return ''

def tryEval(statement, try_float=False):
    
    try:
        statement = eval(statement)
    except:
        statement = ''

    if try_float:
        statement = [tryFloat(i) for i in statement]

    return statement

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

    setup.debug = tryBool(config.get('General', 'debug'))
    setup.fireball_name = tryStr(config.get('General', 'fireball_name'))
    setup.run_mode = tryStr(config.get('General', 'run_mode')).lower()
    setup.difference_filter_all = tryBool(config.get('General', 'difference_filter_all'))
    setup.get_data = tryBool(config.get('General', 'get_data'))

    ### [Files] ###

    setup.working_directory = tryStr(config.get('Files', 'working_directory'))
    setup.arrival_times_file = tryStr(config.get('Files', 'arrival_times_file'))
    setup.sounding_file = tryStr(config.get('Files', 'sounding_file'))
    setup.perturbation_spread_file = tryStr(config.get('Files', 'perturbation_spread_file'))
    setup.station_picks_file = tryStr(config.get('Files', 'station_picks_file'))
    setup.replot_points_file = tryStr(config.get('Files', 'replot_points_file'))

    ### [Parameters] ###
    setup.lat_centre = tryFloat(config.get('Parameters', 'lat_centre'))
    setup.lon_centre = tryFloat(config.get('Parameters', 'lon_centre'))
    setup.deg_radius = tryFloat(config.get('Parameters', 'deg_radius'))
    setup.start_datetime = tryDateTime(config.get('Parameters', 'start_datetime'))
    setup.end_datetime = tryDateTime(config.get('Parameters', 'end_datetime'))
    setup.v_sound = tryFloat(config.get('Parameters', 'v_sound'))

    ### [Ballistic] ###
    setup.t0 = tryFloat(config.get('Ballistic', 't0'))
    setup.v = tryFloat(config.get('Ballistic', 'v'))
    setup.azim = tryFloat(config.get('Ballistic', 'azim'))
    setup.zangle = tryFloat(config.get('Ballistic', 'zangle'))
    setup.lat_i = tryFloat(config.get('Ballistic', 'lat_i'))
    setup.lon_i = tryFloat(config.get('Ballistic', 'lon_i'))
    setup.elev_i = tryFloat(config.get('Ballistic', 'elev_i'))
    setup.lat_f = tryFloat(config.get('Ballistic', 'lat_f'))
    setup.lon_f = tryFloat(config.get('Ballistic', 'lon_f'))
    setup.elev_f = tryFloat(config.get('Ballistic', 'elev_f'))

    setup.show_ballistic_waveform = tryBool(config.get('Ballistic', 'show_ballistic_waveform'))

    ### [Fragmentation] ###

    setup.fragmentation_point = tryEval(config.get('Fragmentation', 'fragmentation_point'))
    setup.show_fragmentation_waveform = tryBool(config.get('Fragmentation','show_fragmentation_waveform'))
    setup.manual_fragmentation_search = tryEval(config.get('Fragmentation', 'manual_fragmentation_search'))

    ### [Restrictions] ### 
    setup.v_fixed = tryFloat(config.get('Restrictions', 'v_fixed'))
    setup.azimuth_min = tryFloat(config.get('Restrictions', 'azimuth_min'))
    setup.azimuth_max = tryFloat(config.get('Restrictions', 'azimuth_max'))
    setup.zangle_min = tryFloat(config.get('Restrictions', 'zangle_min'))
    setup.zangle_max = tryFloat(config.get('Restrictions', 'zangle_max'))
    setup.x_min = tryFloat(config.get('Restrictions', 'x_min'))
    setup.x_max = tryFloat(config.get('Restrictions', 'x_max'))
    setup.y_min = tryFloat(config.get('Restrictions', 'y_min'))
    setup.y_max = tryFloat(config.get('Restrictions', 'y_max'))
    setup.t_min = tryFloat(config.get('Restrictions', 't_min'))
    setup.t_max = tryFloat(config.get('Restrictions', 't_max'))
    setup.v_min = tryFloat(config.get('Restrictions', 'v_min'))
    setup.v_max = tryFloat(config.get('Restrictions', 'v_max'))
    setup.max_error = tryFloat(config.get('Restrictions', 'max_error'))

    setup.enable_restricted_time = tryBool(config.get('Restrictions', 'enable_restricted_time'))
    setup.restricted_time = tryDateTime(config.get('Restrictions', 'restricted_time'))

    setup.min_time = tryFloat(config.get('Restrictions', 'min_time'))
    setup.max_time = tryFloat(config.get('Restrictions', 'max_time'))
    setup.traj_tol = tryFloat(config.get('Restrictions', 'traj_tol'))
    setup.weight_distance_min = tryFloat(config.get('Restrictions', 'weight_distance_min'))
    setup.weight_distance_max = tryFloat(config.get('Restrictions', 'weight_distance_max'))

    setup.restrict_to_trajectory = tryBool(config.get('Restrictions', 'restrict_to_trajectory'))

    setup.search_area = tryEval(config.get('Restrictions', 'search_area'), try_float=True)
    setup.search_height = tryEval(config.get('Restrictions', 'search_height'), try_float=True)

    ### [Atmosphere] ###

    setup.enable_winds = tryBool(config.get('Atmosphere', 'enable_winds'))
    setup.weather_type = tryStr(config.get('Atmosphere', 'weather_type')).lower()
    setup.grid_size = tryFloat(config.get('Atmosphere', 'grid_size'))

    ### [Perturbations] ### 

    setup.perturb = tryBool(config.get('Perturbations', 'perturb'))
    setup.perturb_method = tryStr(config.get('Perturbations', 'perturb_method')).lower()
    setup.perturb_times = tryInt(config.get('Perturbations', 'perturb_times'))
    setup.observe_frag_no = tryInt(config.get('Perturbations', 'observe_frag_no'))

    ### [Speed] ### 

    setup.fast_ballistic = tryBool(config.get('Speed', 'fast_ballistic'))
    setup.fit_type = tryInt(config.get('Speed', 'fit_type'))
    setup.n_theta  = tryInt(config.get('Speed', 'n_theta'))
    setup.n_phi = tryInt(config.get('Speed', 'n_phi'))
    setup.angle_precision = tryFloat(config.get('Speed', 'angle_precision'))
    setup.angle_error_tol = tryFloat(config.get('Speed', 'angle_error_tol'))

    ### [PSO] ###
    
    setup.maxiter = tryInt(config.get('PSO', 'maxiter'))
    setup.swarmsize = tryInt(config.get('PSO', 'swarmsize'))
    setup.run_times = tryInt(config.get('PSO', 'run_times'))
    setup.minstep = tryFloat(config.get('PSO', 'minstep'))
    setup.minfunc = tryFloat(config.get('PSO', 'minfunc'))
    setup.phip = tryFloat(config.get('PSO', 'phip'))
    setup.phig = tryFloat(config.get('PSO', 'phig'))
    setup.omega = tryFloat(config.get('PSO', 'omega'))
    setup.pso_debug = tryStr(config.get('PSO', 'pso_debug'))

    ### [Graphing] ###

    setup.plot_all_stations = tryBool(config.get('Graphing', 'plot_all_stations'))
    setup.colortoggle = tryBool(config.get('Graphing', 'colortoggle'))
    setup.dot_tol = tryInt(config.get('Graphing', 'dot_tol'))
    setup.contour_res = tryInt(config.get('Graphing', 'contour_res'))

    setup.high_f = tryEval(config.get('Graphing', 'high_f'))
    setup.high_b = tryEval(config.get('Graphing', 'high_b'))
    setup.rm_stat = tryEval(config.get('Graphing', 'rm_stat'))

    setup.img_dim = tryInt(config.get('Graphing', 'img_dim'))

    setup.reported_points = tryEval(config.get('Graphing', 'reported_points'))

    ### [Extra Stations] ###

    station_list = tryEval(config.get('Extra Stations', 'stations'))
    setup.stations = []

    for line in station_list:

        pos = position(line[2], line[3], line[4])
        stn = station(line[0], line[1], pos, line[5], line[6], line[7])
        setup.stations.append(stn)

    return setup

def setDefaults(setup):

    config_vars =   [setup.perturb_times, setup.weight_distance_min, setup.weight_distance_max, setup.min_time, setup.max_time,
                    setup.minfunc, setup.minstep, setup.phip, setup.phig, setup.omega, setup.azimuth_min, setup.azimuth_max, setup.zangle_min,
                    setup.zangle_max, setup.x_min, setup.x_max, setup.y_min, setup.y_max, setup.t_min, setup.t_max, setup.v_min, setup.v_max,
                    setup.max_error, setup.search_area]

    defaults = [0, 0, 0, -3600, 3600, 1e-8, 1e-8, 0.5, 0.5, 0.5, 0, 360, 0, 90, -200, 200, -200, 200, -200, 200, 11, 70, 1000, [0, 0, 0, 0]]
    
    for ii, element in enumerate(config_vars):
        if element == '':
            element = defaults[ii]

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

    setup = setDefaults(setup)

    setup.perturb_times += 1

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

    try:
        #turn coordinates into position objects
        setup.traj_i = position(setup.lat_i, setup.lon_i, setup.elev_i)
        setup.traj_f = position(setup.lat_f, setup.lon_f, setup.elev_f)
    except:
        setup.traj_i = position(0, 0, 0)
        setup.traj_f = position(0, 0, 0)
        print("Warning: Unable to build trajectory points")

