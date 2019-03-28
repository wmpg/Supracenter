'''Reads input config files'''

try:
    # Python 2  
    import ConfigParser as configparser

except:
    # Python 3
    import configparser


import os
import sys

from datetime import datetime 



def configRead(input_name):
    ''' Reads .ini files and runs the Supracenter program with the options inside. Information on each of the
        parameters can be found in the user manual: /docs/User Manual.txt

    Arguments:
        input_name: [string] name of the .ini file

    Returns:
        setup: [Object] object containing all of the options for the program to run
    '''

    class Config:

        def __init__(self):
            pass

    setup = Config()
    config = configparser.ConfigParser()

    if not os.path.exists(input_name):
        print('The input file: ', input_name, 'does not exist!')
        sys.exit()

    config.read(input_name)

    config.sections()

    # General
    try:
        #setup.username =        config.get('General', 'username')
        #setup.password =        config.get('General', 'password')
        setup.search_type =     config.get('General', 'search_type')
        setup.output_name =     config.get('General', 'output_name')
        setup.dir_name =        config.get('General', 'dir_name')
        setup.debug =           config.get('General', 'debug')
    except:
        print("INI ERROR: [General] Error in variables. Required variables: search_type, output_name, dir_name, debug")
        sys.exit()

    # Accuracy
    try:
        setup.enable_winds =    config.get('Accuracy', 'enable_winds')
    except:
        print("INI ERROR: [General] Error in enable_winds variable.")
        sys.exit()

    try:
        setup.fit_type =        int(config.get('Accuracy', 'fit_type'))
        setup.n_theta =         int(config.get('Accuracy', 'n_theta'))
        setup.n_phi =           int(config.get('Accuracy', 'n_phi'))
        setup.swarm_size =      int(config.get('Accuracy', 'swarm_size'))
        setup.max_iter =        int(config.get('Accuracy', 'max_iter'))
    except:
        print("INI ERROR: [Accuracy] fit_type, n_theta, n_phi, swarm_size, and max_iter must be strictly integers!")
        sys.exit()

    try:
        setup.precision =       float(config.get('Accuracy', 'precision'))
        setup.tol =             float(config.get('Accuracy', 'tol'))
        setup.minfunc =         float(config.get('Accuracy', 'minfunc'))
        setup.minstep =         float(config.get('Accuracy', 'minstep'))
        setup.phip =            float(config.get('Accuracy', 'phip'))
        setup.phig =            float(config.get('Accuracy', 'phig'))
        setup.omega =           float(config.get('Accuracy', 'omega'))
        setup.max_error =       float(config.get('Restriction', 'max_error'))
    except:
        print("INI ERROR: [Accuracy] precision, tol, minfunc, minstep, phip, phig, omega, and max_error must be floats!")
        sys.exit()

    # Restriction
    try:
        setup.rest_time =       datetime.strptime(config.get('Restriction', 'rest_time'), "%Y-%m-%d %H:%M:%S.%f")
    except:
        print("INI ERROR: [Restriction] rest_time must be in the form of: YYYY-mm-dd HH:MM:SS.ffffff")
        sys.exit()

    try:
        setup.min_time  =       float(config.get('Restriction', 'min_time'))
        setup.max_time  =       float(config.get('Restriction', 'max_time'))
    except:
        print("INI ERROR: [Restriction] min_time and max_time must be a float!")
        sys.exit()

    try:
        setup.traj = config.get('Restriction', 'traj')
    except:
        print("INI WARNING: [Restriction] traj is not set, use none, 1p or 2p. Defaulting to none")
        setup.traj == 'none'

    if setup.traj == '1p':
        try:
            setup.lat_i = float(config.get('Restriction', 'lat_i'))
            setup.lon_i = float(config.get('Restriction', 'lon_i'))
            setup.elev_i = float(config.get('Restriction', 'elev_i'))
            setup.az = float(config.get('Restriction', 'az'))
            setup.ze = float(config.get('Restriction', 'ze'))
            setup.traj_tol = float(config.get('Restriction', 'traj_tol'))
            print("Status: Detecting Trajectory Restriction")
        except:
            print("INI ERROR: [Restriction] 1p is not set up correctly, required variables: lat_i, lon_i, elev_i, az, ze, traj_tol")   

    elif setup.traj == '2p':
        try:
            setup.lat_i = float(config.get('Restriction', 'lat_i'))
            setup.lon_i = float(config.get('Restriction', 'lon_i'))
            setup.elev_i = float(config.get('Restriction', 'elev_i'))
            setup.lat_f = float(config.get('Restriction', 'lat_f'))
            setup.lon_f = float(config.get('Restriction', 'lon_f'))
            setup.elev_f = float(config.get('Restriction', 'elev_f'))
            setup.traj_tol = float(config.get('Restriction', 'traj_tol'))
        except:
            print("INI ERROR: [Restriction] 2p is not set up correctly, required variables: lat_i, lon_i, elev_i, lat_f, lon_f, elev_f, traj_tol")   

    elif setup.traj == 'none':
        pass

    else:
        print("INI WARNING: [Restriction] traj is not set correctly, use none, 1p or 2p. Defaulting to none")
        setup.traj == 'none'

    # Atmospheric
    try:
        setup.atm_hour =        int(config.get('Atmospheric', 'atm_hour'))
    except:
        print("INI ERROR: [Atmospheric] atm_hour must be strictly an integer!")
        sys.exit()

    try:
        setup.grid_size =       float(config.get('Atmospheric', 'grid_size'))
        setup.speed_of_sound =  float(config.get('Atmospheric', 'speed_of_sound'))
    except:
        print("INI ERROR: [Atmospheric] grid_size and speed_of_sound must be floats!")
        sys.exit()

    try:
        setup.file_name =       config.get('Atmospheric', 'file_name')
        setup.weather_type =    config.get('Atmospheric', 'weather_type').lower()
        setup.get_data =        config.get('Atmospheric', 'get_data')
    except:
        print("INI ERROR: [Atmospheric] Error in variables. Required variables: file_name, weather_type, get_data.")
        sys.exit()

    try: 
        setup.perturb = config.get('Perturb', 'perturb')
        setup.perturb_picks = config.get('Perturb', 'perturb_picks')
    except:
        print("INI ERROR: [Perturb] Error in variables. Required variables: perturb, perturb_picks.")

    try:
        setup.fragno = int(config.get('Perturb', 'fragno'))
    except:
        print("INI WARNING: [Perturb] fragno (if used) must be an integer, setting to 1")
        setup.fragno = 1

    # Station
    try:
        setup.station_name =    config.get('Station', 'station_name')
    except:
        print("INI ERROR: [Station] Error in station_name variable")
        sys.exit()
    
    try:
        setup.ref_time =        datetime.strptime(config.get('Station', 'ref_time'), "%Y-%m-%d %H:%M:%S.%f")
    except:
        print("INI ERROR: [Station] ref_time must be in the form of: YYYY-mm-dd HH:MM:SS.ffffff")
        sys.exit()

    try:
        setup.dmin =            float(config.get('Station', 'dmin'))
        setup.dmax =            float(config.get('Station', 'dmax'))
    except:
        print("INI ERROR: [Station] dmin and dmax must be floats!")
        sys.exit()

    # Search Area
    try:
        setup.search_area =     list((config.get('Search Area', 'search_area')).replace(' ', '').split(','))
        setup.search_height =   list((config.get('Search Area', 'search_height')).replace(' ', '').split(','))
        setup.single_point =    list((config.get('Search Area', 'single_point')).replace(' ', '').split(','))
    except:
        print("INI ERROR: [Search Area] Error in variables. Required variables: search_area, search_height, single_point")
        sys.exit()
    
    try:
        setup.search_area =     [float(i) for i in setup.search_area]
        setup.search_height =   [float(i) for i in setup.search_height]
        setup.single_point =    [float(i) for i in setup.single_point]
    except:
        print("INI ERROR: [Search Area] components of search_area, search_height, and single_point must be floats!")
        sys.exit()

    # make sure that list components have the correct number of indices
    if len(setup.search_area) != 4:
        print("INI ERROR: [Search Area] search_area must have 4 components!")
        sys.exit()

    if len(setup.search_height) != 2:
        print("INI ERROR: [Search Area] search_height must have 2 components!")
        sys.exit()

    if len(setup.single_point) != 3:
        print("INI ERROR: [Search Area] single_point must have 3 components!")
        sys.exit()

    # turn reported points into a list
    try:
        l = config.get('Search Area', 'reported_points')
        setup.reported_points = eval(l)
    except:
        print("INI ERROR: [Search Area] Error in reported_points variable.")
        sys.exit()

    return setup
    