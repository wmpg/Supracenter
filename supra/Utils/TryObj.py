import datetime

from supra.Utils.Classes import Angle, Position, Trajectory, Supracenter

def strToBool(my_str):

    return (my_str.lower() == 'true')

def tryFloat(statement):

    try:
        return float(statement)
    except:
        return None

def tryInt(statement):

    try:
        return int(statement)
    except:
        return None

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
        return None

def tryEval(statement, try_float=False):
    
    try:
        statement = eval(statement)
    except:
        statement = None

    if try_float:
        statement = [tryFloat(i) for i in statement]

    return statement

def tryPosition(lat, lon, elev):

    try:
        result = Position(lat, lon, elev)
    except:
        result = None

    return result

def tryAngle(angle):
    try:
        if len(angle) == 0:
            return None
    except:
        return None
    try:
        result = Angle(tryFloat(angle))
    except:
        result = None

    return result

def tryTrajectory(t, v, az, ze, i, f):

    # try:
    result = Trajectory(t, v, azimuth=az, zenith=ze, pos_i=i, pos_f=f)
    # except:
    #     result = None

    return result

def trySupracenter(statement, *t):

    supra_list = []
    if len(t) == 0:
        try:
            for event in statement:
                supra_list.append(Supracenter(Position(event[0], event[1], event[2]), event[3]))
        except:
            supra_list = [None]
    else:
        try:
            supra_list.append(Supracenter(statement, t[0]))
        except:
            supra_list = [None]

    return supra_list

def saveDefaults(setup):

        if setup.fireball_name == None:                     setup.fireball_name = 'Untitled Fireball'
        if setup.get_data == None:                          setup.get_data = 'false'
        if setup.run_mode == None:                          setup.run_mode = 'search'
        if setup.debug == None:                             setup.debug = 'false'

        if setup.deg_radius == None:                        setup.deg_radius = 2
        if setup.v_sound == None:                           setup.v_sound = 310

        if setup.show_ballistic_waveform == None:           setup.show_ballistic_waveform = 'false'
        if setup.show_fragmentation_waveform == None:       setup.show_fragmentation_waveform = 'false'

        if setup.azimuth_min == None:                       setup.azimuth_min = 0
        if setup.azimuth_max == None:                       setup.azimuth_max = 359.99
        if setup.zenith_min == None:                        setup.zenith_min = 0
        if setup.zenith_max == None:                        setup.zenith_max = 89.99
        if setup.x_min == None:                             setup.x_min = -200000
        if setup.x_max == None:                             setup.x_max = 200000
        if setup.y_min == None:                             setup.y_min = -200000
        if setup.y_max == None:                             setup.y_max = 200000
        if setup.z_min == None:                             setup.z_min = 0
        if setup.z_max == None:                             setup.z_max = 100000
        if setup.t_min == None:                             setup.t_min = -200
        if setup.t_max == None:                             setup.t_max = 200
        if setup.v_min == None:                             setup.v_min = 11000
        if setup.v_max == None:                             setup.v_max = 30000 
        if setup.max_error == None:                         setup.max_error = 1000000
        if setup.enable_restricted_time == None:            setup.enable_restricted_time = 'false'
        if setup.weight_distance_min == None:               setup.weight_distance_min = 0
        if setup.weight_distance_max == None:               setup.weight_distance_max = 0

        if setup.enable_winds == None:                      setup.enable_winds = 'false'
        if setup.weather_type == None:                      setup.weather_type = 'none'

        if setup.perturb == None:                           setup.perturb = 'false'
        if setup.perturb_method == None:                    setup.perturb_method = 'none'
        if setup.perturb_times == None:                     setup.perturb_times = 0
        if setup.observe_frag_no == None:                   setup.observe_frag_no = 0

        if setup.n_theta == None:                           setup.n_theta = 45
        if setup.n_phi == None:                             setup.n_phi = 90
        if setup.h_tol == None:                             setup.h_tol = 1e-5
        if setup.v_tol == None:                             setup.v_tol = 1000

        if setup.maxiter == None:                           setup.maxiter = 100
        if setup.swarmsize == None:                         setup.swarmsize = 100
        if setup.run_times == None:                         setup.run_times = 1
        if setup.phip == None:                              setup.phip = 0.5
        if setup.phig == None:                              setup.phig = 0.5
        if setup.omega == None:                             setup.omega = 0.5
        if setup.pso_debug == None:                         setup.pso_debug = 'false'
        if setup.minfunc == None:                           setup.minfunc = 1e-8
        if setup.minstep == None:                           setup.minstep = 1e-8

        if setup.contour_res == None:                       setup.contour_res = 10

        return setup