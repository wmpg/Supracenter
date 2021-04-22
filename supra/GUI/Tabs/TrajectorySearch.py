import os

from supra.GUI.Tools.GUITools import *

from supra.Utils.Classes import Position, Trajectory, Angle, Supracenter
from supra.Utils.pso import pso

from supra.Fireballs.SeismicTrajectory import trajSearch, timeOfArrival

def psoTrajectory(station_list, bam, prefs):

    point_on_traj = None
    # point_on_traj = Supracenter(Position(46.31, -121.80, 50000), -40.0)
    # point_on_traj = Supracenter(Position(45.8, -122.00, 20000), -7.0)

    ref_pos = Position(bam.setup.lat_centre, bam.setup.lon_centre, 0)
    # ref_pos = bam.setup.trajectory.pos_f

    if bam.setup.pos_min.isNone() or bam.setup.pos_max.isNone():
        errorMessage('Search boundaries are not defined!', 2, info='Please define the minimum and maximum parameters in the "Sources" tab on the left side of the screen!')
        return None

    bam.setup.pos_min.pos_loc(ref_pos)
    bam.setup.pos_max.pos_loc(ref_pos)

    if point_on_traj is None:

        bounds = [
            (bam.setup.pos_min.x, bam.setup.pos_max.x), # X0
            (bam.setup.pos_min.y, bam.setup.pos_max.y), # Y0
            (bam.setup.t_min, bam.setup.t_max), # t0
            (bam.setup.v_min, bam.setup.v_max), # Velocity (m/s)
            (bam.setup.azimuth_min.deg, bam.setup.azimuth_max.deg),     # Azimuth
            (bam.setup.zenith_min.deg, bam.setup.zenith_max.deg)  # Zenith angle
            ]

    else:

        bounds = [
            (bam.setup.v_min, bam.setup.v_max), # Velocity (m/s)
            (bam.setup.azimuth_min.deg, bam.setup.azimuth_max.deg),     # Azimuth
            (bam.setup.zenith_min.deg, bam.setup.zenith_max.deg)  # Zenith angle
            ]

    lower_bounds = [bound[0] for bound in bounds]
    upper_bounds = [bound[1] for bound in bounds]
    
    if prefs.debug:
        print('Free Search')

    import matplotlib.pyplot as plt

    plt.ion()

    fig, ax = plt.subplots(2, 3, sharey='row')
    ax[0, 0].set_ylabel("Total Error")
    ax[1, 0].set_ylabel("Total Error")
    ax[0, 0].set_xlabel("Latitude")
    ax[0, 1].set_xlabel("Longitude")
    ax[0, 2].set_xlabel("Time")
    ax[1, 0].set_xlabel("Velocity")
    ax[1, 1].set_xlabel("Azimuth")
    ax[1, 2].set_xlabel("Zenith")



    plot = []

    for i in range(2):
        for j in range(3):

            plot.append(ax[i, j].scatter([], []))

    x, fopt = pso(trajSearch, lower_bounds, upper_bounds, args=(station_list, ref_pos, bam, prefs, plot, ax, fig, point_on_traj), \
        maxiter=prefs.pso_max_iter, swarmsize=prefs.pso_swarm_size, \
        phip=prefs.pso_phi_p, phig=prefs.pso_phi_g, debug=False, omega=prefs.pso_omega, \
        particle_output=False)

    # if point_on_traj is None:
    print('Results:')
    print('X: {:.4f}'.format(x[0]))
    print('Y: {:.4f}'.format(x[1]))
    print('Time: {:.4f}'.format(x[2]))
    print('Velocity: {:.4f}'.format(x[3]))
    print('Azimuth: {:.4f}'.format(x[4]))
    print('Zenith: {:.4f}'.format(x[5]))
    print('Adjusted Error: {:.4f}'.format(fopt))
    # else:
    #     print('Results:')
    #     print('Velocity: {:.4f}'.format(x[0]))
    #     print('Azimuth: {:.4f}'.format(x[1]))
    #     print('Zenith: {:.4f}'.format(x[2]))
    #     print('Adjusted Error: {:.4f}'.format(fopt))

    # if point_on_traj is None:
    geo = Position(0, 0, 0)
    geo.x = x[0]
    geo.y = x[1]
    geo.z = 0
    geo.pos_geo(ref_pos)

    print('Geometric Landing Point:')
    print(geo)

    stat_names = []
    stat_pick = []

    final_traj = Trajectory(x[2], x[3], zenith=Angle(x[5]), azimuth=Angle(x[4]), pos_f=geo)
    points = final_traj.findPoints(gridspace=100, min_p=17000, max_p=50000)

    for stn in station_list:
        stat_names.append("{:}-{:}".format(stn[1], stn[2]))

        t_nom, t_pert = timeOfArrival(np.array([stn[3], stn[4], stn[5]]), final_traj, bam, prefs, points, ref_loc=ref_pos)
        
        stat_pick.append(t_nom - stn[6])


    return [x, fopt, geo, stat_names, stat_pick]

    # else:
    #     return [x, fopt]

def getStationList(file_name): 
    """ Reads station .csv file and produces a list containing the station's position and signal time. 
        Accepts files exported from MakeIRISPicks.py. A custom file can be made in the following form:
        *header
        pick_group, station_network, station_code, latitude (deg), longitude (deg), zangle (m), time_of_arrival (Julian Date)

    Arguments:
        file_name: [string] location of the station.csv file

    Returns:
        data: [ndarray] parsed station location and times
    """   
    with open(file_name) as f:

        # Skip the header
        for i in range(1):
            next(f)

        data = []
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '').replace('\t', '')

            # Split the line by the delimiter
            line = line.split(',')

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = (entry.strip())
                if i in [3, 4, 5, 6]:
                    line[i] = float(line[i])
                if i in [3, 4]:
                    line[i] = np.radians(line[i])
            # Add the contents of the line to the data list
            data.append(line)

        return data


def trajectorySearch(bam, prefs):

    stat_file = os.path.join(prefs.workdir, bam.setup.fireball_name, bam.setup.station_picks_file)
    
    try:
        station_list = getStationList(stat_file)
    except TypeError as e:
        errorMessage('Unexpected station list location!', 2, info="Can not find where 'station_picks_file' is!", detail='{:}'.format(e))
        return None
    except FileNotFoundError as e:
        errorMessage('Station Picks File was not found', 2, info="A .csv station picks file is required!", detail='{:}'.format(e))
        return None

    class StationPick:
        def __init__(self):
            pass

    station_obj_list = []
    for i in range(len(station_list)):
        
        stnp = StationPick()

        stnp.group = int(station_list[i][0])
        stnp.network = station_list[i][1]
        stnp.code = station_list[i][2]
        stnp.position = Position(np.degrees(float(station_list[i][3])), \
                                 np.degrees(float(station_list[i][4])), \
                                 float(station_list[i][5]))
        stnp.position.pos_loc(Position(bam.setup.lat_centre, bam.setup.lon_centre, 0))
        stnp.time = float(station_list[i][7])


        station_obj_list.append([stnp.group, stnp.network, stnp.code, \
                                stnp.position.x, stnp.position.y, stnp.position.z, \
                                stnp.time])



    return psoTrajectory(station_obj_list, bam, prefs)