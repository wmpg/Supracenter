import numpy as np
from supra.Utils.Classes import Position, Trajectory
import subprocess


atm = 'atm_file'
output = 'output_file'

def readDarkflight(output, mass, header=30, partial=False):
    """ Reads output from a darkflight output.txt file, and adds a mass "stamp" to each line of data, so that
        each point can be associated with a mass

    Arguments:
        output: [string] location of output file
        mass: [float] mass stamp to be added to the point
        header: [int] number of header lines in the file
        partial: [Boolean] whether to take all data (False), or every 20th line (True).
    """

    with open(output) as f:

        for i in range(header):
            next(f)

        data = np.array([0]*10)
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '')

            # Split the line by the delimiter
            line = line.split()

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = float(entry.strip())
            line.append(mass)

            # Add the contents of the line to the data list
            data = np.vstack((data, np.array(line)))

        # First row was all zeroes
        data = np.delete(data, 0, 0)

        if partial:
            # Take every 20th line
            return data[1::20, :]

        else:

            # Take every line
            return data

def darkFlight(source_loc, m, az, ze):
    darkflight_dir = "/home/luke/source/Supracenter/supra/Fireballs/"

    p = subprocess.Popen(['./darkflight', '--src', source_loc.lat + ',' + source_loc.lon + ',' + 0,\
                '--vel', 4, '--az', az, '--zn', ze, '--mas', m, \
                '--atm', atm, '-out', output], cwd=darkflight_dir)

 
    p.communicate()

def makeTraj(x, args):

    base_points = args

    traj = Trajectory(0, 20000, zenith=Angle(x[1]), azimuth=Angle(x[0]), pos_f=Position(x[2], x[3], 0))
    P = traj.trajInterp(div=20, write=True)

    masses = []

    for pt in base_points:
        masses.append(pt[1])

    error_pts = []

    for source_loc in P:
        error = 0
        for ii, m in enumerate(masses):
            az = (180 + traj.azimuth.deg)%360
            ze = traj.zenith.deg
            darkFlight(source_loc, m, az, ze)
            data = readDarkflight(output, m)
            lat = data[-1, 0]
            lon = data[-1, 1]
            alt = data[-1, 2]

            error +=    (np.sqrt((lat-base_point[ii][0].lat)**2 + \
                                 (lon-base_point[ii][0].lon)**2 + \
                                 (alt-base_point[ii][0].elev)**2))

        error_pts.append(error)

    return np.nanmin(error_pts)

def psoDark():
    azim_min = 330
    azim_max = 335
    zenith_min = 40
    zenith_max = 45
    lat_min = 45
    lat_max = 46
    lon_min = 15
    lon_max = 16

    base_points = [[Position(45.867428, 15.075969, 0), 0.469],
                   [Position(45.816811, 15.112467, 0), 0.204],
                   [Position(45.811897, 15.131436, 0), 0.048]]

    bounds_min = [azim_min, zenith_min, lat_min, lon_min]
    bounds_max = [azim_max, zenith_max, lat_max, lon_max]


    x_opt_temp, f_opt = pso(makeTraj, bounds_min, bounds_max, args=base_points, processes=multiprocessing.cpu_count(), particle_output=False) 

    # pso calculates darkflight along trajectory and finds releases with least error
    # pso finds trajectory with least error

if __name__ == '__main__':
    psoDark()
