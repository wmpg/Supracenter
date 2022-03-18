"""Reads station data from a file"""

import numpy as np
import sys
sys.path.insert(0,'..')

from supra.Utils.AngleConv import geo2Loc
# from wmpl.Utils.TrajConversions import date2JD
# from supra.Fireballs.Program import position

def readStationDat(station_name):
    """ Reads the station data .txt file for use with convStationDat (made from the template)

    Arguments:
        station_name: [String] station data .txt file name
        
    Returns:
        data: [ndarray[float]] station location and arrival times
        names: [list] name of each station in order given
    """

    with open(station_name) as f:

        # Skip the header
        for i in range(2):
            next(f)

        names = []

        data = np.array([0, 0, 0, 0])
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '')

            # Split the line by the delimiter
            line = line.split(',')
            names.append(line[0])
            line = line[1:]

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = float(entry.strip())

            # Add the contents of the line to the data list
            data = np.vstack((data, line))

        # First row was all zeroes
        data = np.delete(data, 0, 0)

        return data, names

def readStationCSV(station_name):
    """ Reads the station data .csv file for use with convStationDat (exported from MakeIRISPicks.py)

    Arguments:
        station_name: [String] station data .csv file name

    Returns:
        data: [ndarray[float]] station location and arrival times
        names: [list] name of each station in order given
    """

    with open(station_name) as f:

        # Skip the header
        for i in range(1):
            next(f)

        names = []

        data = np.array([0, 0, 0, 0, 0, 0])
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '').replace('\t', '')

            # Split the line by the delimiter
            line = line.split(',')
            names.append(str(line[1]).strip() + '-' + str(line[2]).strip())
            # weights.append(float(line[-1]))
            line = line[3:]

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = float(entry.strip())
            # Add the contents of the line to the data list

            data = np.vstack((data, line))

        # First row was all zeroes
        data = np.delete(data, 0, 0)

        return data, names


def dweight(data, d_min, d_max):
    ''' Optional weighting of far away stations with a cosine taper. The stations weights will be adjusted to 
        fit the model. The weight will change to the model if its weight was not set to zero. If 
        dmax is set to zero, all custom weights will be used. Setting all custom weights to 1 will return a
        purely automatic weighting solution, for example.

        See SUPRACENTER pg 21 for more information
        Using ref_pos (average station position) as the 0 km point.
    
    Arguments:
        data: [ndarray] array of station locations [x, y, z] and arrival times
        dmin: [float] stations up to this distance from the first station have a weight of 1
        dmax: [float] stations between dmin and dmax have weights as a cosine taper. Stations further than dmax 
                    have a weight of 0

    Returns:
        w: [ndarray] array of fit weights for each stations

    '''

    ns = len(data)
    w = [1]*ns

    ### OVERRIDE
    return w

    # full weighting
    if d_max == 0:
        return w

    # first station
    fs = np.argmin(data[:, 3])

    # initialize vector of distances to fs
    dists = np.zeros(ns)

    for i in range(ns):
        dists[i] = np.sqrt((data[i, 0] - data[fs, 0])**2 + (data[i, 1] - data[fs, 1])**2)
        
        # after d_max, no weight
        if dists[i] > float(d_max):
            w[i] = 0

        # between d_min and d_max
        elif dists[i] >= d_min and dists[i] <= d_max and w[i] != 0: 
            
            # cosine taper
            w[i] = (np.cos(np.pi/(d_max - d_min)*(dists[i] - (d_max - d_min))) + 1)/2

        # before d_min, full weight
        elif dists[i] < d_min and w[i] != 0:
            w[i] = 1

    return w


def convStationDat(station_name, ref_pos, d_min=0, d_max=100000):
    """ Reads and converts the station data, and returns it in a usable format
    
    Arguments:
        station_name: [str] File name containing the station data 
        dmin: [float] stations up to this distance from the first station have a weight of 1
        dmax: [float] stations between dmin and dmax have weights as a cosine taper. Stations further than dmax 
                    have a weight of 0
        ref_time: [floats] the time which the station arrival times are in reference to

    Returns:
        data: [array] location data of each station [x, y, z] in local coordinates
        names: [lst] string name of each station
        weights: [lst] weight of each station (0 - 1, 1 - default, 0 - ignore)
        ref_pos: [array] average geographic coordinate between stations to use as a reference for the local coordinate system
    """
    
    # Read station file
    if '.csv' in station_name:
        data, names = readStationCSV(station_name)

        for i in range(0, len(data)):
            data[i, 3] = data[i, 4]
    else:
        data, names = readStationDat(station_name)

    # convert station positions to local coordinate system
    for i in range(len(data)):
        data[i, 0], data[i, 1], data[i, 2] = geo2Loc(ref_pos.lat, ref_pos.lon, ref_pos.elev, data[i,0], data[i,1], data[i,2])       

    # weighting scheme
    if (d_min == None) or (d_max == None):
        d_max = 0

    weights = dweight(data, d_min, d_max)

    return data, names, weights
