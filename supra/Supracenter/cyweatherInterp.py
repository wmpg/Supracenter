"""Creates a merged weather profile between two points with a given grid of atmospheric data"""

import numpy as np

import copy
import time
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Supracenter.angleConv import loc2Geo, roundToNearest
from supra.Supracenter.netCDFconv import findECMWFSound, findMERRASound, findUKMOSound
from supra.Supracenter.cyzInteg import zInteg
from supra.Supracenter.bisearch import bisearch

def nearestPoint2D(points_list, point):
    """ HELPER FUNCTION: finds the index of the closest point to point out of a grid of points (points_list)

    Argurments: 
        points_list: [lists] 2 lists containing the x and y values of all the points. It is assumed that these make 
        up a grid.
        point: [list] the x and y values of the point to optimize the distance to.

    Returns:
        i, j: [int] the index of the minimum point (where i is the minimum x, and j is the minimum y)
    """

    # Initialize grid of points
    dists = np.zeros((len(points_list[0]), len(points_list[1])))

    # x-values
    for i in range(len(points_list[0])):

        # y-values
        for j in range(len(points_list[1])):

            # distance    
            d = ((point[0] - points_list[0][i])**2 + (point[1] - points_list[1][j])**2)**0.5
            dists[i, j] = d

    return np.unravel_index(dists.argmin(), dists.shape)


def divLine(a, b, div, points_x,  points_y):
    ''' Divides a 2-D line formed by points a and b into "div" number of divisions, and "rounds" them to the
        closest list of usable points from points_x and points_y

    Arguments:
        a, b: [list, list] start and endpoint of the line to be divided
        div: [float] number of divisions in the line
        points_x, points_y: [list, list] 2D list of points that are able to be used

    Reutrns:
        p: [list]: list of usable points
    '''
    # divisions of each line
    dx = (b[0] - a[0])/(div - 1)
    dy = (b[1] - a[1])/(div - 1)

    # initialize array
    p = np.array([[0.0, 0.0]])

    x_opt = 0
    y_opt = 0 

    for i in range(div):

        # round each division to the usable points
        # x_opt = bisearch(np.flip(points_x, 0), a[0] + i*dx)
        # y_opt = bisearch(np.flip(points_y, 0), b[0] + i*dy) 
        x_opt, y_opt = nearestPoint2D([points_x, points_y], [a[0] + i*dx, b[0] + i*dy])

        # add points to usable points list
        p = np.vstack((p, [points_x[x_opt], points_y[y_opt]]))

    # remove dupes
    p = np.unique(p, axis=0)

    # First row was all zeroes
    p = np.delete(p, 0, 0)

    return p

def interpWeather(s, d, weather_type, dataset):
    """ Returns a merged atmospheric profile from merged weather data between the points s and d
        using the atmospheric data from grid points close to the path.

    Arguments:
        s: [ndarray] [lat, lon, elev] supracenter location
        d: [ndarray] [lat, lon, elev] detector location
        weather_type: [int] type of atmospheric profile 
        dataset: atmospheric profile of the area

    Returns:
        merged_list: [ndarray] merged atmospheric profile of the path between s and d
    """

    # initialize variables
    sounding_list = []
    merged_list = np.array([0, 0, 0, 0])

    # initial points to be used between s and d
    div = 100

    # Try 100 points first, duplicates will be removed
    points = divLine(s, d, div, dataset[0], dataset[1])

    # remove duplicate points, divisions can round to the same point
    # greatly reduces number of points
    divisions = len(points)

    # build for each division
    for i in range(divisions):

        # MERRA
        if weather_type == 'merra':

            # find atmospheric profile for given lat/lon
            loc_sounding = np.array(findMERRASound(points[i][0], points[i][1], dataset))


        # ECMWF
        elif weather_type == 'ecmwf':
            
            # find atmospheric profile for given lat/lon
            loc_sounding = np.array(findECMWFSound(points[i][0], points[i][1], dataset))


        # UKMO
        elif weather_type == 'ukmo':

            # find atmospheric profile for given lat/lon
            loc_sounding = np.array(findUKMOSound(points[i][0], points[i][1], dataset))

        # Fix heights to limit between source and detector
        loc_sounding = zInteg(d[2], s[2], loc_sounding)

        # Add to list of atmospheric profiles
        sounding_list.append(loc_sounding)

    # Linear approximation
    # use atmospheric profiles on a line from source and detector, then splice the profiles along a line in height
    # between the source and detector, then merge together to a new atmospheric profile
    dz = (s[2] - d[2])/divisions

    # Add sections of each profile to the merged profile
    for i in range(divisions):

        # merge lists of each atmospheric profile
        merged_list = np.vstack((merged_list, zInteg(i*dz + d[2], (i + 1)*dz + d[2], np.array(sounding_list[i]))))

    merged_list = np.delete(merged_list, 0, 0)
    
    ### smooth out boundaries
    a = []

    for row in range(1, len(merged_list)):
        
        # Find top of one section and bottom of another (duplicates)
        if merged_list[row, 0] == merged_list[row - 1, 0]:
            
            for i in [1, 2, 3]:

                # interpolate between boundaries
                merged_list[row, i] = (merged_list[row, i] + merged_list[row - 1, i])/2

            a.append(row - 1)

    # remove duplicate rows
    merged_list = np.delete(merged_list, tuple(a), 0)   

    return merged_list, points


def getWeather(S, D, weather_type, ref_pos, dataset, convert=True):
    """ Returns the weather data based off of the weather type:
        weather_type = 'custom' - the weather profile of the given text file
        weather_type = 'merra' - the weather is interpolated between the detector and supracenter for MERRA data.
                        points are taken from the grid which are closest to the line, and the atmospheric data
                        is split up in chunks of even height along the chosen grid points.
        weather_type = 'ecmwf' - ECMWF data '' '' ''
        weather_type = 'ukmo' - UKMO data '' '' '' 

    Arguments:
        S: [lat, lon, elevation] Supracenter location
        D: [lat, lon, elevation] Detector location
        weather_type: [int] see above
        ref_pos: reference position used for converting between geographic and local coordinates
        dataset: atmospheric profile of the area
        convert: [boolean] set to True if using local coordinates, set to False if using lat/lon

    Returns
        sounding: [ndarray] basic atmospheric profile
    OR
        merged_list: [ndarray] interpolated atmospheric profile
        
        points: [ndarray] list of lats and lons used in the atmospheric profile
    """
    # custom weather

    if weather_type == 'custom' or weather_type == 'none':

        # Only return interval between heights
        sounding = zInteg(D[2], S[2], dataset)
        points = []

        return sounding, points
    
    else:

        if convert:
            # convert to unrounded lat/lon
            s = np.array(loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], S))
            d = np.array(loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], D))
            s[2] = S[2]
            d[2] = D[2]
        else:
            s = S
            d = D

        # s and d are in lat/lon, S and D are in local
        # interp weather from grid between s and d
        merged_list, points = interpWeather(s, d, weather_type, dataset)    

        return merged_list, points

    