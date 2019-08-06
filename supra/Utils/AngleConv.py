"""Various coordinate conversions and utility functions for use with Supracenter"""

import numpy as np
from scipy.special import erfc

from wmpl.Utils.Math import rotateVector
from wmpl.Utils.TrajConversions import latLonAlt2ECEF, ecef2ENU, enu2ECEF, ecef2LatLonAlt

def invertColor(color):
    color[0] = 255 - color[0]
    color[1] = 255 - color[1]
    color[2] = 255 - color[2]
    return color

def angle2NDE(angle):
    """ Converts an angle in degrees from an East due North coordinate system to a North due East coordinate system
        Calling this function again will convert it back, since: x = angle2NDE(angle2NDE(x))

    Arguments:
        angle: [float] angle to be converted in degrees.

    Returns:
        angle: [float] converted angle in degrees.
    """

    # Rotate coordinate system 90 degrees CCW
    angle = (angle - 90)%360

    # Flip coordinate system horizontally
    angle = (360 - angle)%360

    return angle


def geo2Loc(lat0, lon0, elev0, lat, lon, elev):
    """ Converts geographic coordinates to local coordinates, while keeping the elevation the same

    Arguments:
        lat0, lon0, elev0: [float] reference coordinates in latitude, longitude, and elevation
        lat, lon, elev: [float] coordinates in geographic system to be converted into local

    Returns:
        local_coord[0], local_coord[1], elev: [float] converted coordinates in local coordinates 
    """

    # Convert to radians
    lat0, lon0 = np.radians(lat0), np.radians(lon0)
    lat, lon = np.radians(lat), np.radians(lon)

    # Convert to local
        # Calculate the ECEF coordinates of the reference position
    x0, y0, z0 = latLonAlt2ECEF(lat0, lon0, elev0)
    ref_ecef = np.array([x0, y0, z0])

    # Convert the geo coordinates of the station into ECEF coordinates
    coord_ecef = latLonAlt2ECEF(lat, lon, elev)

    ### Convert the ECEF coordinates into to local coordinate system ###
    
    local_coord = coord_ecef - ref_ecef

    # Rotate the coordinates so the origin point is tangent to the Earth's surface
    local_coord = np.array(ecef2ENU(lat0, lon0, *local_coord))


    # Ignore height component transformation
    return local_coord[0], local_coord[1], elev


def loc2Geo(lat0, lon0, elev0, local_coord):
    """ Converts local coordinates to geographic, while keeping the elevation the same

    Arguments:
        lat0, lon0, elev0: [float] reference coordinates in latitude, longitude and elevation
        local_coord: [ndarray] coordinates in local system to be converted to geographic system

    Returns:
        lat, lon, local_coord[2]: [float] converted coordinates in geographic coordinates

    """

    # Convert to radians
    lat0, lon0 = np.radians(lat0), np.radians(lon0)

    # Calculate the ECEF coordinates of the reference position
    x0, y0, z0 = latLonAlt2ECEF(lat0, lon0, elev0)
    ref_ecef = np.array([x0, y0, z0])

    # Convert the coordinates back to ECEF
    coord_ecef = np.array(enu2ECEF(lat0, lon0, *local_coord)) + ref_ecef

    # Convert ECEF coordinates back to geo coordinates
    lat, lon, elev = ecef2LatLonAlt(*coord_ecef)

    # Convert to degrees
    lat, lon = np.degrees(lat), np.degrees(lon)

    # Ignore height component transformation
    return lat, lon, local_coord[2]


def roundToNearest(number, nearest):
    """ Rounds a decimal number to the closest value, nearest, given 

    Arguments:
        number: [float] the number to be rounded
        nearest: [float] the number to be rouned to

    Returns:
        rounded: [float] the rounded number
    """

    A = 1/nearest

    rounded = round(number*A)/A

    return rounded


def geopot2Geomet(h):
    """ Converts geopotential heights in an array to geometrical heights

    Arguments:
        h: [ndarray] matrix of geopotential heights

    Returns:
        h: [ndarray] matrix of geometric heights
    """

    # polar radius
    EARTH_RAD = 6356000

    for i in range(len(h)):

        # Conversion factor: http://www.pdas.com/geopot.pdf
        h[i] = EARTH_RAD*h[i]/(EARTH_RAD - h[i])

    return h

def trajRestriction(setup):
    """ create a trajectory defined by two points
    """
    # a is the initial point given
    a = geo2Loc(setup.ref_pos[0], setup.ref_pos[1], setup.ref_pos[2], setup.lat_i, setup.lon_i, setup.elev_i)
    a = np.array(a)
    a[2] *= 1000

    # if setup.traj == '1p':
    #     ze = np.radians(setup.ze)
    #     az = np.radians(setup.az)

    #     b = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])
    # elif setup.traj == '2p':
    b = geo2Loc(setup.ref_pos[0], setup.ref_pos[1], setup.ref_pos[2], setup.lat_f, setup.lon_f, setup.elev_f)
    b = np.array(b)
    b[2] *= 1000
    # else:
    #     print("ERROR: setup.traj type unknown, use 1p or 2p!")

    return a, b

def point2LineDist2D(x, y, a, b):
    """ distance from point (x, y) to a line defined by a and b
    """
    ax, ay = a[0], a[1]
    bx, by = b[0], b[1]

    d = abs((by - ay)*x - (bx - ax)*y + bx*ay - by*ax)/ ((by - ay)**2 + (bx - ax)**2)**0.5

    return d

def point2LineDist3D(a, b, c):
    """ distance from point c to line defined by a and b
    """
    B = np.cross((c-a), (c-b))
    A = b - a

    A = np.sqrt(A.dot(A))
    B = np.sqrt(B.dot(B))

    d = B/A

    return d



def chauvenet(data):
    try:
        indexes = ~np.isnan(data)
        data = data[indexes]
    except:
        pass
    N = len(data)
    std = np.std(data)
    mean = np.mean(data)

    remove = []
    a = []
    
    for ii, element in enumerate(data):
        D = abs(element - mean)/std
        chance = erfc(D**0.5)
        if chance*N < 0.5:
            a.append(ii)
            remove.append(element)

    data = np.delete(data, a)

    return data, remove



if __name__ == "__main__":
    a = [1, 2, 3]
    print(chauvenet(a))
