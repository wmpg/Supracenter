"""Various coordinate conversions and utility functions for use with Supracenter"""

import numpy as np
from scipy.special import erfc

from wmpl.Utils.Math import rotateVector
from wmpl.Utils.TrajConversions import date2JD, jd2Date, raDec2ECI, geo2Cartesian, cartesian2Geo, raDec2AltAz, eci2RaDec, latLonAlt2ECEF, ecef2ENU, enu2ECEF, ecef2LatLonAlt

def latLon2Local(lat0, lon0, elev0, lat, lon, elev):
    """ Convert geographic coordinates into a local coordinate system where the reference coordinates will be
        the origin. The positive direction of the X axis points towards the south, the positive direction of
        the Y axis points towards the east and the positive direction of the Z axis points to the zenith at
        the reference coordinates.

    Arguments:
        lat0: [float] reference latitude +N (radians).
        lon0: [float] reference longtidue +E (radians).
        elev0: [float] reference zangle above sea level (meters).
        lat: [float] Latitude +N (radians).
        lon: [float] Longtidue +E (radians).
        elev: [float] zangle above sea level (meters).


    Return:
        (x, y, z): [3 element ndarray] (x, y, z) local coordinates.
    """

    # Calculate the ECEF coordinates of the reference position
    x0, y0, z0 = latLonAlt2ECEF(lat0, lon0, elev0)
    ref_ecef = np.array([x0, y0, z0])

    # Convert the geo coordinates of the station into ECEF coordinates
    coord_ecef = latLonAlt2ECEF(lat, lon, elev)


    ### Convert the ECEF coordinates into to local coordinate system ###
    
    local_coord = coord_ecef - ref_ecef

    # Rotate the coordinates so the origin point is tangent to the Earth's surface
    local_coord = np.array(ecef2ENU(lat0, lon0, *local_coord))

    # Rotate the coordinate system so X points towards the south and Y towards the east
    local_coord = rotateVector(local_coord, np.array([0, 0, 1]), np.pi/2)

    ######


    return local_coord




def local2LatLon(lat0, lon0, elev0, local_coord):
    """ Convert local coordinates into geographic coordinates. See latLon2Local for more details.

    Arguments:
        lat0: [float] reference latitude +N (radians).
        lon0: [float] reference longtidue +E (radians).
        elev0: [float] reference zangle above sea level (meters).
        local_coord: [3 element ndarray] (x, y, z):
            - x: [float] Local X coordinate (meters).
            - y: [float] Local Y coordinate (meters).
            - z: [float] Local Z coordinate (meters).

    Return:
        (lat, lon, elev): [3 element ndarray] Geographic coordinates, angles in radians, zangle in meters.
    """


    # Calculate the ECEF coordinates of the reference position
    x0, y0, z0 = latLonAlt2ECEF(lat0, lon0, elev0)
    ref_ecef = np.array([x0, y0, z0])


    # Rotate the coordinate system back to ENU
    local_coord = rotateVector(local_coord, np.array([0, 0, 1]), -np.pi/2)

    # Convert the coordinates back to ECEF
    coord_ecef = np.array(enu2ECEF(lat0, lon0, *local_coord)) + ref_ecef

    # Convert ECEF coordinates back to geo coordinates
    lat, lon, elev = ecef2LatLonAlt(*coord_ecef)


    return lat, lon, elev

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
    lat0, lon0 = lat0*np.pi/180, lon0*np.pi/180
    lat, lon = lat*np.pi/180, lon*np.pi/180

    # Convert to local
    x, y, z = latLon2Local(lat0, lon0, elev0, lat, lon, elev)

    # Rotate to a +X = East, +Y = North grid
    local_coord = np.array([x, y, z])
    local_coord = rotateVector(local_coord, np.array([0, 0, 1]), -np.pi/2)

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
    lat0, lon0 = lat0*np.pi/180, lon0*np.pi/180

    # Rotate to +X = South, +Y = East grid
    local_coord = rotateVector(local_coord, np.array([0, 0, 1]), np.pi/2)

    # Convert to geographic
    lat, lon, elev = local2LatLon(lat0, lon0, elev0, local_coord)

    # Convert to degrees
    lat, lon = lat*180/np.pi, lon*180/np.pi

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

def strToBool(my_str):

    return (my_str.lower() == 'true')

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
    pass
    # # convert 90 east due north to north due east
    # print("angle", 90)
    # print("angle2NDE", angle2NDE(90))
    # # calling again will convert it back
    # print("angle2EDN", angle2NDE(angle2NDE(90)))

    # # convert lat/lon/elev 49.7/-83.3/500m from geographic coordinates to local coordinates, using a reference 
    # # coordinate of 50/-83/0m
    # print("Geographic", 49.7, -83.3, 500)
    # print("Local", geo2Loc(50, -83, 0, 49.7, -83.3, 500))
    # print("Geographic", loc2Geo(50, -83, 0, geo2Loc(50, -83, 0, 49.7, -83.3, 500))) 

    # # round 0.78 to the nearest multiple of 0.25
    # print("round_to_nearest", roundToNearest(0.78, 0.25))
