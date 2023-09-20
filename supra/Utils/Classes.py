# -*- coding: utf-8 -*-
import numpy as np
import copy

from supra.Utils.AngleConv import geo2Loc, loc2Geo, angle2NDE

from wmpl.Utils.TrajConversions import latLonAlt2ECEF, ecef2LatLonAlt, ecef2ENU

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore

try:
    from geopy import distance
except:
    pass

class Config:

    def __init__(self):
        pass


class Constants:

    def __init__(self):

        # Ideal gas constant
        self.R = 8.31432 # J/K mol

        # Kelvins to C
        self.K = 273.15

        # Heat capacity ratio
        self.GAMMA = 1.40

        # Molar mass of the air
        self.M_0 = 0.0289644 #kg/mol

        # standard gravity at mean sea level
        self.g_0 = 9.80665 #m/s^2

        self.ABSORPTION_COEFF = 2e-4 # kg/m^2

        self.SCALE_HEIGHT = 8500 # meters

        self.R_EARTH = 6378100

        self.eps = 7./3 - 4./3 -1

class Pick:
    def __init__(self, t, stn, stn_no, channel, group):
        self.time = t
        self.stn = stn
        self.stn_no = stn_no
        self.channel = channel
        self.group = group


class Angle:

    def __init__(self, deg):
        
        if deg != None:
            self.deg = deg
            self.rad = np.radians(deg)
        else:
            self.deg = None
            self.rad = None

    def __str__(self):
        return '{:6.2f} deg'.format(self.deg) 

    def __add__(self, other):
        return (self.deg + other.deg)

    def __sub__(self, other):
        return (self.deg - other.deg)

    def __gt__(self, other):
        return self.deg > other.deg

    def invert(self):
        # Rotate coordinate system 90 degrees CCW
        angle = (self.deg - 90)%360

        # Flip coordinate system horizontally
        angle = (360 - angle)%360

        self.deg = angle

class Position:
    """ A latitude/longitude/elevation object describing a position in 3D space
    Give lat [-90, 90], lon [-180, 180] in degrees
    Give elev in m 
    """

    def __init__(self, lat, lon, elev):

        self.lat  = lat
        self.lon  = lon
        self.elev = elev

        try:
            self.lat_r = np.radians(lat)
            self.lon_r = np.radians(lon)
        except (AttributeError, TypeError):
            self.lat_r = None
            self.lon_r = None


    def __str__(self):
        degree_sign= 'o'
        try:
            result = "Lat: {:8.4f}{:}N, Lon: {:8.4f}{:}E, Elev: {:10.2f} m".format(self.lat, degree_sign, self.lon, degree_sign, self.elev)
        except TypeError:
            result = "Position is None Type"
        return result

    def __sub__(self, other):

        self.pos_loc(other)
        other.pos_loc(other)

        self.x -= other.x
        self.y -= other.y
        self.z -= other.z

        return self

    def __mul__(self, other):
        
        self.x *= other
        self.y *= other
        self.z *= other


        return self


    def isNone(self):

        if self.lat == None or self.lon == None or self.elev == None:
            return True

        return False
    
    def pos_loc(self, ref_pos):
        """
        Converts a position object to local coordinates based off of a reference position

        Arguments:
            ref_pos: [position Obj] the position of a reference location 
            used to convert another position object to local coordinates
        """

        self.x, self.y, self.z = geo2Loc(ref_pos.lat, ref_pos.lon, ref_pos.elev, self.lat, self.lon, self.elev)
        self.xyz = np.array([self.x, self.y, self.z])

    def pos_geo(self, ref_pos):
        """ Converts a position object to geometric coordinates
        """

        self.lat, self.lon, self.elev = loc2Geo(ref_pos.lat, ref_pos.lon, ref_pos.elev, [self.x, self.y, self.z])
        self.xyz = np.array([self.x, self.y, self.z])

    def pos_distance(self, other):
        """ 3D distance between positions 'self' and 'other'
        """

        self.pos_loc(self)
        other.pos_loc(self)

        A = self.xyz
        B = other.xyz

        return np.sqrt((A[0] - B[0])**2 + (A[1] - B[1])**2 + (A[2] - B[2])**2)

    def pos_distance_components(self, other):
        """ returns dx dy dz between two points from A to B
        """
        self.pos_loc(self)
        other.pos_loc(self)

        A = self.xyz
        B = other.xyz

        dx, dy, dz = B[0] - A[0], B[1] - A[1], B[2] - A[2]

        return dx, dy, dz

    def ground_distance(self, other, geospy=False):
        """ 2D horizontal distance between positions 'self' and 'other'
        """

        self.pos_loc(self)
        other.pos_loc(self)

        if geospy:
            return distance.distance((self.lat, self.lon), (other.lat, other.lon)).km*1000 

        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)

    def ground_latlon_distance(self, other):
        """ longer distances between two points
        """
        consts = Constants()

        a = np.sin((self.lat_r - other.lat_r)/2)**2 + np.cos(self.lat_r)*np.cos(other.lat_r)*np.sin((self.lon_r - other.lon_r)/2)**2
        c = 2*np.arctan2(np.sqrt(a), np.sqrt(1 - a))
        d = consts.R_EARTH*c

        return d

    def line_between_2d(self, other, show_error=False):
        """ Generates a 2D line between two positions 
        """

        m = (other.y - self.y)/(other.x - self.x) 

        b_1 = other.y - m*other.x
        b_2 = self.y - m*self.x

        if show_error:
            print("Intercept Error: {:}".format(b_1 - b_2))

        print("y = mx + b")
        print("m = {:}".format(m))
        print("b = {:}".format((b_1 + b_2)/2))

    def mag(self):

        return (self.x**2 + self.y**2 + self.z**2)**0.5

    def angleBetween(self, other):

        dLon = other.lon_r - self.lon_r

        y = np.sin(dLon)*np.cos(other.lat_r)
        x = np.cos(self.lat_r)*np.sin(other.lat_r) - np.sin(self.lat_r)*np.cos(other.lat_r)*np.cos(dLon)
        brng = np.degrees(np.arctan2(y, x))
        # brng = angle2NDE(np.degrees(brng))

        return brng


class Vector3D:

    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.xyz = np.array([x, y, z])

    def __add__(self, other):
        result = Position(0, 0, 0)
        result.x = self.x + other.x
        result.y = self.y + other.y
        result.z = self.z + other.z
        return result

    def __mul__(self, other):
        result = Vector3D(None, None, None)
        result.x = self.x * other
        result.y = self.y * other
        result.z = self.z * other
        return result

    def __str__(self):
        return '[ {:.4f}, {:.4f}, {:.4f}]'.format(self.x, self.y, self.z)

    def mag(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)


class Annote:

    def __init__(self, title, time, length, group, source, height, notes, color, bandpass):

        self.title = title
        self.time = time
        self.length = length
        self.group = group
        self.source = source
        self.height = height
        self.notes = notes
        self.color = color
        self.bandpass = bandpass


    def __str__(self):

        return 'Annotation Object "{:}" located at {:.2f} s'.format(self.title, self.time)



class Supracenter:
    def __init__(self, position, t):
        self.position = position
        self.time = t

    def __str__(self):
        try:
            A = "Position: Lat: {:8.4f} deg N, Lon: {:8.4f} deg E, Elev: {:10.2f} m | Time: {:8.3f} s".format(self.position.lat, self.position.lon, self.position.elev, self.time)
        except:
            A = "Undefined Supracenter Position"
        return A

    def toList(self):
        return [self.position.lat, self.position.lon, self.position.elev, self.time]       

class Trajectory:
    def __init__(self, t, v, zenith=None, azimuth=None, pos_i=Position(None, None, None), pos_f=Position(None, None, None), v_f=None, verbose=False):
        """ A trajectory object

        Arguements:
        t [float] - time after reference the trajectory intersects with the ground (if continued at velocity v) [s]
        v [float] - average velocity of the meteor [m/s]
        zenith [Angle Obj] - zenith angle (angle between vertical and meteor if meteor was coming towards observer)
        azimuth [Angle Obj] - azimuth angle (angle East from North the trajectory vector points)
        pos_i [Position Obj] - 3D starting position of meteor
        pos_f [Position Obj] - 3D final position of meteor
        v_f -- old variable
        """

        TOL = 2
        self.t = t

        # Find vector from angles
        try:

            self.vector = Vector3D(np.sin(azimuth.rad)*np.sin(zenith.rad), \
                                   np.cos(azimuth.rad)*np.sin(zenith.rad), \
                                                      -np.cos(zenith.rad))
        except:
            self.vector = None

        # Case 1: everything is given
        if pos_i.lat is not None and pos_f.lat is not None and zenith is not None and azimuth is not None:

            if verbose:
                print("Case 1: Everything is given")
            pos_i.pos_loc(pos_f)
            pos_f.pos_loc(pos_f)

            scale = (pos_f.z - pos_i.z) / self.vector.z

        # Case 2: end points are given
        elif pos_i.lat is not None and pos_f.lat is not None:


            if verbose:
                print("Case 2: End points are given")

            pos_i.pos_loc(pos_f)
            pos_f.pos_loc(pos_f)

            dx = pos_f.x - pos_i.x
            dy = pos_f.y - pos_i.y
            dz = pos_f.z - pos_i.z
            dh = np.sqrt(dx**2 + dy**2)


            # A = latLonAlt2ECEF(np.radians(pos_i.lat), np.radians(pos_i.lon), pos_i.elev) 
            # B = latLonAlt2ECEF(np.radians(pos_f.lat), np.radians(pos_f.lon), pos_f.elev)

            # # Note: lat/lon here is incorrect
            # A = ecef2ENU(np.radians(90-pos_i.lat), np.radians(pos_i.lon-90), A[0], A[1], A[2])
            # B = ecef2ENU(np.radians(90-pos_f.lat), np.radians(pos_f.lon-90), B[0], B[1], B[2])


            azimuth = Angle(np.degrees(np.arctan2(dx, dy)))


            zenith = Angle(np.degrees(np.arctan2(dh, -dz)))

            self.vector = Vector3D(np.cos(azimuth.rad)*np.sin(zenith.rad), \
                                   np.sin(azimuth.rad)*np.sin(zenith.rad), \
                                                      -np.cos(zenith.rad))

            scale = dz / self.vector.z

        # Case 3: top point and angles are given
        elif pos_i.lat is not None and zenith is not None and azimuth is not None:

            END_HEIGHT = 11400

            if verbose:
                print("Case 3: Top point and angles are given")

            pos_i.pos_loc(pos_i)

            A = pos_i.xyz
            
            vector = Vector3D(np.cos(azimuth.rad)*np.sin(zenith.rad), \
                              np.sin(azimuth.rad)*np.sin(zenith.rad), \
                                                      -np.cos(zenith.rad))


            scale = (END_HEIGHT - A[2])/vector.z 

            B = A + vector.xyz*scale

            pos_f = Position(0, 0, 0)
            pos_f.x = B[0]
            pos_f.y = B[1]
            pos_f.z = B[2]
            pos_f.pos_geo(pos_i)

            self.vector = vector


        # Case 4: bottom point and angles are given
        elif pos_f.lat is not None and zenith is not None and azimuth is not None:

            START_HEIGHT = 81300

            if verbose:
                print("Case 4: Bottom point and angles are given")

            pos_f.pos_loc(pos_f)

            B = pos_f.xyz

            vector = Vector3D(np.cos(azimuth.rad)*np.sin(zenith.rad), \
                              np.sin(azimuth.rad)*np.sin(zenith.rad), \
                                                      -np.cos(zenith.rad))

            scale = (B[2] - START_HEIGHT)/vector.z

            A = B - vector.xyz*scale

            pos_i = Position(0, 0, 0)
            pos_i.x = A[0]
            pos_i.y = A[1]
            pos_i.z = A[2]
            pos_i.pos_geo(pos_f)

            self.vector = vector


        # Case 5: not enough information is given
        else:
            if verbose:
                print("Case 5: Not enough information is given")
            scale = None

        if pos_i.lat == None:
            pos_i = Position(None, None, None)

        if pos_f.lat == None:
            pos_f = Position(None, None, None)

        if azimuth == None:
            azimuth = Angle(None)

        if zenith == None:
            zenith = Angle(None)

        self.pos_i = pos_i
        self.pos_f = pos_f
        self.azimuth = azimuth
        self.zenith = zenith
        self.scale = scale        
        self.v = v

        if v_f is not None:
            self.v_f = v_f
            self.v_avg = np.mean((v, v_f))
        else:
            self.v_f = v
            self.v_avg = v

        # if using linear decel: speed goes from v -> v_f, v_avg is the mean
        # if using static vel: v == v_f == v_avg

    def __str__(self):
        degree_sign= 'o'

        try:
            A = "Trajectory: Lat: {:6.4f}{:} N, Lon: {:6.4f}{:} E, Elev: {:9.2f} km \n".format(self.pos_i.lat, \
                degree_sign, self.pos_i.lon, degree_sign, self.pos_i.elev/1000)
        except:
            A = "Trajectory: No known starting position \n"
        try:
            B = "    to      Lat: {:6.4f}{:} N, Lon: {:6.4f}{:} E, Elev: {:9.2f} km \n".format(self.pos_f.lat, \
                degree_sign, self.pos_f.lon, degree_sign, self.pos_f.elev/1000)
        except:
            B = "    to      No known ending position \n"

        try:
            C = "            Velocity: Initial {:4.1f} km/s \n".format(self.v_i/1000) + \
                "                        Final {:4.1f} km/s \n".format(self.v_f/1000)
        except:
            C = "            Velocity:  {:4.1f} km/s \n".format(self.v/1000)

        D = "            Time:      {:4.2f} s \n".format(self.t)

        try:
            E = "            Azimuth: {:6.4f}{:} \n".format(self.azimuth.deg, degree_sign)
        except:
            E = "            No known azimuth \n"

        try: 
            F = "            Zenith:  {:6.4f}{:} \n".format(self.zenith.deg, degree_sign)
        except:
            F = "            No known zenith \n"

        return A + B + C + D + E + F
    
    def getVelAtHeight(self, height):

        h_i = self.pos_i.elev
        h_f = self.pos_f.elev

        frac = (height - h_i)/(h_f - h_i)

        return frac*(self.v_f - self.v) + self.v       

    def getTrajVect(self):
        return np.array([np.sin(self.azimuth.rad)*np.sin(self.zenith.rad), \
                         np.cos(self.azimuth.rad)*np.sin(self.zenith.rad), \
                        -np.cos(self.zenith.rad)])




    def trajInterp2(self, div=250, min_p=None, max_p=None, xyz=False):
        """ Returns N points along a trajectory between two specified heights

        Arguments:
        div [int]       - N number of points to return (including both boundaries)
        min_p [float]   - Lower boundary to print out [meters from surface]
        max_p [float]   - upper boundary to print out [meters from surface]    
        xyz [Boolean]   - If True, return ECEF coordinates, else return Lat/Lon/Elev

        Returns:
        [lat, lon, elev, time] - [4, div] array of latitude, longitude, elevation and 
        time of each point given in deg, deg, meters, seconds.

        if xyz is set to True, then lat, lon elev will be replaced by x, y, z in meters
        in the ECEF coordinte system
        """
    
        # Use this to go past the trajectory limits
        SCALER = 1

        # Use trajectory limits if not given
        if min_p is None:
            min_p = self.pos_f.elev

        if max_p is None:
            max_p = self.pos_i.elev

        pos_i = self.findGeo(max_p)
        pos_f = self.findGeo(min_p)


        # Convert trajectory into local coords
        pos_i.pos_loc(pos_f)
        pos_f.pos_loc(pos_f)

        A = pos_i.xyz
        B = pos_f.xyz

        # A = latLonAlt2ECEF(np.radians(self.pos_i.lat), np.radians(self.pos_i.lon), self.pos_i.elev) 
        # B = latLonAlt2ECEF(np.radians(self.pos_f.lat), np.radians(self.pos_f.lon), self.pos_f.elev)

        if div <= 2:
            P = []
            if xyz:
                P.append([A[0], A[1], A[2], self.findTime(A[2])])
                P.append([B[0], B[1], B[2], self.findTime(B[2])])   
            else:
                P.append([self.pos_i.lat, self.pos_i.lon, self.pos_i.elev, self.findTime(self.pos_i.elev)])
                P.append([self.pos_f.lat, self.pos_f.lon, self.pos_f.elev, self.findTime(self.pos_f.elev)])

            return P

        # Points Array
        P = [None]*(int(SCALER*div))

        for i in range(int(SCALER*div)):

            # Initialize position object, lat/lon/elev will be reset
            P[i] = Position(0, 0, 0)

            # interp positions
            P[i].x = A[0] + i*(B[0] - A[0])/(div-1)
            P[i].y = A[1] + i*(B[1] - A[1])/(div-1)
            P[i].z = A[2] + i*(B[2] - A[2])/(div-1)

            # Convert to lat/lon/elev if needed
            P[i].pos_geo(self.pos_f)

        # Organize Arrays
        A = []
        for pt in P:
            if xyz:
                A.append([pt.x, pt.y, pt.z, self.findTime(pt.z)])
            else:
                A.append([pt.lat, pt.lon, pt.elev, self.findTime(pt.elev)])


        return np.array(A)


    def findGeo(self, height):
        """ Returns the coordinates (lat/lon/elev) of the trajectory at a 
        specific height

        Arguments:
        height [float] - height to return lat/lon of [meters]

        Returns:
        Position Object with lat, lon, elev at the height
        """

        # Convert trajectory into ECEF

        self.pos_i.pos_loc(self.pos_f)
        self.pos_f.pos_loc(self.pos_f)

        A = self.pos_i.xyz
        B = self.pos_f.xyz

        A = np.array(A)
        B = np.array(B)

        frac = (height - self.pos_i.elev)/(self.pos_f.elev - self.pos_i.elev) 

        C = A + frac*(B - A)

        geo = Position(0, 0, 0)
        geo.x = C[0]
        geo.y = C[1]
        geo.z = C[2]

        geo.pos_geo(self.pos_f)

        return geo

    def verticalVel(self, v=None):
        """ Returns the z component of the velocity
        """
        if v is None:
            v_h = self.v*np.cos(self.zenith.rad)
        else:
            v_h = v*np.cos(self.zenith.rad)

        return v_h


    def approxHeight(self, time, v_list=[]):
        """ Returns the height of the meteoroid at a time in relation to the same 
        reference time assuming no deceleration
        """

        v_list = np.array([[self.pos_i.elev, 17.600],
                          [60000, 17.600],
                          [40200, 17.503],
                          [35500, 17.411],
                          [26400, 16.780],
                          [20000, 15.090],
                          [17900, 13.731],
                          [16690, 13.000],
                          [14900, 11.118],
                          [13300,  9.598],
                          [12200,  7.345],
                          [11100,  6.000],
                          [ 9500,  3.000],
                          [    0,  3.000]])

        # constant velocity
        if len(v_list) == 0:

            dt = time - self.t

            v_h = self.verticalVel()

            height = self.pos_i.elev - dt*v_h

            return height

        else:

            target_dt = 0
            target_h = self.pos_i.elev

            # time from beginning to target height
            dt = time - self.t

            h_list = v_list[:, 0]
            vs = v_list[:, 1]*1000

            dt_list = []

            for ii in range(1, len(v_list)):

                vh = (self.verticalVel(v=vs[ii-1]) + self.verticalVel(v=vs[ii]))/2
                dh = abs(h_list[ii] - h_list[ii - 1])


                if target_dt + dh/vh > dt:

                    frac = (dt - target_dt)/(dh/vh)

                    final_h = frac*(-dh) + target_h
                    # print("Final Answer: Time: {:}, Height: {:.2f}, Vertical Velocity: {:.2f}, Layer: {:}".format(dt, final_h, vh, ii))
                    return final_h

                else:

                    target_dt += dh/vh
                    target_h -= dh

                # print("Temp: Time: {:}, Height: {:.2f}, Vertical Velocity: {:.2f}, Layer: {:}".format(target_dt, target_h, vh, ii))

            return None



    def findTime(self, height, v_list=[]):
        """ Returns the time that it takes for the meteor to get from its inital position
        to a given height. Assumes that the bolide travels at a constant speed

        Arguments:
        height [float] - the height to find the time at [meters]

        Returns:
        time [float] - the time at that point [seconds]
        """

        # constant time
        if len(v_list) == 0:

            # Vertical component of velocity
            v_h = self.verticalVel()

            dh = self.pos_i.elev - height

            # Get time travelled
            t = dh/v_h + self.t

            return t



    def findLength(self):
        """ Returns the length of the meteor from given initial and final positions
        in meters
        """

        # Convert trajectory into ECEF
        A = latLonAlt2ECEF(self.pos_i.lat_r, self.pos_i.lon_r, self.pos_i.elev) 
        B = latLonAlt2ECEF(self.pos_f.lat_r, self.pos_f.lon_r, self.pos_f.elev)

        length_of_meteor = np.sqrt((A[0] - B[0])**2 + (A[1] - B[1])**2 + (A[2] - B[2])**2)

        return length_of_meteor

    def getRanges(self, stat_pos, div=100, write=False):

        P = self.trajInterp(div=div)

        pos_list = []
        for element in P:
            r = stat_pos.pos_distance(element)
            pos_list.append(r)

            if write:
                print(element, r)

        return pos_list

class Plane:
    def __init__(self, p1, p2, p3):
        
        n = np.cross((p2 - p1), (p3 - p1))
        self.n = n/np.sqrt(np.dot(n, n))
        self.x0 = p3

    def __str__(self):

        return "Plane object with normal: {:}".format(self.n)

    def checkOnPlane(self, x, tol=0):
    
        if np.abs(np.dot(self.n, (x - self.x0))) <= tol:
            return True
        else:
            return False

class Color:
    def __init__(self):
        self.nominal = (255, 0, 238)
        self.both = (255, 0, 0)
        self.ballistic = (0, 0, 255)
        self.fragmentation = [(0, 255, 0)]
        self.perturb = [([(0, 255, 26, 150), (3, 252, 176, 150), (252, 3, 3, 150), (176, 252, 3, 150), (255, 133, 3, 150),
                        (149, 0, 255, 150), (76, 128, 4, 150), (82, 27, 27, 150), (101, 128, 125, 150), (5, 176, 249, 150)])]
        self.BLACK = (0, 0, 0)
        self.WHITE = (255, 255, 255)

    def generate(self):
        r = np.random.randint(0, 255)
        g = np.random.randint(0, 255)
        b = np.random.randint(0, 255)
        return (r, g, b)

    def generateLight(self):
        r = np.random.randint(150, 255)
        g = np.random.randint(150, 255)
        b = np.random.randint(150, 255)
        return (r, g, b)

class RectangleItem(pg.GraphicsObject):
    def __init__(self, data, c_map="forward", alpha=70):
        pg.GraphicsObject.__init__(self)

        self.c_map = c_map

        self.data = data
        
        if self.c_map != "set":
            self.normed_dat = self.normData()
        
        self.generatePicture(alpha=alpha)

    
    def normData(self):
        a = np.array(self.data)

        a = a[:, 4]
        
        max_data = np.nanmax(a)
        min_data = np.nanmin(a)

        r = max_data - min_data 

        a = (a - min_data)/r

        return a

    def weightedNormData(self):
        a = np.array(self.data)
        a = a[:, 4]
        
        max_data = np.nanmax(a)
        min_data = np.nanmin(a)

        r = max_data - min_data 

        a = (a - min_data)/r

        return a


    def gradient(self, weight, alpha):
        if weight == weight:
            if self.c_map == "forward":
                c_ender = (7, 65, 112)
                c_end = (111, 255, 0)
                c_start   = (224, 245, 66)
                c_starter = (232, 72, 9)
            elif self.c_map == "reverse":
                c_starter = (7, 65, 112)
                c_start = (111, 255, 0)
                c_end   = (224, 245, 66)
                c_ender = (232, 72, 9)
            elif self.c_map == "white":
                c_starter = (255, 255, 255)
                c_start =   (255, 255, 255)
                c_end   =   (255, 255, 255)
                c_ender =   (255, 255, 255)                

            if weight > (2/3):
                r = c_end[0] + 3*(weight - (2/3))*(c_ender[0] - c_end[0])
                g = c_end[1] + 3*(weight - (2/3))*(c_ender[1] - c_end[1])
                b = c_end[2] + 3*(weight - (2/3))*(c_ender[2] - c_end[2])
            elif weight > (1/3):
                r = c_start[0] + 3*(weight - (1/3))*(c_end[0] - c_start[0])
                g = c_start[1] + 3*(weight - (1/3))*(c_end[1] - c_start[1])
                b = c_start[2] + 3*(weight - (1/3))*(c_end[2] - c_start[2])
            else:
                r = c_starter[0] + 3*(weight)*(c_start[0] - c_starter[0])
                g = c_starter[1] + 3*(weight)*(c_start[1] - c_starter[1])
                b = c_starter[2] + 3*(weight)*(c_start[2] - c_starter[2])

            return (r, g, b, alpha)
        else:
            return (0, 0, 0, 0)

    def generatePicture(self, alpha=70):
        ## pre-computing a QPicture object allows paint() to run much more quickly, 
        ## rather than re-drawing the shapes every time.
        self.picture = QtGui.QPicture()
        p = QtGui.QPainter(self.picture)

        for ii, (c_x, c_y, h, w, o) in enumerate(self.data):
            
            if self.c_map == "set":
                p.setBrush(pg.mkBrush(o))
                p.setPen(pg.mkPen(o))
            else:
                z = self.normed_dat[ii]
                p.setBrush(pg.mkBrush(self.gradient(z, alpha)))
                p.setPen(pg.mkPen(self.gradient(z, 0)))
            l = c_x - w/2
            t = c_y + h/2
            p.drawRect(QtCore.QRectF(l, t, w, h))
        p.end()
    
    def paint(self, p, *args):
        p.drawPicture(0, 0, self.picture)
    
    def boundingRect(self):
        ## boundingRect _must_ indicate the entire area that will be drawn on
        ## or else we will get artifacts and possibly crashing.
        ## (in this case, QPicture does all the work of computing the bouning rect for us)
        return QtCore.QRectF(self.picture.boundingRect())

if __name__ == '__main__':


    print("### Trajectory Testing ###")
    print("# Check parameters of trajectory are the same")
    pos_i = Position(43.085167, -80.726829, 63910)
    pos_f = Position(43.192440, -79.446173, 21799)

    A = Trajectory(0, 17400, pos_i=pos_i, pos_f=pos_f, verbose=True)
    B = Trajectory(0, 17400, pos_i=pos_i, azimuth=A.azimuth, zenith=A.zenith, verbose=True)
    C = Trajectory(0, 17400, pos_f=pos_f, azimuth=A.azimuth, zenith=A.zenith, verbose=True)
    # times = [0, 1, 2, 3, 4]
    # v_list = np.array([[60000, 17.385],
    #                   [40200, 17.306],
    #                   [35500, 17.214],
    #                   [26400, 16.576],
    #                   [17900, 14.11],
    #                   [14900, 11.97],
    #                   [13300, 10.27],
    #                   [11100, 6.50],
    #                   [ 9000,  2.96],
    #                   [    0,  0.00]])

    # for tt in times:
    #     A.approxHeight(tt, v_list)


    # B = Trajectory(0, 17400, pos_i=pos_i, azimuth=A.azimuth, zenith=A.zenith, verbose=True)
    # C = Trajectory(0, 17400, pos_f=pos_f, azimuth=A.azimuth, zenith=A.zenith, verbose=True)
    

    from tabulate import tabulate
    
    print(tabulate([["Pos_i - lat", A.pos_i.lat, B.pos_i.lat, C.pos_i.lat], \
                    ["Pos_i - lon", A.pos_i.lon, B.pos_i.lon, C.pos_i.lon], \
                    ["Pos_i - elev", A.pos_i.elev, B.pos_i.elev, C.pos_i.elev], \
                    ["Pos_f - lat", A.pos_f.lat, B.pos_f.lat, C.pos_f.lat], \
                    ["Pos_f - lon", A.pos_f.lon, B.pos_f.lon, C.pos_f.lon], \
                    ["Pos_f - elev", A.pos_f.elev, B.pos_f.elev, C.pos_f.elev], \
                    ["Azimuth", A.azimuth.deg, B.azimuth.deg, C.azimuth.deg], \
                    ["Zenith", A.zenith.deg, B.zenith.deg, C.zenith.deg]], \
                    headers=["Case 2", "Case 3", "Case 4"]))


    



  