# -*- coding: utf-8 -*-
import numpy as np

from supra.Utils.AngleConv import geo2Loc, loc2Geo, angle2NDE

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore

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

        self.k = 2e-4

        self.b = 1.19e-4

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
        degree_sign= u'\N{DEGREE SIGN}'
        try:
            result = "Lat: {:8.4f}{:}N, Lon: {:8.4f}{:}E, Elev: {:10.2f} m".format(self.lat, degree_sign, self.lon, degree_sign, self.elev)
        except TypeError:
            result = "Position is None Type"
        return result

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

        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)

    def ground_distance(self, other):
        """ 2D horizontal distance between positions 'self' and 'other'
        """

        self.pos_loc(self)
        other.pos_loc(self)

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


class Annote:

    def __init__(self, title, time, length, group, source, height, notes, color):

        self.title = title
        self.time = time
        self.length = length
        self.group = group
        self.source = source
        self.height = height
        self.notes = notes
        self.color = color

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
    def __init__(self, t, v, zenith=None, azimuth=None, pos_i=Position(None, None, None), pos_f=Position(None, None, None), v_f=None):
        
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

            pos_i.pos_loc(pos_f)
            pos_f.pos_loc(pos_f)

            scale = (pos_f.z - pos_i.z) / self.vector.z

        # Case 2: end points are given
        elif pos_i.lat is not None and pos_f.lat is not None:

            pos_i.pos_loc(pos_f)
            pos_f.pos_loc(pos_f)

            azimuth = Angle(angle2NDE(np.degrees(np.arctan2(pos_f.y - pos_i.y, pos_f.x - pos_i.x))))

            ground_distance = (np.sqrt((pos_f.y - pos_i.y)**2 + (pos_f.x - pos_i.x)**2))

            zenith = Angle(np.degrees(np.arctan2(ground_distance,-(pos_f.z - pos_i.z))))

            self.vector = Vector3D(np.sin(azimuth.rad)*np.sin(zenith.rad), \
                                   np.cos(azimuth.rad)*np.sin(zenith.rad), \
                                                      -np.cos(zenith.rad))

            scale = (pos_f.z - pos_i.z) / self.vector.z

        # Case 3: top point and angles are given
        elif pos_i.lat is not None and zenith is not None and azimuth is not None:

            pos_i.pos_loc(pos_i)

            scale = -pos_i.z / self.vector.z

            pos_f = self.vector * scale + pos_i

            pos_f.pos_geo(pos_i)

        # Case 4: bottom point and angles are given
        elif pos_f.lat is not None and zenith is not None and azimuth is not None:

            pos_f.pos_loc(pos_f)

            scale = -50000 / self.vector.z

            pos_i = self.vector * -scale + pos_f 

            pos_i.pos_geo(pos_f)

        # Case 5: not enough information is given
        else:
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

        try:
            A = "Trajectory: Lat: {:6.2f} deg N, Lon: {:6.2f} deg E, Elev: {:9.2f} m \n".format(self.pos_i.lat, self.pos_i.lon, self.pos_i.elev)
        except:
            A = "Trajectory: No known starting position \n"
        try:
            B = "    to      Lat: {:6.2f} deg N, Lon: {:6.2f} deg E, Elev: {:9.2f} m \n".format(self.pos_f.lat, self.pos_f.lon, self.pos_f.elev)
        except:
            B = "    to      No known ending position \n"

        try:
            C = "            Velocity: Initial {:4.1f} km/s \n".format(self.v_i/1000) + \
                "                        Final {:4.1f} km/s \n".format(self.v_f/1000)
        except:
            C = "            Velocity:  {:4.1f} km/s \n".format(self.v/1000)

        D = "            Time:      {:4.2f} s \n".format(self.t)

        try:
            E = "            Azimuth: {:6.2f} deg \n".format(self.azimuth.deg)
        except:
            E = "            No known azimuth \n"

        try: 
            F = "            Zenith:  {:6.2f} deg \n".format(self.zenith.deg)
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

    def trajInterp(self, div=10, write=False):

        self.pos_i.pos_loc(self.pos_f)
        self.pos_f.pos_loc(self.pos_f)

        v = self.vector * self.scale

        dv = v * (1/(div-1))

        P = [None]*div
        for i in range(div):

            P[i] = Position(0, 0, 0)

            P[i].x = self.pos_i.x + i*dv.x
            P[i].y = self.pos_i.y + i*dv.y
            P[i].z = self.pos_i.z + i*dv.z

            P[i].pos_geo(self.pos_f)
        
        if write:
            for pt in P:
                print(pt, self.findTime(pt.elev))

        return P

    def trajInterp2(self, div=250, min_p=17000, max_p=50000):
        
        ref_loc = self.pos_f
        ref_loc.elev = 0

        self.pos_i.pos_loc(ref_loc)
        self.pos_f.pos_loc(ref_loc)

        v = self.vector * self.scale

        dv = v * (1/(div-1))

        P = [None]*div
        for i in range(div):

            P[i] = Position(0, 0, 0)

            P[i].x = self.pos_i.x + i*dv.x
            P[i].y = self.pos_i.y + i*dv.y
            P[i].z = self.pos_i.z + i*dv.z

            P[i].pos_geo(self.pos_f)
            
        A = []
        for pt in P:
            if min_p <= pt.elev <= max_p:
                A.append([pt.lat, pt.lon, pt.elev, self.findTime(pt.elev)])

        return np.array(A)
    
    def findPoints(self, gridspace=250, min_p=0, max_p=0):    
        
        GRID_SPACE = gridspace
        
        if min_p == 0:
            MIN_HEIGHT = self.pos_f.elev
        else:
            MIN_HEIGHT = min_p
        if max_p == 0:
            MAX_HEIGHT = self.pos_i.elev
        else:
            MAX_HEIGHT = max_p

        u = self.vector.xyz

        ground_point = np.array([0, 0, 0])
        # find top boundary of line given maximum elevation of trajectory
        if self.pos_i.elev != None:
            scale = -self.pos_i.elev/u[2]

        else:
            scale = -50000/u[2]

        # define line top boundary
        top_point = ground_point - scale*u

        ds = scale / (GRID_SPACE)

        points = []

        for i in range(GRID_SPACE + 1):
            points.append(top_point + i*ds*u)

        points = np.array(points)

        offset = np.argmin(np.abs(points[:, 2] - MAX_HEIGHT))
        bottom_offset = np.argmin(np.abs(points[:, 2] - MIN_HEIGHT))

        points = np.array(points[offset:(bottom_offset+1)])

        return points

    def findGeo(self, height):
        
        self.pos_i.pos_loc(self.pos_f)
        self.pos_f.pos_loc(self.pos_f)

        C = Position(0, 0, 0)
        
        C.z = height
        s = (C.z - self.pos_i.z) / self.vector.z

        C.x = s*self.vector.x + self.pos_i.x
        C.y = s*self.vector.y + self.pos_i.y

        C.pos_geo(self.pos_f)

        return C

    def findTime(self, height):

        
        frac = (height - self.pos_f.elev) / (self.pos_i.elev - self.pos_f.elev)

        # v = (self.v_f + self.getVelAtHeight(height))/2
        v = self.v
        d = (1-frac)*(self.findLength())

        vect = self.getTrajVect()
        s = -self.pos_f.elev*vect[2]
        dist_vect = s*vect

        d_ground = np.sqrt(dist_vect[0]**2 + dist_vect[1]**2 + dist_vect[2]**2) 

        t = d/v + d_ground/self.v_f

        return self.t + t



    def findLength(self):

        self.pos_i.pos_loc(self.pos_f)
        self.pos_f.pos_loc(self.pos_f)

        A = self.pos_i
        B = self.pos_f

        length_of_meteor = np.sqrt((A.x - B.x)**2 + (A.y - B.y)**2 + (A.z - B.z)**2)
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

# Using American spelling to avoid confusion with other "colour" related tools, sorry
# This program is proudly Canadian
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


    # S = Position(48.8461, 13.7179, 0)

    # print("Line")
    # D.line_between_2d(S, show_error=True)
    # print("Dist", D.ground_distance(S))    
    # m = 1.483894107652387
    # b = -13062.851161267856
    # R = 73460.63
    # print((-m*b + np.sqrt(m**2*R**2 + b**2 - R**2))/(m**2 + 1))

    # S_p = Position(0, 0, 0)
    # S_m = Position(0, 0, 0)
    # S_p.x = 42122.68
    # S_p.y = 49443.00
    # S_p.z = 0
    # S_m.x = 54260.14
    # S_m.y = 67453.77
    # S_m.z = 0

    # S_p.pos_geo(ref)
    # S_m.pos_geo(ref)

    # print(S_p)
    # print(S_m)Begin point on the trajectory:
    pos_i = Position(48.1724466606, 13.0926245672, 50000)
    A = Trajectory(2.471993030094728, 13913.0, pos_i=pos_i, \
                             zenith=Angle(85), azimuth=Angle(354.67))

    s = Position(48.8461, 13.7179,   1141.10)

    wrp = A.findGeo(43288.5906040269)
    rag = wrp.pos_distance(s)
    print(rag)
    # import numpy as np
    # import matplotlib.pyplot as plt

    # import csv

    # with open('C:\\Users\\lmcfd\\Desktop\\Theoretical\\aaa.csv', newline='') as f:
    #     reader = csv.reader(f)
    #     data = list(reader)

    # xs = []
    # ys = []
    # current_stn = "NORI"
    # for line in data:
    #     if float(line[3]) > 50:
    #         line[0] = line[0].replace(u'\xa0', u'')
    #         if line[0] == current_stn:
                
    #             stn_no = None
    #             for i in range(len(stn)):

    #                 if line[0] == stn[i][0]:
    #                     stn_no = i


    #             zenith = float(line[4])
    #             height = float(line[2])
    #             spread = float(line[1])/2

    #             stat_pos = Position(stn[stn_no][1], stn[stn_no][2], stn[stn_no][3])
    #             pos_i = Position(48.1724466606, 13.0926245672, 50000)
    #             A = Trajectory(0, 11, pos_i=pos_i, zenith=Angle(zenith), azimuth=Angle(354.67))

    #             wrp = A.findGeo(height)
    #             rag = wrp.pos_distance(stat_pos)
                
    #             xs.append(rag)
    #             ys.append(spread)

    # plt.scatter(xs, ys)

    # plt.show()



    # A.trajInterp(div=200, write=True)
    # print(A.findGeo(34500))
    # print(A.findTime(34500))
    # A = Position()
    #A = Trajectory(0, 13913, pos_i=Position(48.2042, 13.0884, 39946.8), pos_f=Position(48.8461, 13.7179, 1098))
    # A = Trajectory(0, 13913, pos_i=Position(48.05977, 13.10846, 85920.0), pos_f=Position(48.3314, 13.0706, 0))
    #A = Trajectory(0, 13913, pos_f=Position(43.6783, -80.0647, 72386), pos_i=Position(43.5327, -80.2335, 95901))
    #A = Trajectory(0, 13913, pos_i=Position(42.0142, -81.2926, 100486), \
    #                         pos_f=Position(41.9168, -81.3655,  71666))
    # A.getRanges(Position(48.846135, 13.71793, 1141.1), write=True, div=1000)
    #A = Trajectory(0, 13913, pos_i=Position(-23.5742415562, 132.712445759, 100000), pos_f=Position(-23.616963, 132.902681, 0))
    # print(A.vector)
    #A = Trajectory(0, 13913, pos_i=Position(42.8408, -81.9617, 119216), pos_f=Position(43.1930, -81.3163, 298))
    #print(A)
    #A = Trajectory(0, 13913, pos_i=Position(42.8408, -81.9617, 119216), pos_f=Position(43.7675, -82.4195, 88458))

    # P = A.trajInterp(div=500, write=True)
    # A = Trajectory(0, 13913, pos_i=Position(43.721, -78.680, 95536), pos_f=Position(44.734, -78.212, 29307))#pos_f=Position(44.828, -78.153, 28866))
    # print(A.azimuth)