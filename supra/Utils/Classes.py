import numpy as np

from supra.Utils.AngleConv import geo2Loc, loc2Geo, angle2NDE

class Config:

    def __init__(self):

        # set configuration defaults
        self.fireball_name = 'Untitled Fireball'
        self.get_data = 'false'
        self.run_mode = 'search'
        self.debug = 'false'

        self.working_directory = None
        self.arrival_times_file = None
        self.sounding_file = None
        self.perturbation_spread_file = None
        self.station_picks_file = None
        self.points_name = None

        self.lat_centre = None
        self.lon_centre = None
        self.deg_radius = 2
        self.fireball_datetime = None
        self.v_sound = 310

        self.t0 = None
        self.v = None
        self.azim = None
        self.zangle = None
        self.lat_i = None
        self.lon_i = None
        self.elev_i = None
        self.lat_f = None
        self.lon_f = None
        self.elev_f = None
        self.show_ballistic_waveform = 'false'

        self.fragmentation_point = None
        self.show_fragmentation_waveform = 'false'
        self.manual_fragmetnation_search = None

        self.v_fixed = 11000
        self.azimuth_min = 0
        self.azimuth_max = 359.99
        self.zenith_min = 0
        self.zenith_max = 89.99
        self.x_min = -200000
        self.x_max = 200000
        self.y_min = -200000
        self.y_max = 200000
        self.z_min = 0
        self.z_max = 100000
        self.t_min = -200
        self.t_max = 200
        self.v_min = 11000
        self.v_max = 30000 
        self.max_error = 1000
        self.restricted_time = None
        self.enable_restricted_time = 'false'
        self.weight_distance_min = 0
        self.weight_distance_max = 0

        self.enable_winds = 'false'
        self.weather_type = 'none'

        self.perturb = 'false'
        self.perturb_method = 'none'
        self.perturb_times = 0
        self.observe_frag_no = 0

        self.n_theta = 45
        self.n_phi = 90
        self.angle_precision = 1e-5
        self.angle_error_tol = 1000

        self.maxiter = 100
        self.swarmsize = 100
        self.run_times = 1
        self.phip = 0.5
        self.phig = 0.5
        self.omega = 0.5
        self.pso_debug = 'false'
        self.minfunc = 1e-8
        self.minstep = 1e-8

        self.contour_res = 10
        self.high_f = None
        self.high_b = None
        self.rm_stat = None

        self.stations = None

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

class Pick:
    def __init__(self, t, stn, stn_no, channel):
        self.time = t
        self.stn = stn
        self.stn_no = stn_no
        self.channel = channel


class Angle:

    def __init__(self, deg):
        
        if deg != None:
            self.deg = deg%360
            self.rad = np.radians(deg)
        else:
            self.deg = None
            self.rad = None

    def __str__(self):
        return '{:6.2f} deg'.format(self.deg) 

    def __add__(self, other):
        return (self.deg + other.deg)%360

    def __sub__(self, other):
        return (self.deg - other.deg)%360

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

        if lat != None:
            self.lat_r = np.radians(lat)
        else:
            self.lat_r = None

        if lon != None:  
            self.lon_r = np.radians(lon)
        else:
            self.lon_r = None

    def __str__(self):
        try:
            result = "Lat: {:8.4f} deg N, Lon: {:8.4f} deg E, Elev: {:10.2f} m".format(self.lat, self.lon, self.elev)
        except:
            result = "Position is None Type"
        return result
    
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

        self.pos_loc(self)
        other.pos_loc(self)

        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2 + (self.z - other.z)**2)

    def ground_distance(self, other):
        
        self.pos_loc(self)
        other.pos_loc(self)

        return np.sqrt((self.x - other.x)**2 + (self.y - other.y)**2)

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
        return '[ {:.2f}, {:.2f}, {:.2f}]'.format(self.x, self.y, self.z)


class Station:

    """
    A station object containing information and position of the station
    """

    def __init__(self, network, code, position, channel, name, file_name):
        """
        Arguments:
            network: [String] the network of the station
            code: [String] the code of the station
            position: [position Obj] a position object of the station containing the lat/lon/elev of the station
            channel: [String] the channel of the station (BHZ, HHZ, BDF, etc.)
            name: [String] A string describing the station's name or location (for display purposes only)
            file_name: [String] Location of the .mseed file inside the working directory
        """

        self.network = network
        self.code = code
        self.position = position
        self.channel = channel
        self.file_name = file_name
        self.name = name
        self.offset = 0

    def __str__(self):
        A = "Station: {:} | Channel {:} \n".format(self.name, self.channel)
        B = str(self.position) + "\n"
        try:
            C = "Ground Distance: {:}".format(self.ground_distance)
        except:
            C = ''
        return A + B + C

    def stn_distance(self, ref_pos):
        self.distance = self.position.pos_distance(ref_pos)

    def stn_ground_distance(self, ref_pos):
        self.ground_distance = self.position.ground_distance(ref_pos)

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
    def __init__(self, t, v, zenith=None, azimuth=None, pos_i=None, pos_f=None):
        
        self.v = v

        self.t = t

        # Find vector from angles
        if zenith != None and azimuth != None:

            self.vector = Vector3D(np.sin(azimuth.rad)*np.sin(zenith.rad), \
                                   np.cos(azimuth.rad)*np.sin(zenith.rad), \
                                                      -np.cos(zenith.rad))
        else:

            self.vector = None

        # Case 1: everything is given
        if pos_i.lat != None and pos_f.lat != None and zenith != None and azimuth != None:

            pos_i.pos_loc(pos_f)
            pos_f.pos_loc(pos_f)

            scale = (pos_f.z - pos_i.z) / self.vector.z

        # Case 2: end points are given
        elif pos_i.lat != None and pos_f.lat != None:

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
        elif pos_i.lat != None and zenith != None and azimuth != None:

            pos_i.pos_loc(pos_i)

            scale = -pos_i.z / self.vector.z

            pos_f = self.vector * scale + pos_i

            pos_f.pos_geo(pos_i)

        # Case 4: bottom point and angles are given
        elif pos_f.lat != None and zenith != None and azimuth != None:

            pos_f.pos_loc(pos_f)

            scale = -100000 / self.vector.z

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


    def __str__(self):

        try:
            A = "Trajectory: Lat: {:6.2f} deg N, Lon: {:6.2f} deg E, Elev: {:9.2f} m \n".format(self.pos_i.lat, self.pos_i.lon, self.pos_i.elev)
        except:
            A = "Trajectory: No known starting position \n"
        try:
            B = "    to      Lat: {:6.2f} deg N, Lon: {:6.2f} deg E, Elev: {:9.2f} m \n".format(self.pos_f.lat, self.pos_f.lon, self.pos_f.elev)
        except:
            B = "    to      No known ending position \n"

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
                print(pt)

        return P

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

        self.pos_i.pos_loc(self.pos_f)
        self.pos_f.pos_loc(self.pos_f)

        frac = height / (self.pos_i.z)

        A = self.pos_i
        B = self.pos_f

        length_of_meteor = np.sqrt((A.x - B.x)**2 + (A.y - B.y)**2 + (A.z - B.z)**2)
        time_of_meteor = length_of_meteor/self.v

        return frac*time_of_meteor

if __name__ == '__main__':
    A = Trajectory(0, 13913, pos_i=Position(48.05977, 13.10846, 85920.0), pos_f=Position(48.3314, 13.0706, 0))
    #A = Trajectory(0, 13913, pos_i=Position(-23.5742415562, 132.712445759, 100000), pos_f=Position(-23.616963, 132.902681, 0))
    A.trajInterp(div=250, write=True)