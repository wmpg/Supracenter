import random
import numpy as np
from netCDF4 import Dataset

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from scipy.interpolate import CubicSpline

from supra.Utils.AngleConv import roundToNearest
from supra.Utils.Classes import Constants
from supra.Supracenter.cyzInteg import zInteg
from supra.GUI.Tools.GUITools import *
from supra.Atmosphere.NRLMSISE import getAtmDensity
from wmpl.Utils.TrajConversions import date2JD
from supra.Atmosphere.HWM93 import getHWM
consts = Constants()

class AtmosType:
    def __init__(self):
        pass

    def interp(self, lat, lon, div=0.25):
        """
        Approximately interpolates grid points of a division between point A and point B

        lat: [lat_initial, lat_final]
        lon: [lon_initial, lon_final]
        """

        # Collect interpolated points
        pts = []

        x0 = lat[0]
        x1 = lat[1]
        y0 = lon[0]
        y1 = lon[1]

        dx = np.abs(x1 - x0)/div
        dy = np.abs(y1 - y0)/div

        # Approximate upper limit on number of steps
        try:
            steps = np.abs(int((dx - 2)*(dy - 2) + 2)//3 + 1)
        except ValueError:
            steps = 10

        x = np.linspace(x0, x1, steps)
        y = np.linspace(y0, y1, steps)

        # Take all steps and lock them to the nearest grid point
        for i in range(steps):
            pts.append([roundToNearest(x[i], div), roundToNearest(y[i], div)])

        if len(pts) == 0:
            return [[roundToNearest(x0, div), roundToNearest(y0, div)], [roundToNearest(x1, div), roundToNearest(y1, div)]]

        return pts


    def convert(self, t, u, v, z):
        """
        Converts temp, u-wind, v-wind, geopotential to
        height, sp of sound, wind magnitude and wind direction
        """

        speed = np.sqrt(consts.GAMMA*consts.R/consts.M_0*t)
        mags = np.sqrt(u**2 + v**2)
        dirs = np.arctan2(u, v)
        level = (z - z[-1])/consts.g_0

        return level, speed, mags, dirs

    def spline(self, level, speed, mags, dirs, interp=100):
        """
        Cubic spline the data for smooth profiles
        """

        speed, _ = self.splineSec(speed, level, interp=interp)
        mags, _ = self.splineSec(mags, level, interp=interp)
        dirs, level = self.splineSec(dirs, level, interp=interp)

        return level, speed, mags, dirs

    def splineSec(self, y, x, interp=100):
        
        if interp is None:
            return y, x

        x, y = np.flip(x), np.flip(y)
        f = CubicSpline(x, y)

        new_x = np.linspace(x[0], x[-1], interp-1)
        new_y = f(new_x)

        return new_y, new_x

class DefaultW(AtmosType):
    def __init__(self):

        pass

    def genProfile(self, lat, lon, heights, prefs, spline, ref_time):

        dh = 1000 #meters

        missing_heights = np.arange(heights[0], heights[1], -dh)
        missing_lats = np.linspace(lat[0], lat[1], len(missing_heights))
        missing_lons = np.linspace(lon[0], lon[1], len(missing_heights))

        if len(missing_heights) == 0:
            if heights[0] > heights[1]:
                missing_heights = np.array(heights)
                missing_lats = np.array(lat)
                missing_lons = np.array(lon)
            else:
                missing_heights = np.array([heights[1], heights[0]])
                missing_lats = np.array([lat[1], lat[0]])
                missing_lons = np.array([lon[1], lon[0]])  

        jd = date2JD(ref_time.year, ref_time.month, ref_time.day, ref_time.hour, ref_time.minute, ref_time.second)

        level = []
        speed = []
        mags = []
        dirs = []


        for hh, la, lo in zip(missing_heights, missing_lats, missing_lons):

            t = getAtmDensity(la, lo, hh, jd)
            speed.append(np.sqrt(consts.GAMMA*consts.R/consts.M_0*t))
            u, v = getHWM(ref_time, la, lo, hh/1000)
            mags.append(np.sqrt(u**2 + v**2))
            dirs.append(np.arctan2(u, v))
            level.append(hh)

        try:
            level, speed, mags, dirs = self.spline(level, speed, mags, dirs, interp=spline)

            sounding = []
            for i in range(len(level)):
                sounding.append([level[i], speed[i], mags[i], dirs[i]])
        except ValueError:

            sounding =         np.array([[    0.0, 310, 0.0, 0.0],
                                         [    0.0, 310, 0.0, 0.0],
                                         [99999.0, 310, 0.0, 0.0]])


        
        sounding = zInteg(heights[0], heights[1], np.array(sounding))

        return sounding


class ECMWF(AtmosType):
    def __init__(self, lat, lon, rng, time, file_name):
        
        # the grided subsection of the world
        self.generateECMWF(lat, lon, rng, time, file_name)

    def readNCDF(self, lat, lon, rng, time, file_name, div):
        dataset = Dataset(file_name, "r+", format="NETCDF4")

        hour = time.time().hour()

        if time.time().minute() > 30:
            hour += 1

        time = np.array(dataset.variables['time'])
        longitude = np.array(dataset.variables['longitude'])
        latitude = np.array(dataset.variables['latitude'])

        lat = roundToNearest(lat, div)
        lon = roundToNearest(lon, div)%360

        try:
            lon_index = int(np.where(longitude==lon)[0])
            lat_index = int(np.where(latitude==lat)[0])
        except TypeError as e:
            # For when we can't find the right longitude and latitude in the data (Thanks Jouse for finding this!)
            errorMessage('Unable to find exact lat/lon in given weather profile, using closest values!', 1, detail='{:}'.format(e))
            lon_index = int(np.argmin(np.abs(longitude-lon)))
            lat_index = int(np.argmin(np.abs(latitude-lat)))

        time_index = int(hour - 1)
        idx_rng = int(np.ceil(rng/div))

        try:
            temperature = np.array(dataset.variables['t'][time_index, :, lat_index-idx_rng:lat_index+idx_rng+1, lon_index-idx_rng:lon_index+idx_rng+1])
            x_wind = np.array(dataset.variables['u'][time_index, :, lat_index-idx_rng:lat_index+idx_rng+1, lon_index-idx_rng:lon_index+idx_rng+1])
            y_wind = np.array(dataset.variables['v'][time_index, :, lat_index-idx_rng:lat_index+idx_rng+1, lon_index-idx_rng:lon_index+idx_rng+1])
            geo = np.array(dataset.variables['z'][time_index, :, lat_index-idx_rng:lat_index+idx_rng+1, lon_index-idx_rng:lon_index+idx_rng+1])
        except IndexError:
            # if only the needed time was downloaded
            temperature = np.array(dataset.variables['t'][0, :, lat_index-idx_rng:lat_index+idx_rng+1, lon_index-idx_rng:lon_index+idx_rng+1])
            x_wind = np.array(dataset.variables['u'][0, :, lat_index-idx_rng:lat_index+idx_rng+1, lon_index-idx_rng:lon_index+idx_rng+1])
            y_wind = np.array(dataset.variables['v'][0, :, lat_index-idx_rng:lat_index+idx_rng+1, lon_index-idx_rng:lon_index+idx_rng+1])
            geo = np.array(dataset.variables['z'][0, :, lat_index-idx_rng:lat_index+idx_rng+1, lon_index-idx_rng:lon_index+idx_rng+1])

        lats = latitude[lat_index-idx_rng:lat_index+idx_rng+1]
        lons = longitude[lon_index-idx_rng:lon_index+idx_rng+1]

        return temperature, x_wind, y_wind, geo, lats, lons

    def generateECMWF(self, lat, lon, rng, time, file_name):
        temperature, x_wind, y_wind, geo, lats, lons = self.readNCDF(lat, lon, rng, time, file_name, 0.25)
        
        self.temperature = temperature
        self.x_wind = x_wind
        self.y_wind = y_wind
        self.geo = geo
        self.lats = lats
        self.lons = lons

    def spread(self, lat, lon, rng, time, file_name):
        temperature, x_wind, y_wind, geo, lats, lons = self.readNCDF(lat, lon, rng, time, file_name, 0.5)

        self.temperature_spr = temperature
        self.x_wind_spr = x_wind
        self.y_wind_spr = y_wind
        self.geo_spr = geo
        self.lats_spr = lats
        self.lons_spr = lons

    def perturb(self, t, u, v, z, t_spr, u_spr, v_spr):

        t_all = []
        u_all = []
        v_all = []

        # Using a different random number for each variable and height
        for level in range(len(t)):
            t_all.append(random.gauss(t[level], t_spr[level]))
            u_all.append(random.gauss(u[level], u_spr[level]))
            v_all.append(random.gauss(v[level], v_spr[level]))


        return np.array(t_all), np.array(u_all), np.array(v_all), z


    def generateProfile(self, t, u, v, z, h, lats, lons, spline=100, ref_time=None):

        level, speed, mags, dirs = self.convert(t, u, v, z)

        # # Temp Hack to get NRL data
        # speed = []

        # lat = lats[0]
        # lon = lons[0]

        # jd = date2JD(ref_time.year, ref_time.month, ref_time.day, ref_time.hour, ref_time.minute, ref_time.second)


        # for hh in level:

        #     t = getAtmDensity(lat, lon, hh, jd)

        #     speed.append(np.sqrt(consts.GAMMA*consts.R/consts.M_0*t)) 

        # speed = np.array(speed)
        # #######################


        # If the region extends above available data
        dh = 1000 #meters
        if h[0] > level[0]:

            missing_heights = np.arange(level[0] + dh, h[0] + dh, dh)

            lat = lats[0]
            lon = lons[0]

            jd = date2JD(ref_time.year, ref_time.month, ref_time.day, ref_time.hour, ref_time.minute, ref_time.second)


            for hh in missing_heights:

                t = getAtmDensity(lat, lon, hh, jd)
                u, v = getHWM(ref_time, lat, lon, hh/1000)
                mag = np.sqrt(u**2 + v**2)
                d = np.arctan2(u, v)

                speed = np.insert(speed, 0, np.sqrt(consts.GAMMA*consts.R/consts.M_0*t)) 
                mags = np.insert(mags, 0, mag)
                dirs = np.insert(dirs, 0, d)
                level = np.insert(level, 0, hh)


        level, speed, mags, dirs = self.spline(level, speed, mags, dirs, interp=spline)

        sounding = []
        for i in range(len(level)):
            sounding.append([level[i], speed[i], mags[i], dirs[i]])
        
        sounding = zInteg(h[0], h[1], np.array(sounding))

        
        return sounding

    def getProfile(self, lat, lon, heights, prefs, spline=100, ref_time=None):
        ''' Interpolates between starting and ending locations, converts, and splines 
            resulting curve into a smooth profile
            lat = [start lat, end lat]
        '''

        spread = True

        perturb = prefs.pert_num

        pts = self.interp(lat, lon, div=0.25)
        
        t = []
        u = []
        v = []
        z = []
        t_spr = []
        u_spr = []
        v_spr = []
        z_spr = []

        num_pts = len(pts)
        num_lvl = len(self.temperature)

        last_frac = 0

        range_lats = lat
        range_lons = lon

        for ii, pt in enumerate(pts):            

            lat = roundToNearest(pt[0], 0.25)
            lon = roundToNearest(pt[1], 0.25)%360

            try:
                lon_index = int(np.where(self.lons==lon)[0])
            except TypeError:
                if np.abs(lon - self.lons[-1]) < np.abs(lon - self.lons[0]):
                    lon_index = 0
                else:
                    lon_index = len(self.lons) - 1
            try:
                lat_index = int(np.where(self.lats==lat)[0])
            except TypeError:
                if np.abs(lat - self.lats[-1]) < np.abs(lat - self.lats[0]):
                    lat_index = 0
                else:
                    lat_index = len(self.lats) - 1
            

            frac = int(np.around((ii+1)/num_pts*num_lvl) + 1)

            t.append(self.temperature[last_frac:frac, lat_index, lon_index])
            u.append(self.x_wind[last_frac:frac, lat_index, lon_index])
            v.append(self.y_wind[last_frac:frac, lat_index, lon_index])
            z.append(self.geo[last_frac:frac, lat_index, lon_index])

            if perturb > 0 and prefs.pert_en:
                lat = roundToNearest(pt[0], 0.5)
                lon = roundToNearest(pt[1], 0.5)%360

                # I'm so sorry
                try:
                    try:
                        lon_index = int(np.where(self.lons_spr==lon)[0])
                    except TypeError:
                        if np.abs(lon - self.lons_spr[-1]) < np.abs(lon - self.lons_spr[0]):
                            lon_index = 0
                        else:
                            lon_index = len(self.lons_spr) - 1
                    try:
                        lat_index = int(np.where(self.lats_spr==lat)[0])
                    except TypeError:
                        if np.abs(lat - self.lats_spr[-1]) < np.abs(lat - self.lats_spr[0]):
                            lat_index = 0
                        else:
                            lat_index = len(self.lats_spr) - 1

                    t_spr.append(self.temperature_spr[last_frac:frac, lat_index, lon_index])
                    u_spr.append(self.x_wind_spr[last_frac:frac, lat_index, lon_index])
                    v_spr.append(self.y_wind_spr[last_frac:frac, lat_index, lon_index])
                    z_spr.append(self.geo_spr[last_frac:frac, lat_index, lon_index])

                except AttributeError:
                    spread = False
                    # no spread has been added
                    pass


            last_frac = frac

        t = np.array([item for sublist in t for item in sublist])
        u = np.array([item for sublist in u for item in sublist])
        v = np.array([item for sublist in v for item in sublist])
        z = np.array([item for sublist in z for item in sublist])

        t_spr = np.array([item for sublist in t_spr for item in sublist])
        u_spr = np.array([item for sublist in u_spr for item in sublist])
        v_spr = np.array([item for sublist in v_spr for item in sublist])
        z_spr = np.array([item for sublist in z_spr for item in sublist])

        sounding = self.generateProfile(t, u, v, z, heights, range_lats, range_lons, spline=spline, ref_time=ref_time)

        perturbations = []

        if perturb > 0 and prefs.pert_en and spread:
            
            for i in range(perturb):

                per_sounding = self.generateProfile(*self.perturb(t, u, v, z, t_spr, u_spr, v_spr), heights, range_lats, range_lons, spline=spline, ref_time=ref_time)

                perturbations.append(np.array(per_sounding))

        perturbations = np.array(perturbations)

        return sounding, perturbations

class Radio:
    def __init__(self):
        pass


if __name__ == '__main__':

    pass
