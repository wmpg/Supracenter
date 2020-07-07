import random
import numpy as np
from netCDF4 import Dataset

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from scipy.interpolate import CubicSpline

from supra.Utils.AngleConv import roundToNearest
from supra.Utils.Classes import Constants
from supra.Supracenter.cyzInteg import zInteg

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

        lon_index = int(np.where(longitude==lon)[0])
        lat_index = int(np.where(latitude==lat)[0])
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


    def generateProfile(self, t, u, v, z, h, spline=100):
        level, speed, mags, dirs = self.spline(*self.convert(t, u, v, z), interp=spline)

        sounding = []
        for i in range(len(level)):
            sounding.append([level[i], speed[i], mags[i], dirs[i]])
        
        sounding = zInteg(h[0], h[1], np.array(sounding))
        
        return sounding

    def getProfile(self, lat, lon, heights, prefs, spline=100):
        ''' Interpolates between starting and ending locations, converts, and splines 
            resulting curve into a smooth profile
            lat = [start lat, end lat]
        '''

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

            last_frac = frac

        t = np.array([item for sublist in t for item in sublist])
        u = np.array([item for sublist in u for item in sublist])
        v = np.array([item for sublist in v for item in sublist])
        z = np.array([item for sublist in z for item in sublist])

        t_spr = np.array([item for sublist in t_spr for item in sublist])
        u_spr = np.array([item for sublist in u_spr for item in sublist])
        v_spr = np.array([item for sublist in v_spr for item in sublist])
        z_spr = np.array([item for sublist in z_spr for item in sublist])

        sounding = self.generateProfile(t, u, v, z, heights, spline=spline)

        perturbations = []

        if perturb > 0 and prefs.pert_en:
            
            for i in range(perturb):

                per_sounding = self.generateProfile(*self.perturb(t, u, v, z, t_spr, u_spr, v_spr), heights, spline=spline)

                perturbations.append(np.array(per_sounding))

        perturbations = np.array(perturbations)

        return sounding, perturbations

class Radio:
    def __init__(self):
        pass