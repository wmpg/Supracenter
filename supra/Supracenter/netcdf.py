import sys

import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
#from pyhdf.SD import SD, SDC

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

import supra.Supracenter.angleConv 
from supra.Supracenter.convLevels import convLevels
from supra.Supracenter.bisearch import bisearch

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


def storeNetCDFECMWF(file_name, lat, lon, consts, start_time=0):
    """ HELPER FUNCTION: Reads ECMWF netCDF file and stores it in memory for reuse later in the program 

    Arguments:
        file_name: [String] name of file
        consts: [object] physical constants

    Returns:
        store_data: [list[ndarray]] atmospheric profile for the search area
    """

    print('Converting weather data. This may take a while...')
    print('AAAAAAAAAA')
    try:
        # Read the file
        dataset = Dataset(file_name, "r+", format="NETCDF4")
    except:
        print("ERROR: Unable to read weather file: ", file_name)
        exit()

    # Check file type
    if not (set(['t', 'u', 'v', 'latitude', 'longitude']) < set(dataset.variables.keys())):
        print('FILE ERROR: File does not contain the correct parameters! Variables Required: t, u, v, latitude, longitude')
        sys.exit()

    # Check against UKMO files
    elif not (set(['level']) < set(dataset.variables.keys())):
        a = input('WARNING: File is an unrecognized .nc! May be old ECMWF or UKMO file. Running may cause program to crash.\
                                                                                                 Bypass? (y/n) ')
        if a.lower() != 'y': 
            sys.exit()

    # print(file_name)
    # print(dataset)
    # print(dataset.variables)
    # dataset.close()
    # exit()

    lon_index = int((lon%360)*4)
    lat_index = int(-(lat+90)*4) -1# - 90*4
    longitude = np.array(dataset.variables['longitude'][:])
    latitude = np.array(dataset.variables['latitude'][:])

    level = np.array(dataset.variables['level'])
    #pressure 1 - 1000 hPa , non-linear
    
    time = np.array(dataset.variables['time'])
    #not known

    start_time = int(start_time)

    # time, (number), level, lat, lon
    #z =  np.array(dataset.variables['z'][start_time, :, lat_index, lon_index])
    temperature = np.array(dataset.variables['t'][start_time, :, lat_index, lon_index])
    x_wind = np.array(dataset.variables['u'][start_time, :, lat_index, lon_index])
    y_wind = np.array(dataset.variables['v'][start_time, :, lat_index, lon_index])

    # Transform Temperature to Speed of Sound (m/s)
    temps = (consts.GAMMA*consts.R/consts.M_0*temperature[:])**0.5

    # Magnitude of winds (m/s)
    mags = np.sqrt(x_wind[:]**2 + y_wind[:]**2)

    #ECWMF Defines winds as +u -> west to east      ^ +v
    #                       +v -> south to north    |
    #                                               o--> +u
    #https://confluence.ecmwf.int/pages/viewpage.action?pageId=111155337
    # Direction the winds are coming from, angle in radians from North due East

    dirs = (np.arctan2(y_wind[:], x_wind[:]))%(2*np.pi)
    dirs = supra.Supracenter.angleConv.angle2NDE(np.degrees(dirs))

    h = convLevels()
    # h = np.flipud(np.array(h))

    # h = np.flipud(np.array(h))
    #level = np.flipud(np.array(level))
    # Store data in a list of arrays
    store_data = [lat, lon, np.array(temps), np.array(mags), np.array(dirs), np.array(h), np.array(level)]
    
    # Store data in a list of ndarrays
    # store_data = [np.array(latitude[:]), np.array(longitude[:]), np.array(temperature[:]), np.array(x_wind[:]), np.array(y_wind[:])]
    dataset.close()

    with open('/home/luke/Desktop/output_profile.txt', 'w') as f:
        f.write('Pressure(hPa), Geopotential Height(m), Temperature(K), u_comp_wind(m/s), v_comp_wind(m/s)\n')
        for i in range(len(level)):
            f.write('{:} , {:} , {:} , {:} , {:} \n'.format(level[i], h[i] ,temps[i], x_wind[i], y_wind[i]))

    plt.plot(temps, h)
    plt.show()
    plt.plot(mags, h)
    plt.show()
    plt.plot(dirs, h)
    plt.show()

if __name__ == '__main__':
    file_name = '/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/StubenbergReanalysis.nc'
    lat = 48.3314
    lon = 13.0706
    consts = Constants()

    storeNetCDFECMWF(file_name, lat, lon, consts)
