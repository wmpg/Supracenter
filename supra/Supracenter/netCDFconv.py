"""Stores atmospheric data in memory, finds weather profile closest to a given lat/lon"""
import sys

import numpy as np
from netCDF4 import Dataset
#from pyhdf.SD import SD, SDC

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

import supra.Supracenter.angleConv 
from supra.Supracenter.convLevels import convLevels
from supra.Supracenter.bisearch import bisearch


def readCustAtm(file_name, consts):
    
    """ Function for reading custom atmospheric data from txt files
    
    Example File: Elev Temp WSpd Wdir
                  0.0 12.10 2.06 0.0
                  etc ...

    Arguments: 
        file_name: [string] file to be converted
        consts: [object] physical constants

    Returns:
        data: [ndarray] Converted atmospheric profile
            [height (m), temperature (K), wind speed (m/s), wind direction (radians from north)]
    """

    with open(file_name) as f:
        next(f)
        data = np.array([0.0, 0.0, 0.0, 0.0])
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '')

            # Split the line by the delimiter
            line = line.split()
            
            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = float(entry.strip())

            # Transform Temperature to Speed of Sound (m/s)
            line[1] += consts.K
            line[1] = (consts.GAMMA*consts.R/consts.M_0*line[1])**0.5
            line[3] = line[3]*np.pi/180
            
            # Add the contents of the line to the data list
            data = np.vstack((data, line))

        # First row was all zeroes
        data = np.delete(data, 0, 0)

        return data

        
def findECMWFSound(lat, lon, dataset):
    """ Finds the atmospheric profile from a given lat/lon, from ECMWF archive data, to the closest grid point 

    Arguments:
        lat, lon: [float] latitude, longitude of the atmospheric profile needed
        grid_size: [float] the grid_size requested to ECMWF
        dataset: [ndarray] array of atmospheric profiles for each lat/lon grid point


    Returns:
        sounding: [ndarray] atmospheric profile to be used by ascan, in the form of 
                [height (m), speed of sound (m/s), wind speed (m/s), wind direction (radians from North due East)] 
    """

    ### Pull variables from the file
    # range of lat/lons used
    latitude = dataset[0]%360
    longitude = dataset[1]%360

    # Temperature, K, with height, lat, lon
    temperature = dataset[2]

    # Wind, positive from W to E, height, lat, lon
    x_wind = dataset[3]

    # Wind, positive from S to N, height, lat, lon
    y_wind = dataset[4]

    # pressure levels in hPa
    height = dataset[5]


    # Find the closest data point to the requested lat, lon
    try:
        lat_i = bisearch(latitude, lat%360)
        lon_i = bisearch(longitude, lon%360)
    except:
        print('ERROR: Cannot find given latitude/longitude in atmospheric profile range! (netCDFconv.py)')
        exit()


    # Create sounding array of all levels from the closest lat and lon
    # Initialize arrays
    sounding = np.array([0, 0, 0])
    row = np.array([0, 0, 0])

    # find total number of height layers
    LAYERS = len(temperature[:, 0])

    # build atmospheric profile
    for i in range(LAYERS):

        # Add each row of the atmospheric profile
        # Negative signs are because ECMWF defines direction from where wind is blowing TO, where cyscan defines it 
        # as where it is blowing FROM, (needed for wind angles) 
        row = np.array((float(temperature[i, lat_i, lon_i]), float(x_wind[i, lat_i, lon_i]), float(y_wind[i, lat_i, lon_i])))
        
        # Append row to sounding
        sounding = np.vstack((sounding, row))

    # First row was all zeroes
    sounding = np.delete(sounding, 0, 0)

    # Levels need to be converted to height. See https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels
    # Each level cooresponds to a height
    levels = np.array(convLevels())

    # Make array 2D so that they can be concatenated
    levels = np.expand_dims(levels, axis=1)

    # Add levels column to data columns
    sounding = np.hstack((levels, sounding))

    # Levels given by ECMWF are from the top down
    sounding = np.flip(sounding, axis=0)
    
    return sounding


def findMERRASound(lat, lon, dataset):
    """ Finds the atmospheric profile from a given lat/lon, from MERRA archive data, to the closest grid point 

    Arguments:
        lat, lon: [float] latitude, longitude of the atmospheric profile needed
        dataset: [ndarray] array of atmospheric profiles for each lat/lon grid point

    Returns:
        sounding: [ndarray] atmospheric profile to be used by ascan, in the form of 
                [height (m), speed of sound (m/s), wind speed (m/s), wind direction (radians from North due East)] 
    """

    ### Pull variables from the file
    # range of lat/lons used
    latitude = dataset[0]
    longitude = dataset[1]

    # Temperature, K, with height, lat/lon
    temperature = dataset[2]

    # Wind, positive from W to E
    x_wind = dataset[3]

    # Wind, positive from S to N
    y_wind = dataset[4]

    # height of level
    height = dataset[5]

    # Find the closest data point to the requested lat, lon
    lat_i = bisearch(latitude, lat)
    lon_i = bisearch(longitude, lon)

    # Create sounding array of all levels from the closest lat and lon
    # Initialize arrays
    sounding = np.array([0, 0, 0, 0])
    row      = np.array([0, 0, 0, 0])

    # find total number of atmospheric layers
    LAYERS = len(height[0])

    # build atmospheric profile
    for i in range(LAYERS):

        # Add each row of the atmospheric profile
        # Negative signs are because ECMWF defines direction from where wind is blowing TO, where ascan defines it 
        # as where it is blowing FROM, (needed for wind angles) 
        row = np.array([float(height[0, i, lat_i, lon_i]), float(temperature[0, i, lat_i, lon_i]), \
                        float(x_wind[0, i, lat_i, lon_i]), float(y_wind[0, i, lat_i, lon_i])])
        
        # Append row to sounding
        sounding = np.vstack((sounding,row))

    # First row was all zeroes
    sounding = np.delete(sounding, 0, 0)

    # Levels given by MERRA are from the top down
    sounding = np.flip(sounding, axis=0)

    return sounding


def findUKMOSound(lat, lon, dataset):
    """ Finds the atmospheric profile from a given lat/lon, from UKMO archive data, to the closest grid point 

    Arguments:
        lat, lon: [float] latitude, longitude of the atmospheric profile needed
        dataset: [ndarray] array of atmospheric profiles for each lat/lon grid point

    Returns:
        sounding: [ndarray] atmospheric profile to be used by ascan, in the form of 
                [height (m), speed of sound (m/s), wind speed (m/s), wind direction (radians from North due East)] 
    """

    ### Pull variables from the file
    # range of lat/lons used
    latitude = dataset[0]
    longitude = dataset[1]

    # Temperature, K, with height, lat/lon
    temperature = dataset[2]

    # Wind, positive from W to E
    x_wind = dataset[3]

    # Wind, positive from S to N
    y_wind = dataset[4]

    # height of level
    height = dataset[5]

    # Find the closest data point to the requested lat, lon
    lat_i = bisearch(latitude, lat)
    lon_i = bisearch(longitude, lon)

    # Create sounding array of all levels from the closest lat and lon
    # Initialize arrays
    sounding = np.array([0, 0, 0, 0])
    row = np.array([0, 0, 0, 0])

    # find total number of atmospheric layers
    LAYERS = len(height[0])

    # build atmospheric profile
    for i in range(0, LAYERS):

        # Add each row of the atmospheric profile
        row = np.array([float(height[0, i, lat_i, lon_i]), float(temperature[0, i, lat_i, lon_i]), \
                        float(x_wind[0, i, lat_i, lon_i]), float(y_wind[0, i, lat_i, lon_i])])
        
        # Append row to sounding
        sounding = np.vstack((sounding,row))

    # First row was all zeroes
    sounding = np.delete(sounding, 0, 0)

    return sounding

def storeNetCDFECMWF(file_name, lat, lon, consts, start_time=0):
    """ HELPER FUNCTION: Reads ECMWF netCDF file and stores it in memory for reuse later in the program 

    Arguments:
        file_name: [String] name of file
        consts: [object] physical constants

    Returns:
        store_data: [list[ndarray]] atmospheric profile for the search area
    """

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
    #print(dataset.variables)

    lon_index = int(np.around((lon%360)*4))
    lat_index = int(np.around(-(lat+90)*4))# - 90*4

    longitude = np.array(dataset.variables['longitude'][lon_index-20:lon_index+21])
    latitude = np.array(dataset.variables['latitude'][lat_index-20:lat_index+21])

    level = np.array(dataset.variables['level'])
    #pressure 1 - 1000 hPa , non-linear
    
    time = np.array(dataset.variables['time'])
    #not known

    start_time = int(start_time)
    
    # time, (number), level, lat, lon
    temperature = np.array(dataset.variables['t'][start_time, :, lat_index-20:lat_index+21, lon_index-20:lon_index+21])
    x_wind = -np.array(dataset.variables['u'][start_time, :, lat_index-20:lat_index+21, lon_index-20:lon_index+21])
    y_wind = -np.array(dataset.variables['v'][start_time, :, lat_index-20:lat_index+21, lon_index-20:lon_index+21])

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
    
    level = convLevels()
    level = np.flipud(np.array(level))

    # Store data in a list of arrays
    store_data = [np.array(latitude[:]), np.array(longitude[:]), np.array(temps), np.array(mags), np.array(dirs), np.array(level)]
    
    # Store data in a list of ndarrays
    # store_data = [np.array(latitude[:]), np.array(longitude[:]), np.array(temperature[:]), np.array(x_wind[:]), np.array(y_wind[:])]
    dataset.close()

    return store_data

    def parseGeneralECMWF(file_name, lat, lon, time, variables):

        try:
        # Read the file
            dataset = Dataset(file_name, "r+", format="NETCDF4")
        except:
            print("ERROR: Unable to read weather file: ", file_name)
            exit()

        lon_index = int(np.around((lon%360)*4))
        lat_index = int(np.around(-(lat+90)*4))
        time_index = int(time)

        sounding = []

        for var in variables:
            if (set([var]) < set(dataset.variables.keys())):
                sounding.append(np.array(dataset.variables[var][time_index, :, lat_index, lon_index]))

        sounding = np.array(sounding)
        return sounding

def storeHDF(file_name, consts):
    """ HELPER FUNCTION: Reads MERRA HDF file and stores it in memory for reuse later in the program

    AS OF AUGUST 2018 this has been changed to .netCDF files due to incompatibilities with obspy 

    Arguments:
        file_name: [String] name of file 
        consts: [object] physical constants

    Returns:
        store_data: [list[ndarray]] atmospheric profile for the search area
    """

    print('Converting weather data. This may take a while...')

    # # # Open file.
    # hdf = SD(file_name, SDC.READ)
   
    # # # Read 
    # dataset = hdf.datasets()
    try:
        dataset = Dataset(file_name, "r+", format="NETCDF4")
    except:
        print("ERROR: Unable to read weather file: ", file_name)
        exit()


    # # Check file type
    # if not (set(['H', 'T', 'U', 'V', 'XDim:EOSGRID', 'YDim:EOSGRID']) < set(dataset.keys())):
    #     print('FILE ERROR: File does not contain the correct parameters! Variables Required: H, T, U, V, XDim:EOSGRID, YDim:EOSGRID')
    #     sys.exit()

    # Read dataset.
    # [time, height, y, x]

    height = np.array(dataset.variables['H'][:, :, :, :])
    temperature = np.array(dataset.variables['T'][:, :, :, :])
    x_wind = np.array(dataset.variables['U'][:, :, :, :])
    y_wind = np.array(dataset.variables['V'][:, :, :, :])

    longitude = np.array(dataset.variables['lon'][:])
    latitude = np.array(dataset.variables['lat'][:])

    # Conversions
    temps = (consts.GAMMA*consts.R/consts.M_0*temperature[:])**0.5

    # Magnitude of winds (m/s)
    mags = np.sqrt(x_wind**2 + y_wind**2)

    # Direction the winds are coming from, angle in radians from North due East
    dirs = (np.arctan2(-y_wind, -x_wind))%(2*np.pi)*180/np.pi
    dirs = supra.Supracenter.angleConv.angle2NDE(dirs)*np.pi/180

    # Store data in a list of arrays
    store_data = [latitude, longitude, temps, mags, dirs, height]

    return store_data


def storeNetCDFUKMO(file_name, area, consts):
    """ HELPER FUNCTION: reads converted UKMO .nc and stores it in memory for reuse later in the program 

    Arguments:
        file_name: [String] name of file
        area: [list] search area [lat min, lat max, lon min, lon max]
        consts: [object] physical constants

    Returns:
        store_data: [list[ndarray]] atmospheric profile for the search area
    """

    print('Converting weather data. Do not terminate while converting! This may take a while...')

    # Only take search area weather, otherwise memory errors
    lat_min = area[0]%360
    lat_max = area[1]%360
    lon_min = area[2]%360
    lon_max = area[3]%360

    # Read the file
    try:
        dataset = Dataset(file_name, "r+", format="NETCDF4")
    except:
        print("ERROR: Unable to read weather file: ", file_name)
        exit()


    # Check file type
    if not (set(['ht', 'u', 'v', 'latitude', 'longitude', 'temp']) < set(dataset.variables.keys())):
        print('FILE ERROR: File does not contain the correct parameters! Variables Required: ht, u, v, latitude, longitude, temp')
        sys.exit()

    # Pull variables from the file

    # lat/lon for winds
    latitude = dataset.variables['latitude'][:]
    longitude = dataset.variables['longitude'][:]

    # indicies of lat/lon search area
    lat_min_i = bisearch(latitude, lat_min)
    lat_max_i = bisearch(latitude, lat_max) + 1
    lon_min_i = bisearch(longitude, lon_min)
    lon_max_i = bisearch(longitude, lon_max) + 1

    # only use lat/lons within search area
    latitude = dataset.variables['latitude'][lat_min_i:lat_max_i]
    longitude = dataset.variables['longitude'][lon_min_i:lon_max_i]

    # temperature = dataset.variables['temp'][:, :, lat_min_i_1:lat_max_i_1, lon_min_i_1:lon_max_i_1]
    temperature = dataset.variables['temp'][:, :, lat_min_i:lat_max_i, lon_min_i:lon_max_i]

    # Wind, positive from W to E
    x_wind = dataset.variables['u'][:, :, lat_min_i:lat_max_i, lon_min_i:lon_max_i]

    # Wind, positive from S to N
    y_wind = dataset.variables['v'][:, :, lat_min_i:lat_max_i, lon_min_i:lon_max_i]

    # height
    ht = dataset.variables['ht'][:, :, lat_min_i:lat_max_i, lon_min_i:lon_max_i]

    # Transform Temperature to Speed of Sound (m/s)
    temps = (consts.GAMMA*consts.R/consts.M_0*temperature)**0.5

    # Magnitude of winds (m/s)
    mags = np.sqrt(x_wind**2 + y_wind**2)

    # Direction the winds are coming from, angle in radians from North due East
    dirs = (np.arctan2(-y_wind, -x_wind))%(2*np.pi)*180/np.pi
    dirs = supra.Supracenter.angleConv.angle2NDE(dirs)*np.pi/180

    # convert heights from geopotential to geometric
    height = supra.Supracenter.angleConv.geopot2Geomet(ht)

    # data is needed in reverse order
    temps = np.flipud(np.array(temps))
    mags = np.flipud(np.array(mags))
    dirs = np.flipud(np.array(dirs))
    height = np.flipud(np.array(height))

    # Store data in a list of arrays
    store_data = [np.array(latitude), np.array(longitude), temps, mags, dirs, height]

    dataset.close()

    return store_data

