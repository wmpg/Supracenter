""" Given the UKMO .nc file, latitude, longitude, get the atmospheric weather profile. """


import os
import argparse
import datetime

import numpy as np
import scipy.interpolate
from netCDF4 import Dataset


import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

import wmpl
from supra.Supracenter.bisearch import bisearch


import matplotlib.pyplot as plt




def saveWeatherProfile(file_path, pressures, heights, temperatures, wind_dirs, wind_speeds):
    """ Save the weather profile into a file. """

    # Sort by descending height
    ht_rev_indices = heights.argsort()[::-1]
    pressures = pressures[ht_rev_indices]
    heights = heights[ht_rev_indices]
    temperatures = temperatures[ht_rev_indices]
    wind_dirs = wind_dirs[ht_rev_indices]
    wind_speeds = wind_speeds[ht_rev_indices]

    with open(file_path, 'w') as f:

        f.write('# Pres(hPa),  Ht (km), Temp (C), Wind dir, Wind v (m/s)\n')

        for p, h, t, wd, ws in zip(pressures, heights, temperatures, wind_dirs, wind_speeds):
            f.write(' {:10.1f} {:9.3f} {:9.1f} {:9.1f} {:9.1f}\n'.format(p, h/1000, t - 273.15, np.degrees(wd), ws))



def extractUKMOData(file_path, lat, lon):
    """ Get the weather profile from the UKMO .nc file. """

    # Load the UKMO file
    dataset = Dataset(file_path, "r+", format="NETCDF4")

    # Get forecast time
    forecast_dt = datetime.datetime(1970, 1, 1, 0, 0, 0) \
        + datetime.timedelta(hours=float(dataset.variables['forecast_reference_time'][0]))
    
    # Load lat/lons
    latitudes = dataset.variables['latitude'][:]
    longitudes = dataset.variables['longitude'][:]


    # Indicies of closest lat/lon in the data
    lat_i = bisearch(latitudes, lat)
    lon_i = bisearch(longitudes, lon)

    # Get the pressure in hPa - this is the independant variable
    pressures = dataset.variables['pressure'][:]

    # Extract geopetential heights
    heights = dataset.variables['geopotential_height'][:, lat_i, lon_i]

    # Extract temperatures
    temperatures = dataset.variables['air_temperature'][:, lat_i, lon_i]

    # Wind, positive from W to E
    x_wind = dataset.variables['x_wind'][:, lat_i, lon_i]

    # Wind, positive from S to N
    y_wind = dataset.variables['y_wind'][:, lat_i, lon_i]

    return forecast_dt, pressures, heights, temperatures, x_wind, y_wind



def getUKMOWeatherProfile(file_path, lat, lon, file_path2=None, dt=None, height_min=None, height_max=None, 
    height_step=None):
    """ 
    Keyword arguments:
        height_min: [float] Minimum heiight in km.
        height_max: [float] Maximum heiight in km.
        height_step: [float] Height step in km, 0.1 km by default.
    """

    # Extract the direcotyr where the file is
    dir_path = os.path.dirname(file_path)

    # Read teh data from the file for the given geo coordinates
    forecast_dt, pressures, heights, temperatures, x_wind, y_wind = extractUKMOData(file_path, lat, lon)

    # If the second file was given, load it and interpolate the weather profile for the given time
    if (file_path2 is not None) and (dt is not None):
        print('Interpolating between 2 files...')

        # Load the data from the second file
        forecast_dt2, pressures2, heights2, temperatures2, x_wind2, y_wind2 = extractUKMOData(file_path2, \
            lat, lon)

        # Make sure the forecasts are ordered in time
        if forecast_dt > forecast_dt2:

            # Swap variables to order them in time
            forecast_dt2, pressures2, heights2, temperatures2, x_wind2, y_wind2, forecast_dt, pressures, \
            heights, temperatures, x_wind, y_wind = forecast_dt, pressures, heights, temperatures, x_wind, \
            y_wind, forecast_dt2, pressures2, heights2, temperatures2, x_wind2, y_wind2


        print('Reference time of the 1st file:', forecast_dt)
        print('Reference time of the 2nd file:', forecast_dt2)
        print('Event time:', dt)

        ### Compute the weight for computation by assuming a linear dependence on time ###

        # Compute the total time difference between the forecasts in seconds
        time_diff = (forecast_dt2 - forecast_dt).total_seconds()

        # Compute the weight for every forecast
        weight1 = (dt - forecast_dt).total_seconds()/time_diff
        weight2 = (forecast_dt2 - dt).total_seconds()/time_diff

        ###

        # Compute the interpolated values by time (assume fixed steps in pressure)
        heights = weight1*heights + weight2*heights2
        temperatures = weight1*temperatures + weight2*temperatures2
        x_wind = weight1*x_wind + weight2*x_wind2
        y_wind = weight1*y_wind + weight2*y_wind2


        # Construct a file name for the plots
        plot_name = dt.strftime('%Y%m%d_%H%M%S')

    else:
        print('Using a single .nc file...')

        print('Forecant time:', forecast_dt)

        # Construct a file name for the plots
        plot_name = forecast_dt.strftime('%Y%m%d_%H%M%S')


    # Interpolate all values by height (heights are negative because the X for interpolation has to be 
    #   increasing)
    pressures_interpol = scipy.interpolate.CubicSpline(-heights, pressures)
    temperatures_interpol = scipy.interpolate.CubicSpline(-heights, temperatures)
    x_wind_interpol = scipy.interpolate.CubicSpline(-heights, x_wind)
    y_wind_interpol = scipy.interpolate.CubicSpline(-heights, y_wind)


    # Construct the height array for interpolation
    if height_min is None:
        height_min = np.min(heights)
    else:
        height_min *= 1000

    if height_max is None:
        height_max = np.max(heights)
    else:
        height_max *= 1000

    if height_step is None:
        height_step = 100 # m
    else:
        height_step *= 1000

    heights_arr = np.arange(height_min, height_max, height_step)


    # Extract interpolated values
    pressures_arr = pressures_interpol(-heights_arr)
    temperatures_arr = temperatures_interpol(-heights_arr)
    x_wind_arr = x_wind_interpol(-heights_arr)
    y_wind_arr = y_wind_interpol(-heights_arr)

    # Magnitude of winds (m/s)
    wind_mags = np.sqrt(x_wind_arr**2 + y_wind_arr**2)

    # Direction the winds are coming from, angle in radians from North due East
    wind_dirs = (np.arctan2(-y_wind_arr, -x_wind_arr))%(2*np.pi)*180/np.pi
    wind_dirs = supra.Supracenter.angleConv.angle2NDE(wind_dirs)*np.pi/180



    # Plot the temperatures
    plt.plot(temperatures_arr, heights_arr/1000)
    plt.xlabel('Temperature (K)')
    plt.ylabel('Height (km)')
    plt.savefig(os.path.join(dir_path, plot_name) + "_temperature.png", dpi=300)
    plt.show()


    # Plot the pressure
    plt.plot(pressures_arr, heights_arr/1000)
    plt.xlabel('Pressure (hPa)')
    plt.ylabel('Height (km)')
    plt.savefig(os.path.join(dir_path, plot_name) + "_pressure.png", dpi=300)
    plt.show()    


    # Plot the wind speed and direction by heights
    plt.scatter(wind_mags, heights_arr/1000, c=np.degrees(wind_dirs), s=1)
    plt.xlabel('Wind speed (m/s)')
    plt.ylabel('Height (km)')
    plt.colorbar(label='Wind direction +E of due N (deg)')
    plt.savefig(os.path.join(dir_path, plot_name) + "_wind.png", dpi=300)
    plt.show()

    saveWeatherProfile(os.path.join(dir_path, plot_name) + "_weather.txt", pressures_arr, heights_arr, temperatures_arr, wind_dirs, \
        wind_mags)


if __name__ == "__main__":


    ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description="Given the UKMO .nc file, latitude, longitude, get the atmospheric weather profile. Optionally, 2 .nc files can be given, in which case the time for linear interpolation can be given with the -t option.")

    arg_parser.add_argument('nc_file_path', nargs='+', metavar='NC_PATH', type=str, \
        help='Path to the UKMO .nc file (or 2 files can be given).')

    arg_parser.add_argument('lat', nargs=1, metavar='LAT', type=float, \
        help='Latitude +N in degrees.')

    arg_parser.add_argument('lon', nargs=1, metavar='LON', type=float, \
        help='Longitude +E in degrees.')

    arg_parser.add_argument('-b', '--bottom', type=float, help='Bottom height (km)')
    arg_parser.add_argument('-u', '--up', type=float, help='Upper height (km)')
    arg_parser.add_argument('-t', '--time', type=str, help='Time for the wind profile in the YYYYMMDD-hhmmss format.')
    arg_parser.add_argument('-s', '--step', type=float, help='Height step (km)')


    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    #########################


    # Check if the second file was given
    if len(cml_args.nc_file_path) > 1:
        file_path2 = cml_args.nc_file_path[1]
    else:
        file_path2 = None


    # If the time was given, parse it
    if cml_args.time is not None:
        dt = datetime.datetime.strptime(cml_args.time, "%Y%m%d-%H%M%S")
    else:
        dt = None


    getUKMOWeatherProfile(cml_args.nc_file_path[0], cml_args.lat[0], cml_args.lon[0], file_path2=file_path2, \
        dt=dt, height_min=cml_args.bottom, height_max=cml_args.up, height_step=cml_args.step)

