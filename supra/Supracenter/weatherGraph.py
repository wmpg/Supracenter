import numpy as np
import matplotlib.pyplot as plt

import sys
import time

from supra.Supracenter.netCDFconv import storeNetCDFECMWF, storeHDF, storeNetCDFUKMO, readCustAtm
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.angleConv import angle2NDE
from supra.Supracenter.SPPT import perturb
from supra.Fireballs.SeismicTrajectory import Constants
from wmpl.Utils.PlotMap import GroundMap

def graphWeather(file_name, weather_type, S, D, method, number, spread_file='', setup=0):
    
    consts = Constants()

    ### Get atmospheric data ###

    # Isotropic weather profile
    if weather_type == 'none':
        dataset = np.array([[0, 330, 0, 0], [1000, 330, 0, 0]])

    # Custom weather data     
    elif weather_type == 'custom':
        dataset = readCustAtm(file_name, consts)

        # Check file name
        if len(dataset) < 2:
            print('FILE ERROR: file must have at least two data points!')
            sys.exit()

        # Check file type
        if '.txt' not in file_name:
            print("FILE ERROR: custom data set must be a .txt file!")
            sys.exit()


    # MERRA-2
    elif weather_type == 'merra':

        # Check file type
        if '.hdf' not in file_name:
            print("FILE ERROR: custom data set must be a .hdf file!")
            sys.exit()

        dataset = storeHDF(file_name, consts)

    # ECMWF
    elif weather_type == 'ecmwf':

        # Check file type
        if '.nc' not in file_name:
            print("FILE ERROR: custom data set must be a .nc file from ECMWF!")
            sys.exit()

        dataset = storeNetCDFECMWF(file_name, D[0], D[1], consts)

    # UKMO
    elif weather_type == 'ukmo':

        # Check file type
        if '.nc' not in file_name:
            print("FILE ERROR: custom data set must be a .nc file from UKMO!")
            sys.exit()

        dataset = storeNetCDFUKMO(file_name, search_area, consts)

    else:
        print("Unrecognized weather type!")
        exit()

    print("Status: Atmospheric Data Loaded")

    
    if weather_type is not 'none':
        ref_pos = 0
        merged_data, _ = getWeather(S, D, weather_type, ref_pos, dataset, convert=False)
    else:
        return None

    h = merged_data[:, 0]
    T = merged_data[:, 1]
    mags = merged_data[:, 2]
    dirs = merged_data[:, 3]

    #convert speed of sound to temp
    t = np.square(T)*consts.M_0/consts.GAMMA/consts.R
    dirs = np.degrees(dirs)
    #dirs = angle2NDE(dirs)
      
    fig = plt.figure(figsize=plt.figaspect(0.5))
    fig.set_size_inches(20.9, 11.7)  

    ax1 = fig.add_subplot(221)
    ax1.set_title('Temperature Profile')
    ax1.set_xlabel('Temperature (K)')
    ax1.set_ylabel('Height (m)')
    ax1.plot(t, h, c='k') 

    ax2 = fig.add_subplot(222)
    ax2.set_title('Wind Magnitude Profile')
    ax2.set_xlabel('Wind Speed (m/s)')
    ax2.set_ylabel('Height (m)')
    ax2.plot(mags, h, c='k')

    ax3 = fig.add_subplot(223)
    ax3.set_title('Wind Direction Profile')
    ax3.set_xlabel('Wind Direction (deg from north)')
    ax3.set_ylabel('Height (m)')
    ax3.plot(dirs, h, c='k')

    ax4 = fig.add_subplot(224)
    ax4 = GroundMap([np.radians(S[0]), np.radians(D[0])], [np.radians(S[1]), np.radians(D[1])], border_size=20, color_scheme='light')
    ax4.scatter(np.radians(S[0]), np.radians(S[1]), s=10, marker='o', c='g')
    ax4.scatter(np.radians(D[0]), np.radians(D[1]), s=10, marker='x', c='r')

    plt.title('Weather Profile')

    sounding = [D[0], D[1], merged_data[:, 1], merged_data[:, 2], merged_data[:, 3], merged_data[:, 0]]

    if method != 'none' and spread_file != '':
        #Perturbations
        for i in range(number):
            sys.stdout.write("\rPerturbation {:}    ".format(i+1))
            sys.stdout.flush()
            time.sleep(0.001)
            sounding_p = perturb(setup, sounding, method, spread_file=spread_file, lat=D[0], lon=D[1], line=True)

            T = sounding_p[2]
            mags = sounding_p[3]
            dirs = sounding_p[4]

            #convert speed of sound to temp
            t = np.square(T)*consts.M_0/consts.GAMMA/consts.R
            #dirs = np.degrees(dirs)
            #dirs = angle2NDE(dirs)


            ax1.plot(t, h, c='g', alpha=0.2)
            ax2.plot(mags, h, c='g', alpha=0.2) 
            ax3.plot(dirs, h, c='g', alpha=0.2) 

    plt.show()

if __name__ == '__main__':
    file_name = '/home/luke/Desktop/StubenbergReanalysis.nc'
    S = np.array([42.0, 42.0, 50000])
    D = np.array([42.0, 41.0,  1000])
    method = 'none'
    number = 100
    weather_type = 'ecmwf'
    graphWeather(file_name, weather_type, S, D, method, number)