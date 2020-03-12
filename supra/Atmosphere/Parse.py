import sys
import numpy as np

from supra.Utils.Classes import Constants
from supra.Supracenter.netCDFconv import storeHDF, storeNetCDFECMWF, storeNetCDFUKMO, readCustAtm, storeAus
from supra.Atmosphere.parseRadio import parseRadio

def parseWeather(setup, t=0, lat=None, lon=None):
    """ Generates a cubic weather profile of all altitudes within a 5 deg lat/lon area from lat/lon centre. 
    """


    consts = Constants()
     # Parse weather type
    setup.weather_type = setup.weather_type.lower()
    if setup.weather_type not in ['custom', 'none', 'ecmwf', 'merra', 'ukmo', 'binary', 'radio']:
        print('INI ERROR: [Atmospheric] weather_type must be one of: custom, radio, none, ecmwf, merra, ukmo, binary')
        return None

    # Custom weather data    
    if setup.weather_type == 'custom':

        # try:
        sounding = readCustAtm(setup.sounding_file, consts)
        # except:
        #     print('ERROR: Could not read custom atmosphere profile!')
        #     exit()

        # Check file name
        if len(sounding) < 2:
            print('FILE ERROR: file must have at least two data points!')
            sys.exit()


    # MERRA-2
    elif setup.weather_type == 'merra':

        # Check file type
        if '.nc' not in setup.sounding_file:
            print("FILE ERROR: custom data set must be a .nc file!")
            sys.exit()

        # # Run data fetching script
        # if setup.get_data == True:
        #     print('Getting data...')
        #     fetchMERRA(setup)

        try:
            sounding = storeHDF(setup.sounding_file, consts)
        except:
            print('ERROR: could not read merra atmosphere profile!')
            exit()

    # ECMWF
    elif setup.weather_type == 'ecmwf':

        # Check file type
        if '.nc' not in setup.sounding_file:
            print("FILE ERROR: custom data set must be a .nc file from ECMWF!")
            return None

        # # Run data fetching script
        # if setup.get_data == True:
        #     fetchECMWF(setup, setup.sounding_file)

        # GetIRISData/MakeIRISPicks

        try:

            #Get closest hour
            start_time = (setup.fireball_datetime.hour + np.round(setup.fireball_datetime.minute/60) + t)%24
            if lat is not None and lon is not None:
                sounding = storeNetCDFECMWF(setup.sounding_file, lat, lon, consts, start_time=start_time)

            else:
                sounding = storeNetCDFECMWF(setup.sounding_file, setup.lat_centre, setup.lon_centre, consts, start_time=start_time)

        # SeismicTrajectory
        except:

            try:
                start_time = (setup.fireball_datetime.hour + t)%24
                sounding = storeNetCDFECMWF(setup.sounding_file, setup.x0, setup.y0, consts, start_time=start_time)
            except:
                print("ERROR: Unable to use weather file, or setup.start_datetime/setup.atm_hour is set up incorrectly. Try checking if sounding_file exists")

        

    # UKMO
    elif setup.weather_type == 'ukmo':
        # Check file type
        if '.nc' not in setup.sounding_file:
            print("FILE ERROR: custom data set must be a .nc file from UKMO!")
            sys.exit()
        
        try:
            sounding = storeNetCDFUKMO(setup.sounding_file, setup.search_area, consts)
        except:
            print('ERROR: Could not read UKMO atmosphere profile!')
            exit()

    elif setup.weather_type == 'binary':

        
        sounding = storeAus(consts)

    elif setup.weather_type == 'radio':

        time_of_sound = [setup.fireball_datetime.year, setup.fireball_datetime.month, setup.fireball_datetime.day, setup.fireball_datetime.hour]
        sounding = parseRadio(setup.sounding_file, time_of_sound)

    else:

        # Sample fake weather profile, the absolute minimum that can be passed
        sounding = np.array([[    0.0, setup.v_sound, 0.0, 0.0],
                             [    0.0, setup.v_sound, 0.0, 0.0],
                             [99999.0, setup.v_sound, 0.0, 0.0]])

    return sounding