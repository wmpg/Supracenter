import os

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

from supra.Fireballs.SeismicTrajectory import Constants
from supra.Supracenter.angleConv import angle2NDE, roundToNearest

def rGen(mode='normal'):
    """ 
    Generates a random number, r, based off of a given distribution

    Arguments:
    mode: [String] distribution used to generate r
            - normal: a normal distribution around 0, stDev = 0.5
            - uniform: a uniform distribution between -0.5 and 0.5
            - spread: a normal distribution around 0, stDev = 1.0

    Returns:
    r: [float] the randomly generated float
    """

    if mode == 'normal':
        r = np.random.normal(0, 0.5)
    elif mode == 'uniform':
        r = np.random.uniform(-0.5, 0.5)
    elif mode == 'spread':
        # r = np.random.normal()
        r = np.random.uniform(-1, 1)
    else:
        print("ERROR [rGen]: Unknown distribution {:}".format(mode))
    return r

def uGen(h):
    """
    Generates a u-value for the SPPT method from a given height. Used because random spread in atmosphere is less in the stratosphere.
    u is a factor that controls how much spread is at a given height

    Arguments:
    h: [float] given height to generate u

    Returns:
    u: [float] the generated u value
    """

    if h <= 300:
        u = 0
    elif 300 < h < 1300:
        u = np.sin((h - 300)*np.pi/2/1000)
    elif 15000 < h < 30000: #TEMPORARY VALUES (should decrease in stratosphere)
        u = abs(np.sin((h - 15000)*np.pi/15000 + np.pi/2))
    else:
        u = 1

    return u

def BMP(x):
    """ 
    BMP Perturbations

    Arguments:
    x: [ndarray] unperturbed variable

    Returns:
    X: [ndarray] perturbed variable
    """

    # Different random number per parameter

    r = rGen(mode='uniform')
    X = (1 + r)*x

    return X

def SPPT(x, h, r):
    """
    SPPT Perturbations

    Arguments:
    x: [ndarray] unperturbed variable
    h: [ndarray] height of the unperturbed variable
    r: [float] random perturbation number (from rGen())

    Returns:
    X: [ndarray] perturbed variable
    """

    # Same random numbers for all parameters

    X = np.zeros_like(x)

    for i in range(len(x)):
    
        u = uGen(h[i])
        X[i] = (1 + r*u)*x[i]

    return X

def temporal(x, x_l, x_u, r):
    """ Returns a variable x, randomly distrubuted between its temporal neighbours
        r = [-0.5, 0.5]

    Arguments:
    x: [ndarray] unperturbed variable
    x_l: [ndarray] unperturbed variable from the last time step
    x_u: [ndarray] unperturbed variable from the next times step
    r: [float] random perturbation number (from rGen())

    Returns:
    X: [ndarray] perturbed variable
    """
    if r <= 0:
        R = (r + 0.5)/0.5
        x_min = np.minimum(x_l, x)
        x_max = np.maximum(x_l, x)
        X = x_min + (x_max - x_min)*R
    else:
        R = r
        x_min = np.minimum(x, x_u)
        x_max = np.maximum(x, x_u)
        X = x_min + (x_max - x_min)*R
    return X

def spread(mean, spread, r):
    """ Returns a perturbed variable based off of a spread array

    Arguments:
    mean: [ndarray] original sounding data
    spread: [ndarray] ensemble spread data
    r: [float] random perturbation number (from rGen())

    Returns:
    X: [ndarray] perturbed variable
    """

    X = mean + r*mean*spread

    return X
def spread_r(mean, spread):
    """ Returns a perturbed variable based off of a spread array

    Arguments:
    mean: [ndarray] original sounding data
    spread: [ndarray] ensemble spread data

    Returns:
    X: [ndarray] perturbed variable
    """
    r = rGen(mode='spread')
    X = mean + r*mean*spread

    return X

def expandSounding(sounding):
    ''' Expands a sounding profile into usable variables
    '''
    if len(sounding) < 6:
      print("ERROR: sounding profile is not long enough. Check if atmosphere is enabled")

    #unpack sounding
    lats = sounding[0]
    lons = sounding[1]
    temps = sounding[2]
    mags = sounding[3]
    dirs = sounding[4]
    level = sounding[5]

    # Init the constants
    consts = Constants()

    #convert speed of sound to temp
    t = np.square(temps)*consts.M_0/consts.GAMMA/consts.R

    #convert to EDN
    dirs = np.radians(angle2NDE(np.degrees(dirs)))

    #convert mags and dirs to u and v
    u = mags*np.cos(dirs)
    v = mags*np.sin(dirs)

    return lats, lons, t, u, v, level

def readSpreadFile(setup, lat, lon, spread_file, line):
    """ reads a ensemble spread file from ECMWF
    """
    # Read the file
    dataset = Dataset(spread_file, "r+", format="NETCDF4")

    lat = roundToNearest(lat, 0.5)
    lon = roundToNearest(lon, 0.5)

    longitude = np.array(dataset.variables['longitude'])
    latitude = np.array(dataset.variables['latitude'])

    lon_index = int(np.where(longitude==lon)[0])
    lat_index = int(np.where(latitude==lat)[0])

    longitude = np.array(dataset.variables['longitude'][lon_index-10:lon_index+11])
    latitude = np.array(dataset.variables['latitude'][lat_index-10:lat_index+11])

    level = np.array(dataset.variables['level'])
    #pressure 1 - 1000 hPa , non-linear
    
    time = np.array(dataset.variables['time'])
    #not known

    start_time = int((setup.start_datetime.hour + np.round(setup.start_datetime.minute/60))%24)//3

    if line:
        # time, (number), level, lat, lon
        temperature = np.array(dataset.variables['t'][start_time, :, lat_index, lon_index])
        x_wind = np.array(dataset.variables['u'][start_time, :, lat_index, lon_index])
        y_wind = np.array(dataset.variables['v'][start_time, :, lat_index, lon_index])
    else:
        # time, (number), level, lat, lon
        temperature = np.array(dataset.variables['t'][start_time, :, lat_index-10:lat_index+11, lon_index-10:lon_index+11])
        x_wind = np.array(dataset.variables['u'][start_time, :, lat_index-10:lat_index+11, lon_index-10:lon_index+11])
        y_wind = np.array(dataset.variables['v'][start_time, :, lat_index-10:lat_index+11, lon_index-10:lon_index+11])
    
        # Repeat axis since spread has a resolution of 0.5 while the base data has a resolution of 0.25
        temperature = np.repeat(temperature, 2, axis=1)
        temperature = np.repeat(temperature, 2, axis=2)
        temperature = temperature[:, 0:-1, 0:-1]

        x_wind = np.repeat(x_wind, 2, axis=1)
        x_wind = np.repeat(x_wind, 2, axis=2)
        x_wind = x_wind[:, 0:-1, 0:-1]

        y_wind = np.repeat(y_wind, 2, axis=1)
        y_wind = np.repeat(y_wind, 2, axis=2)
        y_wind = y_wind[:, 0:-1, 0:-1]

    temperature = np.flip(temperature, axis=0)
    x_wind = np.flip(x_wind, axis=0)
    y_wind = np.flip(y_wind, axis=0)

    dataset.close()

    return temperature, x_wind, y_wind

def perturb(setup, sounding, method, sounding_u=[], sounding_l=[], spread_file='', lat=0, lon=0, line=False):
    ''' Takes in a sounding profile, and a method, and returns a perturbed sounding profile using that method

      inputs:
        sounding [ndarray] input sounding profile
        method [string] method to perturb with ('bmp' or 'sppt')

      outputs:
        sounding_p [ndarray] output perturbed sounding profile
    '''
    # Init the constants
    consts = Constants()
    
    lats, lons, t, u, v, level = expandSounding(sounding)

    if method == 'bmp':
        T = BMP(t)
        U = BMP(u)
        V = BMP(v)

    elif method == 'sppt':
        r = rGen(mode='normal')
        T = SPPT(t, level, r)
        U = SPPT(u, level, r)
        V = SPPT(v, level, r)

    elif method == 'temporal':
        r = rGen(mode='uniform')
        _, _, t_u, u_u, v_u, _ = expandSounding(sounding_u)
        _, _, t_l, u_l, v_l, _ = expandSounding(sounding_l)
        T = temporal(t, t_l, t_u, r)
        U = temporal(u, u_l, u_u, r)
        V = temporal(v, v_l, v_u, r)

    elif method == 'spread':
        r = rGen(mode='spread')
        if setup.debug:
            print("Perturbation Variable: {:}".format(r))
        spread_t, spread_u, spread_v = readSpreadFile(setup, lat, lon, spread_file, line)
        T = spread(t, spread_t, r)
        U = spread(u, spread_u, r)
        V = spread(v, spread_v, r)

    elif method == 'spread_r':
        spread_t, spread_u, spread_v = readSpreadFile(setup, lat, lon, spread_file, line)
        T = spread_r(t, spread_t)
        U = spread_r(u, spread_u)
        V = spread_r(v, spread_v)

    else:
        #print('WARNING: Unrecognized perturbation type, using SPPT')
        return None

    # convert temp to speed of sound
    TEMPS = (consts.GAMMA*consts.R/consts.M_0*T)**0.5

    # convert u and v to mags and dirs
    MAGS = np.sqrt(U**2 + V**2)
    DIRS = np.degrees(np.arctan2(V, U))

    #convert to NDE
    DIRS = angle2NDE(DIRS)

    #pack sounding
    sounding_p = [np.array(lats), np.array(lons), np.array(TEMPS), np.array(MAGS), np.array(DIRS), np.array(level)]
    try:
        np.save(os.path.join(setup.working_directory, setup.fireball_name, 'Perturbation_ID{:}'.format(r)), np.array(sounding_p))
    except:
        print('WARNING: Unable to save perturbed weather data!')
    return sounding_p

if __name__ == '__main__':
    
    mag_mean = np.array([  2.43213713 ,  3.63627388 ,  3.32134008 , 3.33470113 ,  3.79744476,
                   3.79345165,   3.34174992,   2.6870471 ,  2.28980092,   1.94628201,
                   1.43032892,   1.26067371,   2.95085383,  4.36697292,   6.19409193,
                   8.30043649,  10.24117654,  12.2603921 , 15.52678616,  19.59397738,
                  30.23930687,  33.28927378,  27.549643  , 23.02175367,  24.27312759,
                  28.2902178 ,  30.26223542,  26.91776968, 27.14700916,  34.7182345,
                  38.85260626,  31.90175241,  26.86746527, 13.5309995 ,  22.51310121,
                  35.37987001,  48.90065085,])


    mag_spread = np.array([ 0.48710405,  0.67396524,  0.75722182,  0.86092323,  1.0052456,   0.91466391,
                  0.92191542 , 1.06028298 , 1.1002522  , 0.97454604 , 0.8741689 ,  0.81977335,
                  0.61068341 , 0.61134561 , 0.70112861 , 0.83201352 , 1.00792929,  1.12437756,
                  1.32893088 , 1.46271835 , 1.05808837 , 0.76212984 , 0.8016715 ,  0.66790973,
                  0.68817111 , 0.70047501 , 0.73364396 , 0.72988745 , 0.70287717,  0.76249328,
                  1.15044044 , 1.50338121 , 1.06578162 , 1.28234778 , 0.97734093,  1.52042002,
                  0.94903383])


    mag_value = np.array([  2.15434646,   3.62424258,   8.9740961 ,  10.3734482 ,   9.29965552,
                 9.01347774,   9.89274613,  10.66685572,  11.63617126,  12.33539573,
                12.58136009,  11.7338793 ,  11.62526159,  13.5366879 ,  16.05015118,
                16.21558232,  15.59094118,  16.48556897,  20.95095013,  30.78167657,
                30.26922035,  26.40588906,  22.1552933 ,  18.75221037,  12.45458423,
                 7.82096244,   7.37549614,   9.54317454,   8.79280203,   5.4068898,
                 4.52304495,  23.21437718,  44.73954845,  63.28190157,  65.29010523,
                72.6312297 ,  60.9037868 ])

    h = np.array([48413.94 ,43738.55, 39850.56, 35586.89, 33763.05, 31330.96, 26635.56, 23900.02,
    20694.90, 18756.34, 16322.8, 15003.5, 13727.1, 12797.3, 11890.2, 11297.9, 10422.6 ,
    9255.70, 8380.36, 7214.09,6631.66, 5759.30 , 4896.02,4341.73, 3814.82 , 3089.25, 
    2654.69, 2261.80, 2081.09, 1750.63, 1459.91, 1328.70, 987.17, 798.72, 566.54, 334.24, 136.62])
    h = np.flip(h, axis=0)

    members =  np.array([[ 50.3755651 ,  33.04412292,  22.38593205,  14.34172003,  25.76589682,
   32.44457353,  39.17514784,  34.79677219,  27.62273496,  26.47415112,
   29.01602656,  27.7579047 ,  23.66288383,  23.15060518,  28.20283763,
   33.72358893,  30.12331954,  19.60190445,  16.25711325,  13.24941484,
   11.32027165,   9.367364  ,   7.02093133,   4.74308476,   3.04475079,
    1.58752045,   1.33443552,   2.00333505,   2.54830661,   3.50041708,
    4.36775209,   5.0751985 ,   5.18499799,   4.45519811,   4.03910337,
    3.63419532,   2.43930417],
 [ 48.05067532,  36.23836752,  24.11775878 , 14.13087569 , 26.99962126,
   32.88312835,  38.30650462,  35.2376642  , 27.50820468 , 27.40306704,
   30.54930621,  27.91458358,  24.71775279 , 22.85428347 , 26.7336385,
   33.16080322,  30.53993426,  21.3689116  , 16.73403585 , 12.4263575,
   10.25875648,   8.1558534 ,   5.83725371 ,  4.06426382 ,  2.68167936,
    0.91082047,   1.56108105,   2.15062182 ,  2.31637816 ,  2.38093897,
    2.97966223,   3.5322466 ,   3.38233262 ,  2.86717545 ,  3.19792147,
    3.98054793,   2.26104301],
 [ 48.3351335 ,  36.68189352,  22.81195858 , 13.00059081,  25.98425065,
   34.26408968,  37.52771335,  34.88599403 , 28.26034632,  26.46953305,
   30.56231222,  28.48696198,  23.784762   , 22.35158172,  28.62212895,
   33.27210034,  29.53134533,  16.91943237 , 13.96606004,  12.57964853,
   11.38529752,   9.15152124,   6.52969198 ,  4.49920514,   3.06267349,
    1.64229063,   1.3901147 ,   2.54168377 ,  2.87166764,   2.85364355,
    3.13456113,   3.73605123,   4.25377154 ,  3.82727974,   3.82183722,
    3.75841243,   2.39496701],
 [ 48.20907979,  34.79337193,  22.3116969  , 12.6503427 ,  28.11774897,
   31.81135229,  39.69747012,  34.96040368 , 26.94803714,  26.57286976,
   30.37396345,  29.09932962,  25.21109047 , 23.01612329,  27.27028579,
   33.33269802,  30.33046752,  18.94313186 , 14.35605132,  10.96408864,
    8.78483296,   6.89940083,   5.26493596 ,  4.09723825,   3.09987044,
    1.81916902,   1.2822816 ,   1.99020721 ,  2.6529959 ,   3.05542843,
    3.55219965,   4.02201918,   4.24357101 ,  3.43740041,   2.72410192,
    2.90156259,   1.89397374],
 [ 48.0858131 ,  34.84612181,  21.85036099,  13.26327988,  26.65752476,
   31.65184844,  39.72718635,  33.80653366,  26.83446736,  26.67943947,
   30.350472  ,  28.81778741,  23.64032428,  23.41292327,  27.99272513,
   33.64230493,  30.84039858,  20.69805313,  17.65406612,  13.50070042,
   10.39102298,   8.36883096,   6.13093599,   4.42350273,   3.35353922,
    1.58590003,   1.44519189,   1.77181185,   2.10967621,   2.48427426,
    3.33017079,   4.33123217,   5.04225703,   4.63580367,   4.48263164,
    4.04123113,   2.90032242],
 [ 48.39270912,  35.6538849 ,  23.23865938,  13.80327029,  26.88573596,
   31.198348  ,  38.52683977,  33.97206157,  26.87488763,  27.19945056,
   30.36372146,  27.86469821,  23.93897182,  23.13316693,  27.83061435,
   33.69262203,  30.48803712,  18.29219385,  15.31074965,  12.12472609,
   10.1295075 ,   8.11509752,   5.87320522,   3.94822326,   2.33955794,
    0.73086468,   2.35901047,   2.63153184,   3.23883011,   3.5092604,
    3.86174171,   3.99336338,   3.50807348,   2.74428715,   3.01773117,
    3.09316992,   2.53899811],
 [ 48.87993652,  34.765147  ,  23.45558201,  14.52328319,  27.61902453,
   29.83437161,  39.08124905,  35.23761247,  26.79486188,  26.90588725,
   30.79920619,  28.46127225,  24.51616133,  22.7301456 ,  27.38557677,
   32.80799039,  29.8676706 ,  20.06653376,  14.49132993,  11.43756286,
    9.84703277,   8.04574291,   5.9132213 ,   3.75892226,   2.3527966,
    0.44312587,   1.5916963 ,   1.77365341,   2.08590516,   2.50406052,
    3.21055132,   2.89125628,   2.53628014,   2.49592717,   2.92622214,
    3.49826733,   2.61464105],
 [ 49.65821561,  33.49693682,  21.39642327,  14.15795374,  28.10864127,
   30.48553112,  39.14040781,  34.72554918,  26.35494376,  27.75243422,
   30.08393469,  28.72240028,  25.07943815,  23.36123511,  27.16571003,
   34.24766406,  32.4655318 ,  21.1734381 ,  14.32770377,  10.83848335,
    9.33282834,   7.71372529,   6.09550782,   4.60670244,   3.15709253,
    2.03630623,   1.66855948,   1.4239507 ,   0.98447763,   1.57278936,
    3.01875117,   3.8276399 ,   3.57942209,   3.04160475,   2.22670924,
    3.53602657,   2.68928182],
 [ 49.73389235,  37.12174839,  21.72251853 , 13.32708827 , 25.8747371,
   32.25080414,  38.94498958,  34.68485281  ,26.98332144  ,26.96964818,
   30.82232032,  27.25660258,  23.80469383  ,23.08853324  ,27.47491424,
   32.63539051,  28.99480475,  19.27515756  ,15.78328603  ,12.37779163,
   10.56572086,   8.89150401,   7.07023256  , 5.16342055  , 3.52265195,
    1.86736764,   1.86252184,   2.89370332  , 3.39025375  , 3.48831309,
    3.58693201,   3.76646576,   3.75048566,   3.33453915,   3.486224,
    4.39360738,   2.59730872],
 [ 49.31353173,  37.16274269,  21.81390731,  12.57017916,  26.67686569,
   32.37542083,  38.40145708,  34.75834027,  27.19548242,  26.83803706,
   29.84355471,  28.69600103,  24.44487189,  23.20942695,  27.0197289,
   32.62818392,  29.3050529 ,  19.33634085,  16.28565018,  12.98394192,
   10.42209098,   8.34993039,   6.24839145,   4.36403331,   2.72752118,
    0.86885893,   1.87050204,   2.22479481,   2.57088265,   3.15255483,
    3.8616398 ,   3.811982  ,   3.48426201,   3.40261711,   3.69020487,
    4.04642728,   2.50462704]])

    for i in range(1000):
        B = BMP(mag_value)
        if i == 1:
            plt.plot(B, h, c='g', alpha=0.2, linewidth=3, label='BMP Scheme')
        else:
            plt.plot(B, h, c='g', alpha=0.2, linewidth=3)

    for i in range(1000):
        S = SPPT(mag_value, h)
        if i == 1:
            plt.plot(S, h, c='r', alpha=0.2, linewidth=3, label='SPPT Scheme')
        else:
            plt.plot(S, h, c='r', alpha=0.2, linewidth=3)

    for i in range(1000):
        Sp = spread(mag_mean, mag_spread)
        if i == 1:
            plt.plot(Sp, h, c='c', alpha=0.2, linewidth=3, label='Spread of Mean')
        else:
            plt.plot(Sp, h, c='c', alpha=0.2, linewidth=3)

    for i in range(1000):
        Sp2 = spread(mag_value, mag_spread)
        if i == 1:
            plt.plot(Sp2, h, c='m', alpha=0.2, linewidth=3, label='Spread of Original Values')
        else:
            plt.plot(Sp2, h, c='m', alpha=0.2, linewidth=3)

    for i in range(10):
        if i == 1:
            plt.plot(np.flip(members[i], axis=0), h, c='orange', alpha=0.2, linewidth=3, label='Ensemble Members')
        else:
            plt.plot(np.flip(members[i], axis=0), h, c='orange', alpha=0.2, linewidth=3)  

    plt.plot(mag_value, h, c='k', linewidth=3, label='Original Values')
    plt.legend(loc='lower right')

    plt.show()
    # X = BMP(x)
    # print(X)