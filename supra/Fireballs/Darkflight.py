
import subprocess
import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import simplekml

from wmpl.Utils.PlotMap import GroundMap

try:
    # Python 2  
    import ConfigParser as configparser

except:
    # Python 3
    import configparser

def halfMags(low, high, div):
    """ Returns a list of magnitude values and half-magnitude values between low and high

    Example: low = 0.049 and high = 100
    will return: [0.049, 0.050, 0.1, 0.5, 1, 5, 10, 50, 100] 

    Arguments:
        low: [float] lowest value to be in list
        high: [float] highest value to be in list
        div: [int] Spacing between intermediate magnitudes is log10(div)
                    div = 10 - every magnitude
                    div = 5  - every magnitude and half magnitude

    Returns:
        [list] list of values between low and high (inclusive) with half-magnitudes between
    """

    # Limits of the mass
    h = np.log10(high)
    l = np.log10(low)

    # Intermediate masses
    
    # every magnitude
    a = np.round((np.arange(l, h)))

    # every partial magnitude above each magnitude
    b = np.round((np.arange(l, h))) + np.log10(div)

    # every partial magnitude below each magnitude
    c = np.round((np.arange(l, h))) - (1 - np.log10(div))

    # combine all magnitudes and remove duplicate
    s = np.sort(np.concatenate((np.array([l, h]), a, b, c)))

    temp = s[s <= h]
    mass_list = temp[temp >= l]

    return np.sort(list(set(10**mass_list)))


def massString(mass):
    """ HELPER FUNCTION. Turns a float mass into a string with the correct units after
    (between mg and kg)

    Arguments:
        mass: [float] the mass to be converted

    Returns:
        mass_str: [string] the mass as a string with the proper units
    """     

    # make mass into str
    if mass < 1:
        mass *= 1000
        if mass < 1:
            mass *= 1000
            mass_str = str(eval("%.0e" % (mass))) + ' mg'
        else:
            mass_str = str(eval("%.0e" % (mass))) + ' g'
    else:
        mass_str = str(eval("%.0e" % (mass))) + ' kg'

    return mass_str

def runDarkflight(params, mass, uncertainty=False):
    """ Wrapper function for Darkflight
        Predicts the travel path of the meteor in darkflight given the initial conditions. Uses Monte Carlo to 
        find uncertainty in locations. Exports a .txt file to params.output which is read later in the program 
        containing the solution to the darkflight.

    Arguments:
        params: [Object] contains all .ini variables, including all inputs into the Darkflight code
        mass: [float] current mass to be calculated (this function runs multiple times for various masses)
        uncertainty: [boolean] True - Calculate the uncertainty mode
                               False - Calculate the darkflight path
    """

    # Name file by mass
    mass_str = massString(mass)

    # Convert all variables to strings
    lat = str(params.x)
    lon = str(params.y)
    elev = str(params.z)
    v = str(params.v)
    az = str(params.az)
    ze = str(params.ze)
    mass = str(mass)
    tan = str(params.tan)
    sig = str(params.sig)
    den = str(params.den)
    end = str(params.end)
    shape = str(params.shape)
    dlat = str(params.dlat)
    dlon = str(params.dlon)
    dh = str(params.dh)
    dv = str(params.dv)
    daz = str(params.daz)
    dze = str(params.dze)
    max_iter = str(params.max_iter)
    drg = str(params.drag)

    # Directory of Darkflight
    dir1 = "/home/dvida/source/darkflight"
    dir2 = "C:\\Users\\lmcfd\\Desktop\\"
    dir3 = "/srv/meteor/fireballs/darkflight"
    if os.path.isdir(dir1):
        darkflight_dir = dir1
    elif os.path.isdir(dir2):
        darkflight_dir = dir2
    elif os.path.isdir(dir3):
        darkflight_dir = dir3
    else:
        print("The directory with the darkflight code cannot be found!")
        print("Tried these directories:")
        print(dir1)
        print(dir2)
        sys.exit()

    # Run Darkflight through terminal calls

    print()
    print("Running the darkflight binary for the {:.3f} kg mass bin...".format(float(mass)))

    # darkflight.c help menu
    if params.help:
        p = subprocess.Popen(['./darkflight', '--help'], cwd=darkflight_dir)
        exit()

    # Monte-Carlo simulations
    if uncertainty:
        # Fragmentation simulation
        if params.frag:
            p = subprocess.Popen(['./darkflight', '--src', lat + ',' + lon + ',' + elev, '--vel', v, '--az', az, '--zn', ze, '--mas', mass, \
                       '--tan', tan, '--sig', sig, \
                       '--den', den, '--end', end, '--shp', shape, \
                       '--dsp', dlat + ',' + dlon + ',' + dh + ',' + dv, \
                       '--dra', daz + ',' + dze, \
                       '--itr', max_iter, '--drg', drg, '--frag',\
                        '--atm', params.atm,'--out',  params.output + '_uncertainty' + '_' + mass_str], cwd=darkflight_dir)
        
        # no fragmentation
        else:
            p = subprocess.Popen(['./darkflight', '--src', lat + ',' + lon + ',' + elev, '--vel', v, '--az', az, '--zn', ze, '--mas', mass, \
                       '--tan', tan, '--sig', sig, \
                       '--den', den, '--end', end, '--shp', shape, \
                       '--dsp', dlat + ',' + dlon + ',' + dh + ',' + dv, \
                       '--dra', daz + ',' + dze, \
                       '--itr', max_iter, '--drg', drg,\
                        '--atm', params.atm,'--out',  params.output + '_uncertainty' + '_' + mass_str], cwd=darkflight_dir)

    # Best trajectory points
    else:


        p = subprocess.Popen(['./darkflight', '--src', lat + ',' + lon + ',' + elev, '--vel', v, '--az', az, '--zn', ze, '--mas', mass, \
                       #'--tan', tan, '--sig', sig, \
                       '--den', den, '--end', end, '--shp', shape, \
                       #'--dsp', dlat + ',' + dlon + ',' + dh + ',' + dv, \
                       #'--dra', daz + ',' + dze, \
                       '--itr', max_iter, '--drg', drg, \
                        '--atm', params.atm,'--out',  params.output + '_' + mass_str], cwd=darkflight_dir)

    # Run terminal call
    p.communicate()

def readDarkflight(output, mass, header=30, partial=False):
    """ Reads output from a darkflight output.txt file, and adds a mass "stamp" to each line of data, so that
        each point can be associated with a mass

    Arguments:
        output: [string] location of output file
        mass: [float] mass stamp to be added to the point
        header: [int] number of header lines in the file
        partial: [Boolean] whether to take all data (False), or every 20th line (True).
    """

    with open(output) as f:

        for i in range(header):
            next(f)

        data = np.array([0]*10)
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '')

            # Split the line by the delimiter
            line = line.split()

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = float(entry.strip())
            line.append(mass)

            # Add the contents of the line to the data list
            data = np.vstack((data, np.array(line)))

        # First row was all zeroes
        data = np.delete(data, 0, 0)

        if partial:
            # Take every 20th line
            return data[1::20, :]

        else:

            # Take every line
            return data


def plotData(data, del_data, params, mass_list):
    """ Function to plot the darkflight data, and create the .kml file

    Arguments:
        data: [ndarray] array containing parsed data from darkflight output
        del_data: [ndarray] array containing parsed data from darkflight uncertainty output
        params: [Object] object containing all .ini configuration variables
        mass_list: [list] list of masses used in simulation

    """

    # Initialize lists
    lat = []
    lon = []
    alt = []
    vel = []
    l   = []
    l_x = []
    t   = []
    az  = []
    zen = []
    mas = []
    lat_long = []
    lon_long = []
    del_lat  = []
    del_lon  = []
    del_alt  = []
    del_t    = []

    # Each index for each mass
    for i in range(len(mass_list)):

        # Skip NaN values
        if np.isnan(data[i][:, 0].all()) or np.isnan(data[i][:, 1].all()):
            continue

        lat.append(data[i][:, 0])
        lon.append(data[i][:, 1])
        alt.append(data[i][:, 2])
        vel.append(data[i][:, 3])
        l.append(data[i][:, 4])
        l_x.append(data[i][:, 5])
        t.append(data[i][:, 6])
        az.append(data[i][:, 7])
        zen.append(data[i][:, 8])
        mas.append(data[i][:, 9])

        if params.error:
            del_lat.append(del_data[i][:, 0])
            del_lon.append(del_data[i][:, 1])
            del_alt.append(del_data[i][:, 2])
            del_t.append(del_data[i][:, 6])


    for i in range(len(mass_list)):
    
        lat_long.extend(lat[i])
        lon_long.extend(lon[i])


    lat_long = np.array(lat_long)
    lon_long = np.array(lon_long)

    # Set up axes
    # 3D Plot
    fig = plt.figure(figsize=plt.figaspect(0.5))
    fig.set_size_inches(20.9, 11.7)

    ax1 = fig.add_subplot(1, 1, 1, projection='3d')

    for i in range(len(mass_list)):
        ax1.plot(lat[i], lon[i], alt[i])
        ax1.text(lat[i][-1], lon[i][-1], alt[i][-1], '{:10.1E} kg'.format(mas[i][-1]), size=6)
        if params.error:
            ax1.scatter(del_lat[i], del_lon[i], del_alt[i], depthshade=False)

    ax1.set_xlabel('Latitude (N)')
    ax1.set_ylabel('Longitude (E)')
    ax1.set_zlabel('Elevation (km)')
    ax1.set_title('Darkflight Path')

    plt.savefig(params.output + '_darkflight_path.png', dpi=300)

    plt.show()


    #.kml file
    mass = [0]*len(mass_list)
    kml=simplekml.Kml()
    # traj = kml.newfolder(name='Fall Path')
    # uncs = kml.newfolder(name='Monte Carlo')
    # poin = kml.newfolder(name='Best Points')
    pnt = [None]*len(mass_list)

    if params.error:
        # Trajectory points
        for j in range(len(mass_list)):
            mass[j] = kml.newfolder(name=massString(mass_list[j]))
            for i in range(len(del_lat[j])):
            # pnt[j] = traj.newpoint(coords=[(lon[j][i], lat[j][i], alt[j][i]*1000)], altitudemode='clampToGround')

            # if j%10 == 0:
            #     pnt[j].style.iconstyle.color = simplekml.Color.blue
            
            # elif j%10 == 1:
            #     pnt[j].style.iconstyle.color = simplekml.Color.orange

            # elif j%10 == 2:
            #     pnt[j].style.iconstyle.color = simplekml.Color.green

            # elif j%10 == 3:
            #     pnt[j].style.iconstyle.color = simplekml.Color.red
    
            # elif j%10 == 4:
            #     pnt[j].style.iconstyle.color = simplekml.Color.purple

            # elif j%10 == 5:
            #     pnt[j].style.iconstyle.color = simplekml.Color.brown

            # elif j%10 == 6:
            #     pnt[j].style.iconstyle.color = simplekml.Color.pink

            # elif j%10 == 7:
            #     pnt[j].style.iconstyle.color = simplekml.Color.grey

            # elif j%10 == 8:
            #     pnt[j].style.iconstyle.color = simplekml.Color.yellow

            # elif j%10 == 9:
            #     pnt[j].style.iconstyle.color = simplekml.Color.black
            
            # pnt[j].style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'


                # Monte Carlo points
                pnt[j] = mass[j].newpoint(coords=[(del_lon[j][i], del_lat[j][i], del_alt[j][i]*1000)], altitudemode='clampToGround')

                if j%10 == 0:
                    pnt[j].style.iconstyle.color = simplekml.Color.blue
                
                elif j%10 == 1:
                    pnt[j].style.iconstyle.color = simplekml.Color.orange

                elif j%10 == 2:
                    pnt[j].style.iconstyle.color = simplekml.Color.green

                elif j%10 == 3:
                    pnt[j].style.iconstyle.color = simplekml.Color.red
        
                elif j%10 == 4:
                    pnt[j].style.iconstyle.color = simplekml.Color.purple

                elif j%10 == 5:
                    pnt[j].style.iconstyle.color = simplekml.Color.brown

                elif j%10 == 6:
                    pnt[j].style.iconstyle.color = simplekml.Color.pink

                elif j%10 == 7:
                    pnt[j].style.iconstyle.color = simplekml.Color.grey

                elif j%10 == 8:
                    pnt[j].style.iconstyle.color = simplekml.Color.yellow

                elif j%10 == 9:
                    pnt[j].style.iconstyle.color = simplekml.Color.black

                pnt[j].style.iconstyle.icon.href = 'http://maps.google.com/mapfiles/kml/shapes/placemark_circle.png'

    for j in range(len(mass_list)):
        mass[j] = kml.newfolder(name=massString(mass_list[j]))
        
        best_t = t[j][-1]
        
        if params.error:
            min_t = np.nanmin(del_t)
            max_t = np.nanmin(del_t)
            up_t = abs(max_t - best_t)
            down_t = abs(min_t - best_t)
            dt = max(up_t, down_t) 

        #Best ground points
        if params.enable_time and params.error:
            pnt[j] = mass[j].newpoint(coords=[(lon[j][-1], lat[j][-1], alt[j][-1]*1000)], altitudemode='clampToGround',\
                name='{:}, {:5.2f} +/- {:3.2f} s'.format(massString(mas[j][-1]), best_t, dt))
        elif params.enable_time:
            pnt[j] = mass[j].newpoint(coords=[(lon[j][-1], lat[j][-1], alt[j][-1]*1000)], altitudemode='clampToGround',\
                name='{:}, {:5.2f} s'.format(massString(mas[j][-1]), best_t))
        else:
            pnt[j] = mass[j].newpoint(coords=[(lon[j][-1], lat[j][-1], alt[j][-1]*1000)], altitudemode='clampToGround',\
                name='{:}'.format(massString(mas[j][-1])))


    
    kml.save(params.output + '_darkflight_map.kml')



    # 2D map plot
    print(np.radians(lat_long), np.radians(lon_long))
    m = GroundMap(np.radians(lat_long), np.radians(lon_long), border_size=5, color_scheme='light')

    for i in range(len(mass_list)):
        m.plot(np.radians(lat[i]), np.radians(lon[i]))
        x, y = m.m(lon[i][-1], lat[i][-1])
        plt.text(x, y, '{:10.1E} kg'.format(mas[i][-1]), size=6)

    #m.scatter(np.radians(lat[-1]), np.radians(lon[-1]), marker='*', c='r')

    plt.savefig(params.output + '_darkflight_map.png', dpi=300)

    plt.show()


def readInfile(infile):
    """ Helper function to parse input .ini file, and return it in an object params
    """

    class Config:

        def __init__(self):
            pass

    params = Config()

    config = configparser.ConfigParser()

    if not os.path.exists(infile):
        print('The input file: ', infile, 'does not exist!')
        sys.exit()

    config.read(infile)

    try:
        params.output =  config.get('General', 'output')
        params.atm    =  config.get('General', 'atm')
        params.error  =  config.get('General', 'error')
    except:
        print('INI ERROR: [General] expected variables: output, atm, error')
        sys.exit()

    try:
        params.help = config.get('General', 'help')
        params.help = (params.help.lower() == 'true')
    except:
        params.help = False

    try: 
        params.enable_time = config.get('Tweaks', 'enable_time')
        params.enable_time = (params.enable_time.lower() == 'true')
    except:
        params.enable_time = False

    try:
        params.mags = float(config.get('Tweaks', 'mags'))
    except:
        params.mags = 5

    try:
        params.x =       float(config.get('Variables', 'lat'))
        params.y =       float(config.get('Variables', 'lon'))
        params.z =       float(config.get('Variables', 'z'))
        params.v =       float(config.get('Variables', 'v'))
        params.az =      float(config.get('Variables', 'az'))
        params.ze =      float(config.get('Variables', 'ze'))

    except:
        print('INI ERROR: [Variables] expected variables: x, y, z, v, az, ze - must be strictly floats!')
        sys.exit()

    try:
        params.mass_min =    float(config.get('Variables', 'mass_min'))
    except:
        # 1 gram
        params.mass_min = 0.001
        print('Using default minimum mass: {:.3f} kg'.format(params.mass_min))

    try:
        params.mass_max =    float(config.get('Variables', 'mass_max'))
    except:
        params.mass_max = 25
        print('Using default mass: {:.3f} kg'.format(params.mass_max))
    
    try:
        params.tan =     float(config.get('Variables', 'tan'))
    except:
        params.tan  = 0.0
        print('Using default meteoroid velocity perpendicular to radiant: {:.1f} km/s'.format(params.tan))

    try:
        params.sig =     float(config.get('Variables', 'sig'))
    except:
        params.sig  = 0.0
        print('Using default perpendicular velocity gaussian spread: {:.1f} km/s'.format(params.sig))

    try:
        params.frag =    config.get('Variables', 'frag')
        params.frag =    (params.frag.lower() == 'true')
    except:
        print('Using default: no fragmentation solution')    
        params.frag = False

    try:
        params.den =     float(config.get('Variables', 'den'))
    except:
        params.den = 3700
        print('Using default density: {:.1f} kg/m^3 (Stoney Chrondrite)'.format(params.den))

    try:
        params.end =     float(config.get('Variables', 'end'))
    except:
        params.end = 0
        print('Using default stopping elevation: {:.1f} km'.format(params.end))

    try:
        params.shape =   config.get('Variables', 'shape')
        if params.shape.lower() == 'sphere':
            params.shape = 0
            print('SHAPE: SPHERE')
        elif params.shape.lower() == 'brick':
            params.shape = 1
            print('SHAPE: BRICK')
        elif params.shape.lower() == 'cone':
            params.shape = 2
            print('SHAPE: CONE')
        else:
            print('Unknown shape, use Sphere, Brick or Cone. Using default shape: Brick')
            params.shape = 1
    except:
        print('Using default shape: Brick')
        params.shape = 1

    try:
        params.dlat = float(config.get('Variables', 'dlat'))
    except:
        params.dlat = 0.001
        print('Using default latitude uncertainty: {:.3f} deg'.format(params.dlat))

    try:
        params.dlon = float(config.get('Variables', 'dlon'))
    except:
        params.dlon = 0.001
        print('Using default longitude uncertainty: {:.3f} deg'.format(params.dlon))

    try:
        params.dh = float(config.get('Variables', 'dh'))
    except:
        params.dh = 1
        print('Using default height uncertainty: {:.3f} km'.format(params.dh))

    try:
        params.dv = float(config.get('Variables', 'dv'))
    except:
        params.dv = 1
        print('Using default velocity uncertainty: {:.3f} km/s'.format(params.dv))

    try:
        params.daz = float(config.get('Variables', 'daz'))
    except:
        params.daz = 0.1
        print('Using default azimuth uncertainty: {:.3f} deg'.format(params.daz))

    try:
        params.dze = float(config.get('Variables', 'dze'))
    except:
        params.dze = 0.1
        print('Using default zenith uncertainty: {:.3f} deg'.format(params.daz))

    try:
        params.max_iter = float(config.get('Variables', 'max_iter'))
    except:
        params.max_iter = 200

    try:
        params.drag = int(config.get('Variables', 'drag'))
    except:
        params.drag = 0


    return params

if __name__ == "__main__":

    ### COMMAND LINE ARGUMENTS

    # Init the command line arguments parser
    arg_parser = argparse.ArgumentParser(description="""
                ~~Darkflight~~ 
    Estimate darkflight path of meteor given
    initial conditions.

    """,
        formatter_class=argparse.RawTextHelpFormatter)

    arg_parser.add_argument('input_file', type=str, help='Path to Supracenter input file.')

    # Parse the command line arguments
    cml_args = arg_parser.parse_args()

    params = readInfile(cml_args.input_file)

    # What masses are chosen
    mass_list = halfMags(params.mass_min, params.mass_max, params.mags)
    #mass_list = [1, 0.1, 0.05, 0.02, 0.01, 0.001]

    print("Masses used:", mass_list)

    params.error = (params.error.lower() == 'true')

    if not os.path.exists(params.output):
        os.makedirs(params.output)

    print('start')
    data = []
    del_data = []
    for i in range(len(mass_list)):

        mass_str = massString(mass_list[i])

        runDarkflight(params, mass_list[i], uncertainty=False)
        data.append(readDarkflight(params.output + '_' + mass_str, mass_list[i], header=34))

        if params.error:
            runDarkflight(params, mass_list[i], uncertainty=True)
            del_data.append(readDarkflight(params.output + '_uncertainty' + '_' + mass_str, mass_list[i], header=34, partial=False))
    print("end")

    plotData(data, del_data, params, mass_list)