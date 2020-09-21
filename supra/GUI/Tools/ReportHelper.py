import numpy as np
from supra.Utils.AngleConv import chauvenet
def latitudify(lat):
    
    if lat > 0:
        lat = '{:}°N'.format(lat)
    else:
        lat = '{:}°S'.format(np.abs(lat))

    return lat

def longitudify(lon):

    if lon > 0:
        lon = '{:}°E'.format(lon)
    else:
        lon = '{:}°W'.format(np.abs(lon))

    return lon

def byteify(b):

    count = 0

    while b > 1024:
        b /= 1024
        count += 1

    # Tb is optimistic
    si_fix = ['b', 'kb', 'Mb', 'Gb', 'Tb']

    return '{:.3f} {:}'.format(b, si_fix[count])

def obtainPerts(data, frag):
    data_new = []

    for i in range(len(data[frag][1])):
        data_new.append(data[frag][1][i][0])
    data, remove = chauvenet(data_new)

    return data, remove