"""Splices atmospheric profile to only include data between their heights"""

import cython
cimport cython
import numpy as np
cimport numpy as np

FLOAT_TYPE = np.float64
ctypedef np.float64_t FLOAT_TYPE_t

from supra.Supracenter.bisearch import bisearch

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef zInteg(long d, long s, region):
    # Author: Wayne Edwards

    """ Returns a subset of region, interval, that the sound waves propegate through, from S to d. The region
        will be extended to s and d if it is not large enough. The layers are kept he same, no interpolation is
        done

    Arguments:
        d: [int] height of detector
        s: [int] height of source
        region: [ndarray] full atmospheric data as given by the user

    Returns:
        interval: [ndarray] atmospheric data from between s and d
    """
    
    # flips heights if there is an error
    if d > s:
        d, s = s, d

    # initialize variables
    cdef:
        int n = 0
        int m = 0
    
    # Dimensions of region
    n, m = region.shape

    # Find position of S & D in the integration region
    cdef:
        int i = bisearch(region[:, 0], d)
        int j = bisearch(region[:, 0], s) - 1

    # interval is subset of region between D and S
    interval = region[i:j + 2, 0:m]

    # Dimensions of interval
    n, m = interval.shape

    # Adjust the interval to properly include the source & destination
    # Change the upper and lower bounds to match S and D
    # Make detector lowest height
    # Detector is higher than lowest height
    if (d >= region[i, 0]):
        interval[0, 0] = d

    # Detector is lowest height
    else: 
        # Copy the lowest row if detector is lower than range
        interval[1:n, 0:m] = interval[0:n - 1, 0:m]
        interval[0, 0] = d
        n += 1;

    # Source is the heighest height
    if (s >= region[j, 0]):

        # Add source location at end of profile
        interval = np.vstack((interval, np.array([s, *interval[-1, 1:]])))

        # If the top region is already where the source is, then get rid of the duplicate
        if (s == region[j, 0]) and (j + 1 != len(region)):
            interval=np.delete(interval, -1, axis=0)

    return interval


def zInterp(d, s, region, div=100):
    f = 0
    max_indx = 0
    min_indx = 0
    r = d - s
    dr = r/(div)
    column_list = []
    
    for i in range(div):
        h = d - i*dr
        row_list = []
        
        for j in range(len(region[:, 0])):
            if h - region[j, 0] < 0:
                max_indx = j
                min_indx = j - 1
                f = (h - region[j - 1, 0])/(region[j, 0] - region[j - 1, 0])
                break

        for e in range(len(region[0])):
            if e > 0:
                row_list.append(f*(region[max_indx, e] - region[min_indx, e]) + region[min_indx, e])
            else:
                row_list.append(h)
        
        column_list.append(row_list)

    region = np.concatenate((region, np.array(column_list)), axis=0)
    region = np.unique(region, axis=0)
    region = region[region[:,0].argsort()]
    return region

if __name__ == "__main__":

    # Test case

    # Source
    s_h = 2001

    # Detector
    d_h = 1.2

    # Sample weather profile
    
                       #[Height (m), Temp (K), wSpd (m/s), wDir (deg fN)]
    zProfile =np.array([[     0    , 293.2   , 2.058     , 0.000        ],
                        [   500    , 296.7   , 3.892     , 0.463        ],
                        [  1000    , 290.3   , 4.630     , 0.783        ],
                        [  1500    , 289.3   , 4.760     , 1.245        ],
                        [  2000    , 289.4   , 9.577     , 1.543        ]])

    # Sample function call
    region = zInteg(d_h, s_h, zProfile)

    print(region)
