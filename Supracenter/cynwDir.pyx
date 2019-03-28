"""Finds component of wind normal to direction of travel"""

import cython
cimport cython
import numpy as np
cimport numpy as np

FLOAT_TYPE = np.float64
ctypedef np.float64_t FLOAT_TYPE_t

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef np.ndarray[FLOAT_TYPE_t, ndim=1] nwDir(np.ndarray[FLOAT_TYPE_t, ndim=1] w, np.ndarray[FLOAT_TYPE_t, ndim=1] wd, np.ndarray[FLOAT_TYPE_t, ndim=1] phi):
    # Author: Wayne Edwards

    """ Calculates the component of the wind vector in the direction of phi for each layer in the atmosphere

    Arguments:
        w: [ndarray] Wind Magnitude for each level in the atmosphere
        wd: [ndarray] Wind direction (angle). Direction (from the north) the wind is blowing to
        phi: [ndarray] Azimuthal angles (NDE)

    Returns:
        x: [ndarray] Component of the wind in the direction of phi
    """

    cdef:
        int m = len(w)
        int n = len(phi)

        # Calculate unit normals in direction phi
        np.ndarray[FLOAT_TYPE_t, ndim=1] nx = np.cos(phi) #changed from cos
        np.ndarray[FLOAT_TYPE_t, ndim=1] ny = np.sin(phi) #changed from sin
        np.ndarray[FLOAT_TYPE_t, ndim=1] wx = w*np.sin(wd) #changed back to positive, cyscan equations are defined on the direction wind is blowing from
        np.ndarray[FLOAT_TYPE_t, ndim=1] wy = w*np.cos(wd) 

        np.ndarray[FLOAT_TYPE_t, ndim=2] nx_t = np.tile(nx, (m, 1))
        np.ndarray[FLOAT_TYPE_t, ndim=2] ny_t = np.tile(ny, (m, 1))

        # Calculate wind components x-y
        np.ndarray[FLOAT_TYPE_t, ndim=2] wx_t = np.tile(wx, (n, 1)).T   
        np.ndarray[FLOAT_TYPE_t, ndim=2] wy_t = np.tile(wy, (n, 1)).T

        # Find component of wind velocity in the direction of phi for each level z
        np.ndarray[FLOAT_TYPE_t, ndim=2] x = wx_t*nx_t + wy_t*ny_t
    
    return x