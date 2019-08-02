
# Based off of Wayne Edwards' ray-tracing algorithm (2003)               
# Finds path between two locations with an atmospheric profile in between #
###########################################################################

import cython
cimport cython
import warnings

import time

import numpy as np
cimport numpy as np
from libc.math cimport sqrt, M_PI, M_PI_2, atan2
 
# Define cython types for numpy arrays
FLOAT_TYPE = np.float64
ctypedef np.float64_t FLOAT_TYPE_t

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def mag(vect):
    return np.sqrt(np.dot(vect, vect))

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def norm(vect):
    return vect/mag(vect)

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
def speedscan(x, *args):
  # unpack all values

    
    S, D, z_profile = args

    phi = x[0]
    theta = x[1]

    propagating = True

    current_layer = 0
    max_layer = len(z_profile)
    
    p_0 = S
    a_theta = np.sin(theta)
    b_theta = np.cos(theta)
    
    z = z_profile[:, 0]
    c = z_profile[:, 1]

    while propagating:

        if phi > M_PI_2:
            direction = 1
        elif phi < M_PI_2:
            direction = -1

        a_phi = np.sin(phi)
        b_phi = np.cos(phi)

        dz = z[current_layer + direction] - z[current_layer]

        n = np.array([a_theta*a_phi, b_theta*a_phi, b_phi])

        u = z_profile[current_layer, 2]
        v = z_profile[current_layer, 3] 

        V = c[current_layer]*n + np.array([u, v, 0])

        p_0 = dz*V/V[2] + p_0

        c_0 = mag(V)

        c_s = z_profile[current_layer + direction, 1]
        u = z_profile[current_layer + direction, 2]
        v = z_profile[current_layer + direction, 3]

        c_1 = mag(c_s*n + np.array([u, v, 0]))

        if (c_1/c_0*np.sin(phi) > 1): # wave reflects
            error = mag(D - p_0)

            return error
        else: # wave refracts
            phi_1 = M_PI - np.arcsin(c_1/c_0*np.sin(phi))
            current_layer += direction


        error = mag(D - p_0)
        
        if current_layer + 1 == max_layer or p_0[2] < 0:
            propagating = False

    return error
