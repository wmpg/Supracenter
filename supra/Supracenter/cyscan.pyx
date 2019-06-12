#    ____                                __         
#   / __/_ _____  _______ ________ ___  / /____ ____
#  _\ \/ // / _ \/ __/ _ `/ __/ -_) _ \/ __/ -_) __/
# /___/\_,_/ .__/_/  \_,_/\__/\__/_//_/\__/\__/_/   
#         /_/   
#
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

from supra.Supracenter.cynwDir import nwDir
 
# Define cython types for numpy arrays
FLOAT_TYPE = np.float64
ctypedef np.float64_t FLOAT_TYPE_t

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef np.ndarray[FLOAT_TYPE_t, ndim=1] cyscan(np.ndarray[FLOAT_TYPE_t, ndim=1] supra_pos, np.ndarray[FLOAT_TYPE_t, ndim=1] detec_pos, 
    np.ndarray[FLOAT_TYPE_t, ndim=2] z_profile, wind=True, int n_theta=45, int n_phi=90, float precision=1e-5, float tol=1000):
    # switched positions (Jun 2019)

    # This function should be called for every station
    # Original Author: Wayne Edwards

    """ Finds an optimal ray from the source of the Supracenter to the detector, by making a guess, 
        and checking for the angle of minimum error. This is repeated with better guesses by taking angles around the best angle of the previous step.

    Arguments:  
        supra_pos: [array] 3-D local coordinates of the source of the sound. Given as: [x, y, Elevation]
        detec_pos: [array] 3-D local coordinates of the detector. Given as: [x, y, Elevation]
        zProfile: [array] The given atmospheric profile between the source and the detector. Given as: [Height, Temperature, Wind Speed, Wind Direction] for different values of height
        wind: [int] 0 - disable winds, 1 - enable winds. Temperature is still used.
        n_theta, n_phi: [int] angle grid spacing of the takeoff, azimuth angles
        tol: [float] Tolerance on the error of the takeoff angle
        precision: [double] minimum resolution of angles

    Returns:    
        t_arrival: [float] Direct arrival time between source and detector
        Azimuth: [float] The initial azimuthal angle for the source to the detector in degrees
        Takeoff: [float] The initial takeoff angle from the source to the detector in degrees

    See diagram on pg 34 of SUPRACENTER for more information
    """

    # Azimuths and Wind directions are measured as angles from north, and increasing clockwise to the East
    
    
    # Switch to turn off winds
    if not wind:
        z_profile[:, 2] = 0

    found = False

    ### Initialize variables ###
    cdef:
        # Initial grid spacing
        # Azimuth angle (0 - 360 degrees from North due East)
        float dtheta = M_PI/n_theta
        
        # Takeoff angle (90 - 180 degrees from vertical)
        float dphi = M_PI_2/n_phi   

        # Horizontal distance between source and a station.
        float dx = detec_pos[0] - supra_pos[0]
        float dy = detec_pos[1] - supra_pos[1]

        # azth - initial guess for azimuth
        float azth = atan2(dy, dx)

        # The number of layers in the integration region
        int n_layers = len(z_profile)

        # Slowness, as defined in SUPRACENTER on pg 35, s = 1/c
        np.ndarray[FLOAT_TYPE_t, ndim=1] s = 1.0/z_profile[0:n_layers, 1]

        # Elevation for that layer
        np.ndarray[FLOAT_TYPE_t, ndim=1] z  = z_profile[0:n_layers, 0]

        # Set up grid of angles
        np.ndarray[FLOAT_TYPE_t, ndim=1] phi = np.linspace(azth-M_PI_2, azth+M_PI_2, n_phi)
        np.ndarray[FLOAT_TYPE_t, ndim=2] Phi = np.tile(phi, (n_theta, 1))
    
        # Component of wind vector in the direction of phi and phi + pi/2 respectively
        np.ndarray[FLOAT_TYPE_t, ndim=2] u = nwDir(z_profile[:, 2], z_profile[:, 3], phi)
        np.ndarray[FLOAT_TYPE_t, ndim=2] v = nwDir(z_profile[:, 2], z_profile[:, 3], phi+M_PI_2)

        # Construct ray parameter net
        # Theta - the take-off angle, of the ray from the source to the station
        np.ndarray[FLOAT_TYPE_t, ndim=1] theta = np.linspace(M_PI_2, M_PI, n_theta)

    # move theta off of the singularity at pi/2
    theta[0] += 1e-6

    s_val = s[n_layers-1]

    cdef:
        np.ndarray[FLOAT_TYPE_t, ndim=2] Theta = np.tile(theta, (n_phi, 1)).T

        # Component of wind along the phi direction (azimuth)
        np.ndarray[FLOAT_TYPE_t, ndim=2] u0 = np.tile(u[n_layers-1, :], (n_theta, 1))

        # ray parameter
        np.ndarray[FLOAT_TYPE_t, ndim=2] p = s_val*np.sin(Theta)/(1 + s_val*u0*np.sin(Theta))

        # Transformed x and y
        np.ndarray[FLOAT_TYPE_t, ndim=2] X = np.zeros((n_theta, n_phi))
        np.ndarray[FLOAT_TYPE_t, ndim=2] Y = np.zeros((n_theta, n_phi))

        # Transformed wind componenets
        np.ndarray[FLOAT_TYPE_t, ndim=2] U = np.empty((n_theta, n_phi))
        np.ndarray[FLOAT_TYPE_t, ndim=2] V = np.empty((n_theta, n_phi))

        #Travel time
        float t_arrival = 0

        #azimuth angle
        float azimuth = 0

        #takeoff angle
        float takeoff = 0

    # ignore negative roots

    np.seterr(divide='ignore', invalid='ignore')

    ### Scan Loop ###
    while not found:
        
        a = np.cos(Phi)
        b = np.sin(Phi)

        for i in range(n_layers - 1):

            s2 = s[i]**2
            delz = z[i + 1] - z[i]

            # Winds Enabled
            if wind:
                # clear old variables
                
                # Wind transformation variables
                U = np.tile(u[i, :], (n_theta, 1))
                V = np.tile(v[i, :], (n_theta, 1))
                
                p2 = p/(1 - p*U)
                # This term produces nans
                A = delz/np.sqrt(s2 - p2**2)

                # Equation (10)
                X += (p2 + s2*U)*A

                # Equation (11)
                Y += s2*V*A

                # Calculate true destination positions (transform back)
                #0.0016s


            # Winds Disabled
            else:

                # Equation (3)
                X += p*(delz)/(np.sqrt(s2 - p**2))

                # Calculate true destination positions (transform back)


        # Compare these destinations with the desired destination, all imaginary values are "turned rays" and are ignored
        E = np.sqrt(((a*X - b*Y - dx)**2 + (b*X + a*Y - dy)**2)) 

        # Ignore all nan slices - not interested in these points
        # Find the position of the smallest value
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            k, l = np.where(E == np.nanmin(E))

        # Check for all nan error function
        if k.shape == (0, ):
            # As handled in original Supracenter
            return np.array([np.nan, np.nan, np.nan])

        # If there are mulitple, take one closest to phi (in the middle)
        if len(k > 1):
            k, l = k[len(k)//2], l[len(l)//2]


        # If the error is within tolerance
        if E[k, l] < tol or dtheta < precision or dphi < precision:

            # pass on the azimuth & ray parameter information for use in traveltime calculation
            found = True

        else:
            ### FAST PART ###
            
            # reduce evenly in both directions
            n_phi = n_theta

            # General Case: central take off angle is between 0 & 90 degrees
            if ((theta[k] != M_PI_2) and (theta[k] != M_PI)):

                # Respace net around best value
                phi = np.linspace(phi[l] - dphi, phi[l] + dphi, n_phi)    
                dphi = 2*dphi/n_phi

                # Check: theta must be > 90 degrees
                if theta[k] - dtheta < M_PI_2:
                    theta = np.linspace(M_PI_2, theta[k] + 2*dtheta, n_theta)
                else: 
                    theta = np.linspace(theta[k] - dtheta, theta[k] + dtheta, n_theta) 
                dtheta = 2*dtheta/n_theta  

            # Case: central takeoff angle is at 180 degrees (vertical)
            elif (theta[k] == M_PI):

                # Respace net around best value
                # Higher accuracy in n_phi helps here
                phi = np.linspace(0, 2*M_PI - dphi, n_phi) 
                dphi = dphi/n_phi
                
                theta = np.linspace(theta[k] - dtheta, theta[k], n_theta)                      
                dtheta = dtheta/n_theta

            # Case: central takeoff angle is at 90 degrees (horizontal)
            elif (theta[k] == M_PI_2):

                # Respace net around best value
                phi = np.linspace(phi[l] - dphi, phi[l] + dphi, n_phi)         
                dphi = 2*dphi/n_phi  

                theta = np.linspace(M_PI_2, theta[k] + 2*dtheta, n_theta) 
                dtheta = dtheta/n_theta/2
            
            # Update values, and try again
            u = nwDir(z_profile[:, 2], z_profile[:, 3], phi)
            v = nwDir(z_profile[:, 2], z_profile[:, 3], phi + M_PI_2)

            n_theta = len(theta)
            n_phi = len(phi)

            # redefine variables
            Phi = np.tile(phi, (n_theta, 1))
            Theta = np.tile(theta, (n_phi, 1)).T
            u0 = np.tile(u[0, :], (n_theta, 1))


            ### END FAST PART ###
            #####

        # The minimum becomes the center of the next net with borders one spacing away

        X = np.zeros((n_theta, n_phi))
        Y = np.zeros((n_theta, n_phi))

        p = s_val*np.sin(Theta)/(1 + u0*np.sin(Theta)*s_val)



    ######################
    ### 0.0033s
    
    
    # Final solution for initial azimuth angle

    # Rotate coordinate system 90 degrees CCW
    # Flip coordinate system horizontally
    azimuth = (450 - phi[l]*180/M_PI)%360

    # Final solution for intial takeoff angle
    takeoff = (theta[k]*180/M_PI)%360

    p1 = p[k, l]
    p2 = p[k, l]**2
    # Find sum of travel times between layers (z)
    for i in range(n_layers - 1):

        s2 = s[i]**2
        # Equation (9)
        t_arrival += (s2/np.sqrt(s2 - p2/(1 - p1*u[i, l])**2))*(z[i + 1] - z[i])
    
    ##########################
    return np.array([t_arrival, azimuth, takeoff])

