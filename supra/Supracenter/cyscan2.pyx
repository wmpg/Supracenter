# Cyscan V 2.0
# Based off of Wayne Edwards' ray-tracing algorithm (2003)               
# Finds path between two locations with an atmospheric profile in between 
###########################################################################

###########
# IMPORTS #
###########
import cython
cimport cython
import warnings

import time

import numpy as np
cimport numpy as np
from libc.math cimport sqrt, M_PI, M_PI_2
 
# Define cython types for numpy arrays
FLOAT_TYPE = np.float64
ctypedef np.float64_t FLOAT_TYPE_t


############################
# FAST APPROXIMATE ARCTAN2 #
############################
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef FLOAT_TYPE_t appatan(float z):

    cdef:
        FLOAT_TYPE_t n1 = -3.10715
        FLOAT_TYPE_t n2 = 9.99042
        FLOAT_TYPE_t a = z*z

    # Approximation of arctan 
    return 0.1*z*((a + n1)*a + n2)

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef FLOAT_TYPE_t appatan2(float y, float x):
    '''
    Return arctan with the correct sign
    '''
    cdef:
        FLOAT_TYPE_t z = y/x

    if x > 0.0:
        return appatan(z)
    elif x < 0.0:
        if y >= 0.0:
            return appatan(z) + np.pi
        else:
            return appatan(z) - np.pi
    else:
        if y > 0:
            return 1.5707963
        elif y < 0:
            return -1.5707963
        else:
            #undefined
            return 0.0

################
# WIND VECTORS #
################
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef np.ndarray[FLOAT_TYPE_t, ndim=1] nwDir(np.ndarray[FLOAT_TYPE_t, ndim=1] w, np.ndarray[FLOAT_TYPE_t, ndim=1] wd, np.ndarray[FLOAT_TYPE_t, ndim=1] phi):
    # Author: Wayne Edwards

    """ Calculates the component of the wind vector in the direction of phi for each layer in the atmosphere

    Arguments:
        w: [ndarray] Wind Magnitude for each level in the atmosphere
        wd: [ndarray] Wind direction (angle). Direction (from the north) the wind is blowing to (Radians)
        phi: [ndarray] Azimuthal angles (NDE)

    Returns:
        x: [ndarray] Component of the wind in the direction of phi
    """

    A = np.expand_dims(w*np.sin(wd), axis=0)
    B = np.expand_dims( np.cos(phi), axis=0)
    C = np.expand_dims(w*np.cos(wd), axis=0)
    D = np.expand_dims( np.sin(phi), axis=0)

    return np.matmul(A.T, B) + np.matmul(C.T, D)


##########
# CYSCAN #
##########
@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef np.ndarray[FLOAT_TYPE_t, ndim=1] cyscan(np.ndarray[FLOAT_TYPE_t, ndim=1] supra_pos, np.ndarray[FLOAT_TYPE_t, ndim=1] detec_pos, 
    np.ndarray[FLOAT_TYPE_t, ndim=2] z_profile, wind=True, int n_theta=90, int n_phi=90, float h_tol=330, float v_tol=1000, debug=False, trace=False):

    # Original Author: Wayne Edwards

    """ Finds an optimal ray from the source of the Supracenter to the detector, by making a guess, 
        and checking for the angle of minimum error. This is repeated with better guesses by taking angles around the best angle of the previous step.

    Arguments:  
        supra_pos: [array] 3-D local coordinates of the source of the sound. Given as: [x, y, Elevation]
        detec_pos: [array] 3-D local coordinates of the detector. Given as: [x, y, Elevation]
        zProfile: [array] The given atmospheric profile between the source and the detector. Given as: [Height, Temperature, Wind Speed, Wind Direction] for different values of height
        wind: [int] 0 - disable winds, 1 - enable winds. Temperature is still used.
        n_theta, n_phi: [int] angle grid spacing of the takeoff, azimuth angles
        h_tol, v_tol: [float] minimum horizontal and vertical distance away from reciever to be considered an arrival

    Returns:    
        t_arrival: [float] Direct arrival time between source and detector
        Azimuth: [float] The initial azimuthal angle for the source to the detector in degrees
        Takeoff: [float] The initial takeoff angle from the source to the detector in degrees

    """

    # Switch to turn off winds
    if not wind:
        z_profile[:, 2] = 0

    found = False

    ########################
    # Initialize Variables #
    ########################
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
        float azth = appatan2(dy, dx)

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
        np.ndarray[FLOAT_TYPE_t, ndim=2] u0 = np.tile(u[n_layers - 1, :], (n_theta, 1))

        # ray parameter
        np.ndarray[FLOAT_TYPE_t, ndim=2] p = s_val*np.sin(Theta)/(1 + s_val*u0*np.sin(Theta))

        # Transformed x and y
        np.ndarray[FLOAT_TYPE_t, ndim=2] X = np.zeros((n_theta, n_phi))
        np.ndarray[FLOAT_TYPE_t, ndim=2] Y = np.zeros((n_theta, n_phi))

        #Travel time
        float t_arrival = 0

        #azimuth angle
        float azimuth = 0

        #takeoff angle
        float takeoff = 0

    # ignore negative roots
    last_error = 1e20
    np.seterr(divide='ignore', invalid='ignore')



#############
# Scan Loop #
#############

    while not found:
        if trace:
            trace_list = []
            trace_list.append([supra_pos[0], supra_pos[1], supra_pos[2]])
        count = 0
        a, b = np.cos(Phi), np.sin(Phi)
        last_z = 0
        
        if n_layers <= 1:
            if debug:
                print('CYSCAN ERROR: No layers')
            if trace:
                return np.array([np.nan, np.nan, np.nan, np.nan, trace_list])
            return np.array([np.nan, np.nan, np.nan, np.nan])

        for i in range(n_layers - 1):

            s2 = s[i]**2
            delz = z[i + 1] - z[i]

            # Winds Enabled
            if wind:
                # clear old variables
                
                p2 = p/(1 - p*u[i, :])

                # This term produces nans
                A = delz/np.sqrt(s2 - p2**2)

                # If this is true, the ray has reflected upward
                if np.isnan(A).all():
                    break

                # Equation (10)
                X += (p2 + s2*u[i, :])*A

                # Equation (11)
                Y += s2*v[i, :]*A

                # Calculate true destination positions (transform back)
                last_z = i + 1

            # Winds Disabled
            else:

                # Equation (3)
                X += p*(delz)/(np.sqrt(s2 - p**2))
                last_z = i + 1
            

            # Calculate true destination positions (transform back)
            horizontal_error = np.sqrt(((a*X - b*Y - dx)**2 + (b*X + a*Y - dy)**2))
            vertical_error = z[n_layers - last_z - 1] - detec_pos[2]

            if trace:
                k, l = np.where(horizontal_error == np.nanmin(horizontal_error))
                trace_list.append([(a*X - b*Y)[k, l], (b*X + a*Y)[k, l], z[n_layers - 1 - last_z]])

            # If the ray is close enough to the station, stop calculating
            # done in order to reduce NaNs -> Theoretically it should go further if it can,
            # but checking for that takes time

            if np.nanmin(horizontal_error) < h_tol and vertical_error < v_tol:
                k, l = np.where(horizontal_error == np.nanmin(horizontal_error))
                found = True
                break

        # Compare these destinations with the desired destination, all imaginary values are "turned rays" and are ignored
        E = np.sqrt(((a*X - b*Y - dx)**2 + (b*X + a*Y - dy)**2 + (z[n_layers - last_z - 1] - detec_pos[2])**2)) 


        # Fast way to check for all nan
        if np.isnan(E).all():
            # RETURN 1 of 3:
            # Failure
            # Caused by ray reflecting upwards
            if debug:
                print('CYSCAN ERROR: All NaNs - Rays reflect upwards')
            if trace:
                return np.array([np.nan, np.nan, np.nan, np.nan, trace_list])
            return np.array([np.nan, np.nan, np.nan, np.nan])

        # Ignore all nan slices - not interested in these points
        # Find the position of the smallest value

        k, l = np.unravel_index(np.nanargmin(E), (n_theta, n_phi))

        new_error = E[k, l]

        if dphi < 1e-10:

            # RETURN 2 of 3:
            # Failure
            # Caused by getting stuck in a loop, cannot find solution
            if debug:
                print('CYSCAN ERROR: Cannot find solution - Angle precision less than tolerance')
            if trace:
                return np.array([np.nan, np.nan, np.nan, np.nan, trace_list])
            return np.array([np.nan, np.nan, np.nan, np.nan])

        if horizontal_error[k, l] < h_tol and vertical_error < v_tol:
        # if E[k, l] < v_tol or dtheta < h_tol or dphi < h_tol:

            # pass on the azimuth & ray parameter information for use in traveltime calculation
            found = True


        else:
            last_error = E[k, l]
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


        # The minimum becomes the center of the next net with borders one spacing away

        X = np.zeros((n_theta, n_phi))
        Y = np.zeros((n_theta, n_phi))

        p = s_val*np.sin(Theta)/(1 + u0*np.sin(Theta)*s_val)

#################
# Return Values #
#################
    
    azimuth = (450 - phi[l]*180/M_PI)%360

    # Final solution for intial takeoff angle
    takeoff = (theta[k]*180/M_PI)%360

    p1 = p[k, l]
    p2 = p[k, l]**2

    for i in range(last_z):

        s2 = s[i]**2
        # Equation (9)
        t_arrival += (s2/np.sqrt(s2 - p2/(1 - p1*u[i, l])**2))*(z[i + 1] - z[i])

    # Add in rest of travel time as an approximation
    t_arrival += np.sqrt(vertical_error**2 + horizontal_error[k, l]**2)/310

    # Return 3 of 3:
    # Success
    if trace:
        return np.array([t_arrival, azimuth, takeoff, E[k, l], trace_list])
    return np.array([t_arrival, azimuth, takeoff, E[k, l]])
