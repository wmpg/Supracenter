import cython
cimport cython
import warnings

import time

import numpy as np
cimport numpy as np

FLOAT_TYPE = np.float64
ctypedef np.float64_t FLOAT_TYPE_t

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef float appatan(float z):
    
    # 0.29 deg error
    # cdef:
    #     float n1 = 0.97239411
    #     float n2 = -0.19194795
    # return (n1 + n2 * z * z) * z

    cdef:
        float n1 = -3.10715
        float n2 = 9.99042
        float a = z*z
    return 0.1*z*((a + n1)*a + n2)

    # cdef:
    #     float n1 = 0.2447
    #     float n2 = 0.0663
    #     float n3 = M_PI_2/2

    # if z >= 0.0:
    #     return n3*z - z*(z - 1)*(n1 + n2*z)
    # else:
    #     return n3*z + z*(z + 1)*(n1 - n2*z)

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef float appatan2(float y, float x):
    
    cdef:
        float z = y/x

    if x > 0.0:
        return appatan(z)
    elif x < 0.0:
        if y >= 0.0:
            return appatan(z) + np.pi
        else:
            return appatan(z) - np.pi
    else:
        if y > 0:
            return np.pi/2
        elif y < 0:
            return -np.pi/2
        else:
            #undefined
            return 0.0

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

    return A.T*B + C.T*D

@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
@cython.nonecheck(False)
cpdef np.ndarray[FLOAT_TYPE_t, ndim=1] cyscan(np.ndarray[FLOAT_TYPE_t, ndim=1] supra_pos, np.ndarray[FLOAT_TYPE_t, ndim=1] detec_pos, 
    np.ndarray[FLOAT_TYPE_t, ndim=2] z_profile, int wind=1, int n_theta=90, int n_phi=90, float h_tol=1e-5, float v_tol=1000):


    cdef:

        float pi_2 = np.pi/2
        int found = 0

        float dtheta = np.pi/n_theta
        
        float dphi = pi_2/n_phi   

        float dx = detec_pos[0] - supra_pos[0]
        float dy = detec_pos[1] - supra_pos[1]

        # Approximate atan2 function, accurate within 0.29 deg
        float azth = appatan2(dy, dx)

        int n_layers = len(z_profile)

        np.ndarray[FLOAT_TYPE_t, ndim=1] s = 1.0/z_profile[0:n_layers, 1]

        np.ndarray[FLOAT_TYPE_t, ndim=1] z  = z_profile[0:n_layers, 0]

        np.ndarray[FLOAT_TYPE_t, ndim=1] phi = np.linspace(azth-pi_2, azth+pi_2, n_phi)
        np.ndarray[FLOAT_TYPE_t, ndim=2] Phi = np.tile(phi, (n_theta, 1))
    
        np.ndarray[FLOAT_TYPE_t, ndim=2] u = nwDir(z_profile[:, 2], z_profile[:, 3], phi)
        np.ndarray[FLOAT_TYPE_t, ndim=2] v = nwDir(z_profile[:, 2], z_profile[:, 3], phi+pi_2)

        np.ndarray[FLOAT_TYPE_t, ndim=1] theta = np.linspace(pi_2 + 1e-6, np.pi, n_theta)

        float s_val = s[n_layers-1]

        np.ndarray[FLOAT_TYPE_t, ndim=2] Theta = np.tile(theta, (n_phi, 1)).T

        # Component of wind along the phi direction (azimuth)
        np.ndarray[FLOAT_TYPE_t, ndim=2] u0 = np.tile(u[n_layers - 1, :], (n_theta, 1))

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
        last_z = 0
        
        for i in range(n_layers - 1):
              
            s2 = s[i]**2
            delz = z[i + 1] - z[i]

            # Winds Enabled
            if wind:
                # clear old variables
                t1 = time.time()
                # Wind transformation variables
                U = np.tile(u[i, :], (n_theta, 1))
                V = np.tile(v[i, :], (n_theta, 1))
                t2 = time.time()
                print('TILE: {:E}'.format(t2-t1))

                t1 = time.time()
                p2 = p/(1 - p*U)
                t2 = time.time()
                print('P2: {:E}'.format(t2-t1))
                # This term produces nans
                t1 = time.time()
                A = delz/np.sqrt(s2 - p2**2)
                t2 = time.time()
                print('A: {:E}'.format(t2-t1))
                # Equation (10)
                t1 = time.time()
                X += (p2 + s2*U)*A

                # Equation (11)
                Y += s2*V*A
                t2 = time.time()
                print('XY: {:E}'.format(t2-t1))
                last_z = i + 1

            # Winds Disabled
            else:

                # Equation (3)
                X += p*(delz)/(np.sqrt(s2 - p**2))
                last_z = i + 1
                # Calculate true destination positions (transform back)


        # Compare these destinations with the desired destination, all imaginary values are "turned rays" and are ignored
        if wind:
            x_err = a*X - b*Y - dx
            y_err = b*X + a*Y - dy
        else:
            x_err = a*X - dx
            y_err = b*X - dy

        z_err = z[n_layers - last_z - 1] - detec_pos[2]

        E = np.sqrt((x_err*x_err + y_err*y_err + z_err*z_err))

        # Find the position of the smallest value
        
        k, l = np.where(E == np.nanmin(E))

        # Check for all nan error function
        if k.shape == (0, ):

            return np.array([np.nan, np.nan, np.nan, np.nan])

        # If there are mulitple, take one closest to phi (in the middle)
        if len(k > 1):
            k, l = k[len(k)//2], l[len(l)//2]

        if E[k, l] < v_tol:
            found = True

        elif dtheta < h_tol or dphi < h_tol:

            return np.array([np.nan, np.nan, np.nan, np.nan])

        else:
            ### FAST PART ###

            # reduce evenly in both directions
            n_phi = n_theta

            # General Case: central take off angle is between 0 & 90 degrees
            if ((theta[k] != pi_2) and (theta[k] != np.pi)):

                # Respace net around best value
                phi = np.linspace(phi[l] - dphi, phi[l] + dphi, n_phi)    
                dphi = 2*dphi/n_phi

                # Check: theta must be > 90 degrees
                if theta[k] - dtheta < pi_2:
                    theta = np.linspace(pi_2, theta[k] + 2*dtheta, n_theta)
                else: 
                    theta = np.linspace(theta[k] - dtheta, theta[k] + dtheta, n_theta) 
                dtheta = 2*dtheta/n_theta  

            # Case: central takeoff angle is at 180 degrees (vertical)
            elif (theta[k] == np.pi):

                # Respace net around best value
                # Higher accuracy in n_phi helps here
                phi = np.linspace(0, 2*np.pi - dphi, n_phi) 
                dphi = dphi/n_phi
                
                theta = np.linspace(theta[k] - dtheta, theta[k], n_theta)                      
                dtheta = dtheta/n_theta

            # Case: central takeoff angle is at 90 degrees (horizontal)
            elif (theta[k] == pi_2):

                # Respace net around best value
                phi = np.linspace(phi[l] - dphi, phi[l] + dphi, n_phi)         
                dphi = 2*dphi/n_phi  

                theta = np.linspace(pi_2, theta[k] + 2*dtheta, n_theta) 
                dtheta = dtheta/n_theta/2
            
            # Update values, and try again
            u = nwDir(z_profile[:, 2], z_profile[:, 3], phi)
            v = nwDir(z_profile[:, 2], z_profile[:, 3], phi + pi_2)

            # redefine variables
            Phi = np.tile(phi, (n_theta, 1))
            Theta = np.tile(theta, (n_phi, 1)).T
            u0 = np.tile(u[0, :], (n_theta, 1))


        X = np.zeros((n_theta, n_phi))
        Y = np.zeros((n_theta, n_phi))

        p = s_val*np.sin(Theta)/(1 + u0*np.sin(Theta)*s_val)



    ######################
    ### 0.0033s

    # Final solution for initial azimuth angle

    # Rotate coordinate system 90 degrees CCW
    # Flip coordinate system horizontally
    azimuth = (450 - np.degrees(phi[l]))%360

    # Final solution for intial takeoff angle
    takeoff = np.degrees(theta[k])

    p1 = p[k, l]
    p2 = p[k, l]**2
    # Find sum of travel times between layers (z)
    for i in range(n_layers - 1):# - n_layers + last_z):

        s2 = s[i]**2
        # Equation (9)
        t_arrival += (s2/np.sqrt(s2 - p2/(1 - p1*u[i, l])**2))*(z[i + 1] - z[i])
    
    ##########################
    return np.array([t_arrival, azimuth, takeoff, E[k, l]])
