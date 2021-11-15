import numpy as np

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})


def anglescan(S, phi, theta, z_profile, wind=True, debug=True, trace=False):
    # switched positions (Jun 2019)

    # This function should be called for every station
    # Original Author: Wayne Edwards
    # Theta is vertical
    # phi is horizontal

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
    phi = (phi - 90)%360

    # Flip coordinate system horizontally
    phi = (360 - phi)%360

    phi = np.radians(phi)
    theta = np.radians(theta)
    # Switch to turn off winds
    if not wind:
        z_profile[:, 2] = 0
        z_profile[:, 1] = 330

    # The number of layers in the integration region
    n_layers = len(z_profile)

    # Slowness, as defined in SUPRACENTER on pg 35, s = 1/c
    s = 1.0/z_profile[0:n_layers, 1]

    # Elevation for that layer
    z  = z_profile[0:n_layers, 0]

    # Component of wind vector in the direction of phi and phi + pi/2 respectively

    u = z_profile[:, 2]*np.sin(z_profile[:, 3])*np.cos(phi) + z_profile[:, 2]*np.cos(z_profile[:, 3])*np.sin(phi)
    v = z_profile[:, 2]*np.sin(z_profile[:, 3])*np.cos(phi+np.pi/2) + z_profile[:, 2]*np.cos(z_profile[:, 3])*np.sin(phi+np.pi/2)
    
    s_val = s[n_layers-1]

    # ray parameter
    p = s_val*np.sin(theta)/(1 + s_val*u[n_layers - 1]*np.sin(theta))

    X = 0
    Y = 0
    #Travel time
    t_arrival = 0

    if trace:
        T = []
        T.append([S[0], S[1], S[2], t_arrival])

    # ignore negative roots
    np.seterr(divide='ignore', invalid='ignore')

    ### Scan Loop ###
    a, b = np.cos(phi), np.sin(phi)
    last_z = 0 
    for i in range(n_layers - 1, 0, -1):

        s2 = s[i]**2
        delz = z[i] - z[i - 1]

        # Winds Enabled
        if wind:
            # clear old variables

            # Wind transformation variables
            U = u[i]
            V = v[i]
            
            p2 = p/(1 - p*U)
            # This term produces nans
            A = delz/np.sqrt(s2 - p2**2)

            if np.isnan(A).all():

                if debug:
                    print("ANGLESCAN ERROR: All NaNs - rays reflect upwards")

                if trace:
                    return np.array([[np.nan, np.nan, np.nan, np.nan]])
                else:
                    return np.array([np.nan, np.nan, np.nan, np.nan])

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
        last_z = i - 1


        t_arrival += (s2/np.sqrt(s2 - p**2/(1 - p*u[i])**2))*delz
        
        if trace:
            T.append([S[0] + (a*X - b*Y), S[1] + (b*X + a*Y), z[last_z], t_arrival])
        # if v_tol is not None and h_tol is not None:
        #     dh = z[last_z] - target[2]
        #     dx = np.sqrt((S[0] + (a*X - b*Y) - target[0])**2 + (S[1] + (b*X + a*Y) - target[1])**2)
        #     if dh <= v_tol and dx <= h_tol:
        #         t_arrival += np.sqrt(dh**2 + dx**2)/310

        # Compare these destinations with the desired destination, all imaginary values are "turned rays" and are ignored
    # E = np.sqrt(((a*X - b*Y)**2 + (b*X + a*Y)**2 + (z[n_layers - last_z - 1])**2)) 

    D = [S[0] + (a*X - b*Y), S[1] + (b*X + a*Y), z[last_z], t_arrival]

    ##########################
    if trace:
        return np.array(D), np.array(T)
    else:
        return np.array(D)

if __name__ == '__main__':
    S = np.array([0, 0, 1000])

    #takeoff
    theta = 135

    #azimuth
    phi = 0

    z_profile = np.array([[    0, 330, 0, 0],
                          [500, 330, 0, 0],
                          [1000, 330, 0, 0]])
    D = anglescan(S, phi, theta, z_profile, trace=True)
    print(D)