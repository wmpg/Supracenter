import numpy as np


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def anglescan(S, phi, theta, z_profile, wind=True, debug=True, trace=False, plot=False):
    # Originally by Wayne Edwards (Supracenter)

    """ Ray-traces from a point given initial launch angles
    
    Arguments:  
        S: [list] [x, y, z] of initial launch point (Supracenter or Wave-Release point)
        phi: [float] initial azimuthal angle of launch [deg] with 0 deg being North and 90 deg being East
        theta: [float] initial takeoff angle of launch [deg] with 90 deg being horizontal and 180 deg being vertically down
        z_profile: [list] weather profile (n_layers * 4)
                [[heights (increasing order) [m], speed of sound [m/s], wind speed [m/s], wind direction [rad] (same angle definition as phi)],
                    ... ]
    Keyword Arguments:
        wind: [Boolean] if False sets all wind speeds to 0
        debug: [Boolean] if True outputs print messages of program status
        trace: [Boolean] if True returns (x, y, z, t) coordinates of the ray trace
        plot: [Boolean] if True plots the ray trace 


    Returns:    
        D: [list] (x, y, z, t) final position and travel time of the raytrace
        T: [list] returned if trace is set to True, (x, y, z, t) of all points along the ray-trace
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
        # z_profile[:, 1] = 330

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
    
    if trace and plot:
        tr = np.array(T)

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(tr[:, 0], tr[:, 1], tr[:, 2], c='b')
        ax.plot(tr[:, 0], tr[:, 1], tr[:, 2], c='k')
        ax.scatter(S[0], S[1], S[2], c='r', marker="*")
        ax.scatter(S[0] + (a*X - b*Y), S[1] + (b*X + a*Y), z[last_z], c='g', marker="^")
        plt.show()
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
    D = anglescan(S, phi, theta, z_profile, trace=True, plot=True)
    print(D)