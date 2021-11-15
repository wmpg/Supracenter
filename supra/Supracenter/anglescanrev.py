import numpy as np

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})


def anglescanrev(S, phi, theta, z_profile, wind=True, trace=False, fix_phi=True):
    #This code is basically cheating, it flips the vertical scale, ray-traces, then flips back.
    # Imagine that it shoots the ray into a mirror world into the ground

    z_profile = np.flipud(z_profile)
    # z_profile[:, 3] = z_profile[:, 3] + np.pi

    if fix_phi:
        # Azimuths and Wind directions are measured as angles from north, and increasing clockwise to the East
        phi = (phi - 90)%360

        # Flip coordinate system horizontally
        phi = (360 - phi)%360

    # Theta is given as angle from ground up to beam shot upward. 
    # Here we use a beam shot downward into the ground and flip it up
    theta = 180 - theta

    phi = np.radians(phi)
    theta = np.radians(theta)
    # Switch to turn off winds
    if not wind:
        z_profile[:, 2] = 0
        z_profile[:, 1] = 330

    if trace:
        T = []
        T.append([S[0], S[1], S[2], 0.0])

    # The number of layers in the integration region
    n_layers = len(z_profile)

    # Slowness, as defined in SUPRACENTER on pg 35, s = 1/c
    s = 1.0/z_profile[0:n_layers, 1]

    # Elevation for that layer
    z  = -z_profile[0:n_layers, 0]

    # Component of wind vector in the direction of phi and phi + pi/2 respectively
    # Backwards Wind
    u = z_profile[:, 2]*np.sin(z_profile[:, 3])*np.cos(phi) + z_profile[:, 2]*np.cos(z_profile[:, 3])*np.sin(phi)
    v = z_profile[:, 2]*np.sin(z_profile[:, 3])*np.cos(phi+np.pi/2) + z_profile[:, 2]*np.cos(z_profile[:, 3])*np.sin(phi+np.pi/2)
    
    s_val = s[n_layers-1]

    # ray parameter
    p = s_val*np.sin(theta)/(1 + s_val*u[n_layers - 1]*np.sin(theta))

    X = 0
    Y = 0
    #Travel time
    t_arrival = 0

    # ignore negative roots
    np.seterr(divide='ignore', invalid='ignore')


    ### Scan Loop ###
    a, b = np.cos(phi), np.sin(phi)
    last_z = 0 
    for i in range(n_layers - 1, 0, -1):

        s2 = s[i]**2
        delz = z[i] - z[i-1]

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

                if trace:
                    return np.array([[np.nan, np.nan, np.nan, np.nan]])
                else:
                    return np.array([np.nan, np.nan, np.nan, np.nan])

            # Equation (10)
            dx = (p2 + s2*U)*A
            X += dx

            # Equation (11)
            dy = s2*V*A
            Y += dy

            # Calculate true destination positions (transform back)
            #0.0016s
            

        else:

            # Equation (3)
            dx = p*(delz)/(np.sqrt(s2 - p**2))
            X += dx
            Y = 0
            dy = 0

        last_z = i - 1
        t_arrival += s2/np.sqrt(s2 - p**2/(1 - p*u[i-1])**2)*delz
        if trace:
            T.append([S[0] + (a*X - b*Y), S[1] + (b*X + a*Y), np.abs(z[last_z]), t_arrival])
        


        

    if not trace:

        D = [S[0] + (a*X - b*Y), S[1] + (b*X + a*Y), np.abs(z[last_z]), t_arrival]
        return np.array(D)
    else:
        return np.array(T)

if __name__ == '__main__':
    S = np.array([0, 0, 33000])

    #takeoff
    theta = 175

    #azimuth
    phi = 0

    z_profile = np.array([[    0, 330, 0, 0],
                          [11500, 330, 0, 0],
                          [33000, 330, 0, 0]])
    D = anglescan(S, phi, theta, z_profile)
    print(D)