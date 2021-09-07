import numpy as np
from wmpl.Utils.TrajConversions import latLonAlt2ECEF, ecef2LatLonAlt
from supra.Utils.Classes import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def findLayer(h, z_profile):

    found_layer = None
    for ll, layer in enumerate(z_profile):
        if layer[0] <= h:
            found_layer = layer
            found_index = ll
        else:
            break

    return found_layer, ll

def propTime(pos, z_profile, dz, u, v, p):

    # x is East, y is North

    current_height = pos.elev + dz

    layer, i = findLayer(current_height, z_profile)

    if layer is None:
        return None

    s2 = (layer[1])**(-2)
    U = u[i]
    V = v[i]

    p2 = p/(1 - p*U)

    A = dz/np.sqrt(s2 - p2**2)

    if np.isnan(A):
        return None

    # Equation (10)
    dx = (p2 + s2*U)*A

    # Equation (11)
    dy = s2*V*A

    dt = -s2/np.sqrt(s2 - p**2/(1 - p*u[i])**2)*dz

    return np.array([dx, dy, dz, dt])


def anglescan(S, phi, theta, z_profile, wind=True, debug=True, trace=False, plot=False, dz=-50, last_z=0):
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


    u = z_profile[:, 2]*np.sin(z_profile[:, 3])*np.cos(phi) + z_profile[:, 2]*np.cos(z_profile[:, 3])*np.sin(phi)
    v = z_profile[:, 2]*np.sin(z_profile[:, 3])*np.cos(phi+np.pi/2) + z_profile[:, 2]*np.cos(z_profile[:, 3])*np.sin(phi+np.pi/2)
    
    s_val = 1/z_profile[-1, 1]

    # ray parameter
    p = s_val*np.sin(theta)/(1 + s_val*u[-1]*np.sin(theta))
 
    S_init = S

    t_arrival = 0
    if trace:
        T = []
        T.append([S.lat, S.lon, S.elev, t_arrival])

    # ignore negative roots
    np.seterr(divide='ignore', invalid='ignore')

    done = False

    while not done:

        S_ref = latLonAlt2ECEF(S.lat_r, S.lon_r, S.elev)

        diff = propTime(S, z_profile, dz, u, v, p)

        if diff is None:
            return None
        else:
            x, y, z, t = diff

        t_arrival += t

    
        new_pos = [S_ref[0] + x, S_ref[1] + y, S_ref[2] + z]


        new_geo_pos = ecef2LatLonAlt(new_pos[0], new_pos[1], new_pos[2])

        S = Position(np.degrees(new_geo_pos[0]), np.degrees(new_geo_pos[1]), new_geo_pos[2])

        if trace:
            T.append([S.lat, S.lon, S.elev, t_arrival])

        if S.elev <= last_z:
            done = True



    
    if trace and plot:
        tr = np.array(T)

        fig = plt.figure()
        ax = Axes3D(fig)
        ax.scatter(tr[:, 1], tr[:, 0], tr[:, 2], c='b')
        ax.plot(tr[:, 1], tr[:, 0], tr[:, 2], c='k')
        ax.scatter(S_init.lon, S_init.lat, S_init.elev, c='r', marker="*")
        ax.scatter(S.lon, S.lat, S.elev, c='g', marker="^")
        plt.show()


    D = [S.lat, S.lon, S.elev, t_arrival]


    
    ##########################
    if trace:
        return np.array(D), np.array(T)
    else:
        return np.array(D)

if __name__ == '__main__':
    S = Position(45, 45, 10000)

    #takeoff
    theta = 135

    #azimuth
    phi = 0

    z_profile = np.array([[   0.0, 330.0, 4.0, 0.0],
                          [1000.0, 330.0, 4.0, 0.0],
                          [2020.0, 330.0, 4.0, 0.0],
                          [3023.0, 350.0, 4.0, 0.0],
                          [4000.0, 350.0, 4.0, 0.0],
                          [5400.0, 330.0, 4.0, 0.0],
                          [6000.0, 330.0, 4.0, 0.0],
                          [8500.0, 330.0, 4.0, 0.0],
                          [8900.0, 330.0, 4.0, 90.0],
                          [10000.0, 330.0, 4.0, 0.0]])
    D = anglescan(S, phi, theta, z_profile, trace=True, plot=True)
    print(D)