"""Finds component of wind normal to direction of travel"""


import numpy as np



def nwDir(w, wd, phi):
    # Author: Wayne Edwards

    """ Calculates the component of the wind vector in the direction of phi for each layer in the atmosphere

    Arguments:
        w: [ndarray] Wind Magnitude for each level in the atmosphere
        wd: [ndarray] Wind direction (angle). Direction (from the north) the wind is blowing to (Radians)
        phi: [ndarray] Azimuthal angles (NDE)

    Returns:
        x: [ndarray] Component of the wind in the direction of phi
    """


    m = len(w)
    n = len(phi)

    return np.tile(w*np.sin(wd), (n, 1)).T*np.tile(np.cos(phi), (m, 1)) + np.tile(w*np.cos(wd), (n, 1)).T*np.tile(np.sin(phi), (m, 1))