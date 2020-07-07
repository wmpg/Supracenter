
import numpy as np
import matplotlib.pyplot as plt
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})
from supra.Supracenter.cyscan2 import cyscan
from supra.Supracenter.anglescan import anglescan
from supra.Supracenter.cyzInteg import zInteg

def eigenError(S, D, z_profile, n_angle=90):
    ''' Used to check the error, and verify that the eigenray solution works properly'''

    z_profile = zInteg(S[2], D[2], z_profile)

    f_time, frag_azimuth, frag_takeoff, frag_err = cyscan(S, D, z_profile, wind=True, \
                    n_theta=n_angle, n_phi=n_angle, h_tol=1e-15, v_tol=330)

    D_new = anglescan(S, frag_azimuth, frag_takeoff, z_profile, wind=True)

    d_loc_err = np.sqrt((D[0]-D_new[0])**2 + (D[1]-D_new[1])**2)
    d_time_err = D_new[3] - f_time
    d_err_err = d_loc_err - frag_err


    print("Spatial Error: ", d_loc_err)
    print("Temporal Error: ", d_time_err)
    print("Error Error: ", d_err_err)

def eigenAngleError(S, az, tf, z_profile, n_angle=90):

    D = anglescan(S, az, tf, z_profile, wind=True)

    z_profile = zInteg(S[2], D[2], z_profile)

    f_time, frag_azimuth, frag_takeoff, frag_err = cyscan(S, D[:3], z_profile, wind=True, \
                n_theta=n_angle, n_phi=n_angle, h_tol=1e-15, v_tol=330)

    tf_error = tf - frag_takeoff
    az_error = az - frag_azimuth
    time_err = D[3] - f_time

    print("Takeoff Error [deg]: ", tf_error)
    print("Azimuth Error [deg]: ", az_error)
    print("Temporal Error: ", time_err)
    print("Spatial Error: ", frag_err)

def eigenConsistancy(S, D, az, tf, z_profile, n_angle=90, h_tol=None, v_tol=None):
    '''
    Checks if a given source-detector, atmosphere, initial launch angle, solution is consistant with itself
    by running through both anglescan and cyscan 
    '''
    n_angle = 360
    D_new = anglescan(S, az, tf, z_profile, wind=True, debug=True)

    if np.isnan(D_new[0]):
        print('Inconsistant - bad anglescan')
        return 0

    T_new, az_new, tf_new, err = cyscan(S, D_new, z_profile, wind=True, \
                n_theta=n_angle, n_phi=n_angle, h_tol=330, v_tol=330, debug=True)

    d_loc_err = np.sqrt((D[0]-D_new[0])**2 + (D[1]-D_new[1])**2)
    d_time_err = D_new[3] - T_new
    tf_error = tf - tf_new
    az_error = az - az_new
    d_err_err = d_loc_err - err

    if np.isnan(T_new):
        print('Inconsistant - nan solution')
        return 0

    max_error = np.sqrt(h_tol**2 + v_tol**2)
    
    if not (d_loc_err < max_error and d_time_err < max_error/310 and tf_error < 1 and az_error < 1):
        print("Inconsistant Solution at height {:.2f} km".format(S[2]/1000))
        print('Spatial Error: ', d_loc_err)
        print('Temporal Error: ', d_time_err)
        print("Takeoff Error [deg]: ", tf_error)
        print("Azimuth Error [deg]: ", az_error)
        print("Error Error: ", d_err_err)
        return 0

    return 1