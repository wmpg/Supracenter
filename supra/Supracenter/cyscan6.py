import numpy as np
import matplotlib.pyplot as plt

from supra.Utils.Classes import *
from supra.Supracenter.anglescan3 import anglescan
from mpl_toolkits.mplot3d import Axes3D
from wmpl.Utils.TrajConversions import latLonAlt2ECEF, ecef2LatLonAlt
from supra.Utils.pso import pso
import multiprocessing



def angleErr(x, *cyscan_inputs):
    S, z_profile, D, wind, debug, h_tol, v_tol = cyscan_inputs
    r = anglescan(S, x[0], x[1], z_profile, trace=False, debug=debug, wind=wind, last_z=1000)

    if r is None:
        return np.inf

    r_xyz = np.array(latLonAlt2ECEF(np.radians(r[0]), np.radians(r[1]), r[2]))
    D_xyz = np.array(latLonAlt2ECEF(D.lat_r, D.lon_r, D.elev))

    err = np.sqrt((D_xyz[0] - r_xyz[0])**2 + (D_xyz[1] - r_xyz[1])**2  + (D_xyz[2] - r_xyz[2])**2 )

    return err

SWARM_SIZE = 100
MAXITER = 100
PHIP = 0.5
PHIG = 0.5
OMEGA = 0.5
MINFUNC = 1e-8
MINSTEP = 1e-8


def cyscan(S, D, z_profile, trace=False, plot=False, particle_output=False, debug=False, wind=False, h_tol=330, v_tol=3000, \
        print_times=False):
    
    
    # phi, theta
    search_min = [0, 90]
    search_max = [360, 180]

    cyscan_inputs = [S, z_profile, D, wind, debug, h_tol, v_tol]

    if particle_output:
        f_opt, x_opt, f_particle, x_particle = pso(angleErr, search_min, search_max, \
            args=cyscan_inputs, processes=multiprocessing.cpu_count()-1, particle_output=True, swarmsize=SWARM_SIZE,\
                 maxiter=MAXITER, phip=PHIP, phig=PHIG, \
                 debug=False, omega=OMEGA, minfunc=MINFUNC, \
                 minstep=MINSTEP)

    else:
        f_opt, x_opt = pso(angleErr, search_min, search_max, \
            args=cyscan_inputs, processes=1, particle_output=False, swarmsize=SWARM_SIZE,\
                 maxiter=MAXITER, phip=PHIP, phig=PHIG, \
                 debug=False, omega=OMEGA, minfunc=MINFUNC, \
                 minstep=MINSTEP)


    if trace:
        r = anglescan(S, f_opt[0], f_opt[1], z_profile, trace=True)
        tr = np.array(r[1])

        if plot:
            
            fig = plt.figure()
            ax = Axes3D(fig)
            ax.scatter(S.lon, S.lat, S.elev, c='r', marker='*')
            ax.scatter(D.lon, D.lat, D.elev, c='g', marker='^')
            ax.scatter(tr[:, 0], tr[:, 1], tr[:, 2], c='b')
            ax.plot(tr[:, 0], tr[:, 1], tr[:, 2], c='k')

            # Plot all missed angles
            if particle_output:
                for particle in range(len(f_particle)):
                    r = anglescan(S, f_particle[particle][0], f_particle[particle][1], z_profile, trace=True)
                    tr = np.array(r[1])
                    ax.plot(tr[:, 0], tr[:, 1], tr[:, 2], alpha=0.3)

            plt.show()
    else:
        r = anglescan(S, f_opt[0], f_opt[1], z_profile, trace=False)


    if debug:
        print("Final Solution: {:.2f} {:.2f}".format(f_opt[0], f_opt[1]))
        print("Final Error:    {:.2f}".format(x_opt))
        if trace:
            phi_list = []
            the_list = []
            for i in range(len(tr[:, 2]) - 1):
                f_v = np.abs(tr[i+1, 2] - tr[i, 2])
                f_h = np.sqrt((tr[i+1, 0] - tr[i, 0])**2 + (tr[i+1, 1] - tr[i, 1])**2)
                f_x = tr[i+1, 0] - tr[i, 0]
                f_y = tr[i+1, 1] - tr[i, 1]
                tr_phi = np.degrees(np.arctan2(f_x, f_y))
                tr_theta = np.degrees(np.arctan2(f_v, f_h)) + 90
                phi_list.append(tr_phi)
                the_list.append(tr_theta)
                if i == 0:
                    # Note that the trace solution may not equal the actual solution if there are winds
                    print("Trace Solution: {:.2f} {:.2f}".format(tr_phi, tr_theta))
            print("Mean Angles: {:.2f} {:.2f}".format(np.nanmean(phi_list), np.nanmean(the_list)))
    
    if trace:
        x, y, z, T = r[0]
    else:
        x, y, z, T = r

    # if print_times:
    #     alltimes = []
    #     for particle in range(len(f_particle)):
    #         r = anglescan(S, f_particle[particle][0], f_particle[particle][1], z_profile, trace=False)
    #         x_t, y_t, z, T = r
    #         h_err = np.sqrt((D[0] - x)**2 + (D[1] - y)**2)
    #         v_err = np.sqrt(D[2] - z)
    #         if h_err <= h_tol and v_err <= v_tol:
    #             alltimes.append(T)

    #     if len(alltimes) == 0:
    #         print("No times within error!")
    #     else:
    #         print("Times:")
    #         print("Minimum Time: {:.4f} s".format(np.nanmin(alltimes)))
    #         print("Maximum Time: {:.4f} s".format(np.nanmax(alltimes)))
    #         print("Time of Minimum Error: {:.4f} s".format(T))


    h_err = np.sqrt((x - D[0])**2 + (y - D[1])**2)
    v_err = np.abs(z - D[2])

    if h_err <= h_tol and v_err <= v_tol:

        # Good Arrival
        R = [T, f_opt[0], f_opt[1], x_opt]



    else:
        R = [np.nan, np.nan, np.nan, x_opt]
        tr = [[np.nan, np.nan, np.nan, np.nan]]

    if trace:
        if particle_output:
            return R, tr, f_particle 
        return R, tr
    else:
        if particle_output:
            return R, f_particle
        return R


if __name__ == '__main__':

    # pass  
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
    D = Position(46, 45.5, 0)
    r = cyscan(S, D, z_profile, trace=True, plot=True)
