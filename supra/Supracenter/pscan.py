from supra.Utils.pso import pso
import multiprocessing
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Supracenter.pcyscan import speedscan

from supra.Supracenter.cyzInteg import zInterp 

def mag(vect):
    return np.sqrt(np.dot(vect, vect))

def norm(vect):
    return vect/mag(vect)

# def speedscan(x, *args):
#   # unpack all values
#     import time
    
#     S, D, z_profile = args

#     phi = x[0]
#     theta = x[1]

#     propagating = True

#     current_layer = 0
#     max_layer = len(z_profile)
    
#     p_0 = S
#     a_theta = np.sin(theta)
#     b_theta = np.cos(theta)
    
#     while propagating:
#         a = time.time()
#         if phi > np.pi/2:
#             direction = 1
#         elif phi < np.pi/2:
#             direction = -1

#         a_phi = np.sin(phi)
#         b_phi = np.cos(phi)

#         dz = z_profile[current_layer + direction, 0] - z_profile[current_layer, 0]

#         n = np.array([a_theta*a_phi, b_theta*a_phi, b_phi])*z_profile[current_layer, 1]

#         u = z_profile[current_layer, 2]
#         v = z_profile[current_layer, 3] 

#         V = n + np.array([u, v, 0])

#         p_0 = dz*V/V[2] + p_0

#         c_0 = mag(V)

#         c_s = z_profile[current_layer + direction, 1]
#         u = z_profile[current_layer + direction, 2]
#         v = z_profile[current_layer + direction, 3]

#         c_eff = c_s*n + np.array([u, v, 0])

#         c_1 = mag(c_eff)

#         try:
#             phi = np.pi - np.arcsin(c_1/c_0*a_phi)
#             current_layer += direction
#         except:
#             error = mag(D - p_0)
#             return error


#         error = mag(D - p_0)
        
#         if current_layer + 1 == max_layer or p_0[2] < 0:
#             propagating = False
#         b = time.time()
#     print((b - a)*100000)
#     return error

def pscan(x, *args, **kwargs):

    more_data = kwargs['more_data']

    # unpack all values
    S, D, z_profile = args

    phi = x[0]
    theta = x[1]

    propagating = True

    current_layer = 0
    max_layer = len(z_profile)
    
    p_0 = S
    t = 0
    a = []

    while propagating:
        
        if phi > np.pi/2:
            direction = 1
        elif phi < np.pi/2:
            direction = -1

        dz = z_profile[current_layer + direction, 0] - z_profile[current_layer, 0]

        n = np.array([np.sin(theta)*np.sin(phi), np.cos(theta)*np.sin(phi), np.cos(phi)])

        c_s = z_profile[current_layer, 1]

        u = z_profile[current_layer, 2]
        v = z_profile[current_layer, 3] 

        w = np.array([u, v, 0])

        V = n*c_s + w

        m = norm(V)

        s = dz/m[2]

        p = s*m + p_0

        if more_data:
            t += mag(s*m)/mag(V)
            a.append(p)

        c_0 = mag(V)

        c_s = z_profile[current_layer + direction, 1]
        u = z_profile[current_layer + direction, 2]
        v = z_profile[current_layer + direction, 3]

        w = np.array([u, v, 0])

        c_eff = c_s*n + w

        c_1 = mag(c_eff)

        if c_1/c_0*np.sin(phi) > 1: # wave reflects
            error = mag(D - p)

            if more_data:
                return t, a
            else:
                return error

        else: # wave refracts
            phi_1 = np.pi - np.arcsin(c_1/c_0*np.sin(phi))
            current_layer += direction

        p_0 = p
        phi = phi_1

        error = mag(D - p)
        
        if current_layer + 1 == max_layer or p[2] < 0:
            propagating = False

    if more_data:
        return t, a
    else:
        return error

def psoRayTrace(S, D, z_profile, more_data=True, all_data=False):

    # z_profile [height, speed of sound, u comp, v comp]
    phi_min = 90.1
    phi_max = 179.9
    
    dx = D[0] - S[0]
    dy = D[1] - S[1]
    theta_approx = np.arctan2(dx, dy) # this works in NDE system

    theta_min = theta_approx - 90
    theta_max = theta_approx + 90

    z_profile = np.flipud(z_profile)
    args = (S, D, z_profile)

    x_opt, f_opt, x, f = pso(speedscan, [np.radians(phi_min), np.radians(theta_min)], [np.radians(phi_max), np.radians(theta_max)], args=args, \
        swarmsize=100, maxiter=100, processes=multiprocessing.cpu_count(), particle_output=True, phip=0.5, phig=0.5, omega=0.5, minstep=1e-4, minfunc=1e-4)

    t, a = pscan(x_opt, S, D, z_profile, more_data=True)
    a = np.array(a)
    x_opt = [np.degrees(i)%360 for i in x_opt]

    # fig = plt.figure(figsize=plt.figaspect(0.5))
    # fig.set_size_inches(5, 5)
    # ax = fig.add_subplot(111, projection='3d')
    # ax.plot3D(a[:, 0], a[:, 1], a[:, 2])
    # plt.show()

    if all_data:
        return t, x_opt[1], x_opt[0], a, f_opt
    elif more_data:
        return t, x_opt[1], x_opt[0], f_opt
    else:
        return t, f_opt


if __name__ == '__main__':

    S = np.array([0.0, 0.0, 10000.0])
    D = np.array([10000.0, 0.0, 0.0])

    z_profile = np.array([[    0.0, 330.0, 0.0, 0.0],
                      [ 1000.0, 330.0, 0.0, 0.0],
                      [ 2000.0, 330.0, 0.0, 0.0],
                      [ 3000.0, 330.0, 0.0, 0.0],
                      [ 4000.0, 330.0, 0.0, 0.0],
                      [ 5000.0, 330.0, 0.0, 0.0],
                      [ 6000.0, 330.0, 0.0, 0.0],
                      [ 7000.0, 330.0, 0.0, 0.0],
                      [ 8000.0, 330.0, 0.0, 0.0],
                      [ 9000.0, 330.0, 0.0, 0.0],
                      [10000.0, 330.0, 0.0, 0.0]])

    z_profile = zInterp(D[2], S[2], z_profile, div=100)

    x, a, b, c = psoRayTrace(S, D, z_profile)
    print(90-a, 135-b)
