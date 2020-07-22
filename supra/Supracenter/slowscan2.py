
# Based off of Wayne Edwards' ray-tracing algorithm (2003)               
# Finds path between two locations with an atmospheric profile in between #
###########################################################################

import warnings

import numpy as np
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})
from supra.Supracenter.cynwDir import nwDir
 


def cyscan(supra_pos, detec_pos, z_profile, wind=True, n_theta=180, n_phi=180, h_tol=330, v_tol=2000):
    # switched positions (Jun 2019)

    # This function should be called for every station
    # Original Author: Wayne Edwards

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
    
    
    # Switch to turn off winds
    if not wind:
        z_profile[:, 2] = 0

    found = False

    ### Initialize variables ###

    # Initial grid spacing
    # Azimuth angle (0 - 360 degrees from North due East)
    dtheta = np.pi/n_theta
    
    # Takeoff angle (90 - 180 degrees from vertical)
    dphi = np.pi/2/n_phi   

    # Horizontal distance between source and a station.
    dx = detec_pos[0] - supra_pos[0]
    dy = detec_pos[1] - supra_pos[1]

    # azth - initial guess for azimuth
    azth = np.arctan2(dy, dx)

    # The number of layers in the integration region
    n_layers = len(z_profile)

    # Slowness, as defined in SUPRACENTER on pg 35, s = 1/c
    s = 1.0/z_profile[0:n_layers, 1]

    # Elevation for that layer
    z  = z_profile[0:n_layers, 0]

    # Set up grid of angles
    phi = np.linspace(azth-np.pi/2, azth+np.pi/2, n_phi)
    Phi = np.tile(phi, (n_theta, 1))

    # Component of wind vector in the direction of phi and phi + pi/2 respectively
    u = nwDir(z_profile[:, 2], z_profile[:, 3], phi)
    v = nwDir(z_profile[:, 2], z_profile[:, 3], phi+np.pi/2)

    # Construct ray parameter net
    # Theta - the take-off angle, of the ray from the source to the station
    theta = np.linspace(np.pi/2, np.pi, n_theta)

    # move theta off of the singularity at pi/2
    theta[0] += 1e-6

    s_val = s[n_layers-1]

    Theta = np.tile(theta, (n_phi, 1)).T

    # Component of wind along the phi direction (azimuth)
    u0 = np.tile(u[n_layers - 1, :], (n_theta, 1))

    # ray parameter
    p = s_val*np.sin(Theta)/(1 + s_val*u0*np.sin(Theta))

    # Transformed x and y
    X = np.zeros((n_theta, n_phi))
    Y = np.zeros((n_theta, n_phi))

    # Transformed wind componenets
    U = np.empty((n_theta, n_phi))
    V = np.empty((n_theta, n_phi))

    #Travel time
    t_arrival = 0

    #azimuth angle
    azimuth = 0

    #takeoff angle
    takeoff = 0

    # ignore negative roots
    last_error = 1e20
    np.seterr(divide='ignore', invalid='ignore')

    ### Scan Loop ###
    while not found:
        trace=[]
        a, b = np.cos(Phi), np.sin(Phi)
        last_z = 0
        for i in range(n_layers - 1):
            
            s2 = s[i]**2
            delz = z[i + 1] - z[i]

            # Winds Enabled
            if wind:
                # clear old variables
                
                # Wind transformation variables
                U = np.tile(u[i, :], (n_theta, 1))
                V = np.tile(v[i, :], (n_theta, 1))
                
                p2 = p/(1 - p*U)
                # This term produces nans
                A = delz/np.sqrt(s2 - p2**2)

                if np.isnan(A).all():

                    break

                # Equation (10)
                X += (p2 + s2*U)*A

                # Equation (11)
                Y += s2*V*A

                # Calculate true destination positions (transform back)
                #0.0016s
                last_z = i + 1

            # Winds Disabled
            else:

                # Equation (3)
                X += p*(delz)/(np.sqrt(s2 - p**2))
                last_z = i + 1
                # Calculate true destination positions (transform back)

            # x = supra_pos[0] + a*X - b*Y
            # y = supra_pos[1] + b*X + a*Y

        E = np.sqrt(((a*X - b*Y - dx)**2 + (b*X + a*Y - dy)**2 + (z[n_layers - last_z - 1] - detec_pos[2])**2)) 

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            k, l = np.where(E == np.nanmin(E))

        # Check for all nan error function
        if k.shape == (0, ):
            # As handled in original Supracenter
            return np.array([np.nan, np.nan, np.nan, np.nan])

        # If there are mulitple, take one closest to phi (in the middle)
        if len(k > 1):
            k, l = k[len(k)//2], l[len(l)//2]
            

        # Compare these destinations with the desired destination, all imaginary values are "turned rays" and are ignored
        

        # Ignore all nan slices - not interested in these points
        # Find the position of the smallest value

        
        new_error = E[k, l]
        # #print(abs(new_error - last_error)/(last_error + 1e-6))
        # if (abs(new_error - last_error)/(last_error + 1e-25) < 1) or (new_error > last_error):
        #     # pass on the azimuth & ray parameter information for use in traveltime calculation
        #     if new_error <= tol or dtheta < precision or dphi < precision:
        #         found = True
        #     else: 
        #         # As handled in original Supracenter
        #         return np.array([np.nan, np.nan, np.nan])

        if E[k, l] < v_tol or dtheta < h_tol or dphi < h_tol:

            # pass on the azimuth & ray parameter information for use in traveltime calculation
            if E[k, l] < v_tol:
                # trace =
                trace=[[supra_pos[0], supra_pos[1], supra_pos[2]]]
                a, b = np.cos(Phi), np.sin(Phi)
                last_z = 0
                for i in range(n_layers - 1):
                    
                    s2 = s[i]**2
                    delz = z[i + 1] - z[i]

                    # Winds Enabled
                    if wind:
                        # clear old variables
                        
                        # Wind transformation variables
                        U = np.tile(u[i, :], (n_theta, 1))
                        V = np.tile(v[i, :], (n_theta, 1))
                        
                        p2 = p/(1 - p*U)
                        # This term produces nans
                        A = delz/np.sqrt(s2 - p2**2)

                        if np.isnan(A).all():

                            break

                        # Equation (10)
                        X += (p2 + s2*U)*A

                        # Equation (11)
                        Y += s2*V*A

                        # Calculate true destination positions (transform back)
                        #0.0016s
                        last_z = i + 1

                    # Winds Disabled
                    else:

                        # Equation (3)
                        X += p*(delz)/(np.sqrt(s2 - p**2))
                        last_z = i + 1
                        # Calculate true destination positions (transform back)

                    x = supra_pos[0] + np.cos(Phi)*X + np.cos(Phi + np.pi/2)*Y
                    y = supra_pos[1] + np.sin(Phi)*X + np.sin(Phi + np.pi/2)*Y
                    trace.append([x[k, l], y[k, l], z[n_layers - 1 - last_z]])

                    #trace.append([x[k, l], y[k, l], z[n_layers - last_z]])
                # print("Azimuth = {:}".format(np.degrees(np.arctan2(-trace[-1][0], -trace[-1][1]))))
                # print("Elevation = {:}".format(np.degrees(np.arctan2(trace[-1][2], np.sqrt((trace[-1][0])**2 + (trace[-1][1])**2)))))
                found = True
            else:
                return np.array([np.nan, np.nan, np.nan, np.nan, np.nan])

        else:
            ### FAST PART ###
            last_error = E[k, l]
            # reduce evenly in both directions
            n_phi = n_theta

            # General Case: central take off angle is between 0 & 90 degrees
            if ((theta[k] != np.pi/2) and (theta[k] != np.pi)):

                # Respace net around best value
                phi = np.linspace(phi[l] - dphi, phi[l] + dphi, n_phi)    
                dphi = 2*dphi/n_phi

                # Check: theta must be > 90 degrees
                if theta[k] - dtheta < np.pi/2:
                    theta = np.linspace(np.pi/2, theta[k] + 2*dtheta, n_theta)
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
            elif (theta[k] == np.pi/2):

                # Respace net around best value
                phi = np.linspace(phi[l] - dphi, phi[l] + dphi, n_phi)         
                dphi = 2*dphi/n_phi  

                theta = np.linspace(np.pi/2, theta[k] + 2*dtheta, n_theta) 
                dtheta = dtheta/n_theta/2
            
            # Update values, and try again
            u = nwDir(z_profile[:, 2], z_profile[:, 3], phi)
            v = nwDir(z_profile[:, 2], z_profile[:, 3], phi + np.pi/2)

            n_theta = len(theta)
            n_phi = len(phi)

            # redefine variables
            Phi = np.tile(phi, (n_theta, 1))
            Theta = np.tile(theta, (n_phi, 1)).T
            u0 = np.tile(u[0, :], (n_theta, 1))


            ### END FAST PART ###
            #####

        # The minimum becomes the center of the next net with borders one spacing away

        X = np.zeros((n_theta, n_phi))
        Y = np.zeros((n_theta, n_phi))

        p = s_val*np.sin(Theta)/(1 + u0*np.sin(Theta)*s_val)



    ######################
    ### 0.0033s

    # Final solution for initial azimuth angle

    # Rotate coordinate system 90 degrees CCW
    # Flip coordinate system horizontally
    azimuth = (450 - phi[l]*180/np.pi)%360

    # Final solution for intial takeoff angle
    takeoff = (theta[k]*180/np.pi)%360

    p1 = p[k, l]
    p2 = p[k, l]**2
    # Find sum of travel times between layers (z)
    for i in range(n_layers - 1 - n_layers + last_z):

        s2 = s[i]**2
        # Equation (9)
        t_arrival += (s2/np.sqrt(s2 - p2/(1 - p1*u[i, l])**2))*(z[i + 1] - z[i])
    
    ##########################
    return t_arrival, azimuth, takeoff, E[k, l], trace

if __name__ == '__main__':
    s = np.array([0, 0, 2])
    d = np.array([0, 1.4, 0])
    z_profile = np.array([ [0, 330, 1, 45],
                        [0.33, 320, 2, 45],
                        [0.67, 315, 2, 45],
                        [1.00, 310, 1, 42],
                        [1.33, 300, 4, 35],
                        [1.67, 293, 5, 30],
                        [2.00, 295, 10, 15]])

    a, b, c, d, e = cyscan(s, d, z_profile)
    print(e)