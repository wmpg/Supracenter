""" Analysis tool for fireballs.

Features:
- loads the trajectory pickle
- plots the dynamic mass calculated from the deceleration as seen from every station (or a given station)

"""

from __future__ import print_function, division, absolute_import

import os
import sys

import scipy.optimize
import numpy as np
import matplotlib.pyplot as plt

from wmpl.Utils.Pickling import loadPickle
from wmpl.Utils.Physics import dynamicMass



def expFunc(x, a, b, c, d):
    """ Jacchia deceleration function. """
    
    return a + b*x - abs(c)*np.exp(d*x)



def expTwoFunc(x, a, b, c, d, e, f):
    """ Two term Jacchia deceleration function. """
    
    return a + b*x - abs(c)*np.exp(d*x) - abs(e)*np.exp(f*x)



def expFuncDer(x, a, b, c, d):
    """ Derivative of Jacchia deceleration function. """

    return b - abs(c)*abs(d)*np.exp(d*x)



def expTwoFuncDer(x, a, b, c, d, e, f):
    """ Derivative of the two term Jacchia deceleration function. """
    
    return b - abs(c*d)*np.exp(d*x) - abs(e*f)*np.exp(f*x)



def expFunc2ndDer(x, a, b, c, d):
    """ 2nd derivative of Jacchia deceleration function. """

    return -abs(c*d*d)*np.exp(d*x)



def expTwoFunc2ndDer(x, a, b, c, d, e, f):
    """ 2nd derivative of the two term Jacchia deceleration function. """

    return -abs(c*d*d)*np.exp(d*x) - abs(e*f*f)*np.exp(f*x)




def polyFunc(x, a, b, c, d):
    """ 3rd order polynomial. """

    return a + x*b + c*x**2 + d*x**3



def polyFuncDer(x, a, b, c, d):

    return b + 2*c*x + 3*d*x**2



def polyFunc2ndDer(x, a, b, c, d):

    return 2*c + 6*d*x



def residuals(params, x, y, fitFunc=expTwoFunc):
    """ Returns the residuals between the predicted and input values of the 2 term exonential.

    Arguments:
        params: [ndarray] function parameters
        x: [ndarray] independant variable
        y: [ndarray] prediction

    Return:
        residuals: [ndarray]
    """

    return fitFunc(x, *params) - y



def residuals_minimize(params, x, y, fitFunc=expTwoFunc):
    """ Wrapper function for calculating fit residuals for minimization. """

    # Squared value of each residual
    z = residuals(params, x, y, fitFunc=fitFunc)**2

    # Smooth approximation of l1 (absolute value) loss
    return np.sum(2*((1 + z)**0.5 - 1))



if __name__ == "__main__":

    ### INPUT FILE ######

    # Trajectory path
    #dir_path = os.path.abspath("../MILIG files/20170923_053525 meteorite dropping/Monte Carlo")
    #dir_path = os.path.abspath("../MILIG files/20171127_meteorite_dropping/Monte Carlo")
    dir_path = os.path.abspath("../MILIG files/20180125_meteorite_dropping/Monte Carlo")

    # Trajectory pickle file
    # traj_file = '20170923_053524_mc_trajectory.pickle'
    # traj_file = '20171127_010546_mc_trajectory.pickle'
    traj_file = '20180125_002305_mc_trajectory.pickle'

    # No. of input station (if -1, the program will loop over all stations)
    station_id = -1


    ### INPUT METEOROID PARAMETERS ###
    
    # Bulk density of meteoroid
    bulk_density = 3500 # kg/m^3

    # Shape factor
    # - sphere      = 1.21
    # - hemisphere  = 1.92
    # - cube        = 1.0
    # - brick 2:3:5 = 1.55
    shape_factor = 1.21

    # Drag coefficient
    gamma = 1.0


    # Fit function
    # - 'exp' - single term exponential
    # - '2exp' - two term exponential
    # - 'poly' - 3rd order polynomial
    fit_func = 'poly'


    ##########################################################################################################



    # Choose the fitting function, velocity, deceleration
    if fit_func == 'exp':

        # Fit function
        fitFunc = expFunc
        fitFuncDer = expFuncDer
        fitFunc2ndDer = expFunc2ndDer

        # Number of function arguments
        fit_func_args = 4

    elif fit_func == '2exp':

        # Fit function
        fitFunc = expTwoFunc
        fitFuncDer = expTwoFuncDer
        fitFunc2ndDer = expTwoFunc2ndDer

        # Number of function arguments
        fit_func_args = 6

    else:

        # Fit function
        fitFunc = polyFunc
        fitFuncDer = polyFuncDer
        fitFunc2ndDer = polyFunc2ndDer

        # Number of function arguments
        fit_func_args = 4


    # Load the trajectory file
    traj = loadPickle(dir_path, traj_file)


    # Check that the given station ID exists
    if (station_id < -1) and (station_id > len(traj.observations) - 1):
        print('station_id should be an integer in the range [0 to {:d}]!'.format(len(traj.observations) - 1))
        sys.exit()


    # List of all time data, fits and dynamic masses, for comparison at the end
    time_data_list = []
    height_data_list = []
    fit_list = []
    length_data_list = []
    dynamic_mass_list = []

    for stat_i, obs in enumerate(traj.observations):

        # If the ID of the station was not given, loop over all observations
        if station_id == -1:
            pass

        elif station_id != stat_i:
            continue


        # Extract time vs. length from the best station
        print(traj.observations[stat_i].station_id)
        time_data = traj.observations[stat_i].time_data
        length_data = traj.observations[stat_i].length
        height_data = traj.observations[stat_i].meas_ht


        ### DATA SELECTION ###

        # # Take only the last quarter of observations
        # last_quart_len = int(0.75*len(time_data))

        # time_data = time_data[last_quart_len:]
        # length_data = length_data[last_quart_len:]
        # height_data = height_data[last_quart_len:]


        ######################




        x0 = np.ones(fit_func_args)


        # Treat the fit as a minimization problem, but use basinhopping for minimizing residuals
        fit_robust_mini = scipy.optimize.basinhopping(residuals_minimize, x0, \
            minimizer_kwargs={'args':(time_data, length_data, fitFunc)}, niter=500)

        # Print the fit status message
        print('Fit result:')
        print(fit_robust_mini)


        # Extract the fitted parameters
        popt = fit_robust_mini.x


        # Plot original data
        plt.plot(time_data, length_data/1000, marker='x', label='Observations')

        # Plot the fit
        t_plot = np.linspace(np.min(time_data), np.max(time_data), 1000)
        plt.plot(t_plot, fitFunc(t_plot, *popt)/1000, label='Model')

        plt.xlabel('Time (s)')
        plt.ylabel('Length (km)')

        plt.legend()

        plt.title('Station ' + str(obs.station_id))

        plt.show()

        plt.clf()
        plt.close()


        mass_dyn_array = []

        # Calcualte the dynamic mass at every point in time
        for t, lat, lon, ele, jd in zip(time_data, traj.observations[stat_i].meas_lat, \
            traj.observations[stat_i].meas_lon, traj.observations[stat_i].meas_ht, \
            traj.observations[stat_i].JD_data):

            # Calculate the velocity and deceleration from the fit
            vel = fitFuncDer(t, *popt)
            decel = np.abs(fitFunc2ndDer(t, *popt))

            # Calculate the dynamic mass
            mass_dyn = dynamicMass(bulk_density, lat, lon, ele, jd, vel, decel, gamma=gamma, \
                shape_factor=shape_factor)

            mass_dyn_array.append(mass_dyn)

        mass_dyn_array = np.array(mass_dyn_array)



        # Add data to cummulative lists
        time_data_list.append(time_data)
        length_data_list.append(length_data)
        height_data_list.append(height_data)
        fit_list.append(popt)
        dynamic_mass_list.append(mass_dyn_array)


        ### PLOT THE DYNAMIC MASS ###
        ######################################################################################################

        fig = plt.figure()

        ax1 = fig.add_subplot(1, 4, 1)
        ax2 = fig.add_subplot(1, 4, 2, sharey=ax1)
        ax3 = fig.add_subplot(1, 4, 3, sharey=ax1)
        ax4 = fig.add_subplot(1, 4, 4, sharey=ax1)
            
        # Plot the fit residuals
        ax1.plot(length_data - fitFunc(time_data, *popt), height_data/1000, zorder=3)

        ax1.set_xlabel('Residuals (m)')
        ax1.set_ylabel('Height (km)')


        # Plot the dynamic mass
        ax2.plot(mass_dyn_array*1000, height_data/1000, zorder=3)

        ax2.set_xlabel('Dynamic mass (g)')


        # Plot Velocity
        ax3.plot(fitFuncDer(time_data, *popt), height_data/1000, zorder=3)

        ax3.set_xlabel('Velocity (m/s)')

    
        # Plot decelration
        ax4.plot(fitFunc2ndDer(time_data, *popt), height_data/1000, zorder=3)

        ax4.set_xlabel('Deceleration (m/s^2)')



        plt.subplots_adjust(wspace=0)

        # Turn of Y ticks on the other plots
        ax2.tick_params(
            axis='y',          # changes apply to the y-axis
            which='both',      # both major and minor ticks are affected
            left='off',        # ticks along the left edge are off
            labelleft='off')  # labels along the left edge are off

        ax3.tick_params(
            axis='y',          # changes apply to the y-axis
            which='both',      # both major and minor ticks are affected
            left='off',        # ticks along the left edge are off
            labelleft='off')  # labels along the left edge are off

        ax4.tick_params(
            axis='y',          # changes apply to the y-axis
            which='both',      # both major and minor ticks are affected
            left='off',        # ticks along the left edge are off
            labelleft='off')  # labels along the left edge are off



        ax1.grid(color='0.9')
        ax2.grid(color='0.9')
        ax3.grid(color='0.9')
        ax4.grid(color='0.9')

        plt.suptitle('Station ' + str(obs.station_id))


        plt.show()

        plt.clf()
        plt.close()

        ######################################################################################################



    # If all dynamic masses were calculated, plot a graph showing all
    fig = plt.figure()

    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2, sharey=ax1)

    for obs, time_data, popt, length_data, height_data, mass_dyn_array in zip(traj.observations, \
        time_data_list, fit_list, length_data_list, height_data_list, dynamic_mass_list):
            
        # Plot the fit residuals
        ax1.plot(length_data - fitFunc(time_data, *popt), height_data/1000, zorder=3)

        # Plot the dynamic mass
        ax2.plot(mass_dyn_array*1000, height_data/1000, zorder=3, label=str(obs.station_id))

        
    plt.subplots_adjust(wspace=0)

    # Turn of Y ticks on the second plot
    ax2.tick_params(
        axis='y',          # changes apply to the y-axis
        which='both',      # both major and minor ticks are affected
        left='off',        # ticks along the left edge are off
        labelleft='off')  # labels along the left edge are off



    ax1.set_xlabel('Residuals (m)')
    ax1.set_ylabel('Height (km)')

    ax2.set_xlabel('Dynamic mass (g)')

    

    ax1.grid(color='0.9')
    ax2.grid(color='0.9')

    plt.suptitle('All stations')

    plt.legend()


    plt.show()
if __name__ == '__main__':

    pass