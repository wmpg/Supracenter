
from supra.Geminus.overpressure2 import overpressureihmod_Ro
from supra.GUI.Tools.GUITools import *
from supra.Utils.pso import pso


SWARM_SIZE = 5
MAX_ITER = 25

def overpressureErr(Ro, *args):

    """ Function to optimize for period or pressure
    """

    source_list, stat, v, theta, dphi, sounding_pres, sw, wind, dopplershift, target, mode, regime = args


    tau, tauws, Z, sR, inc, talt, dpws, dp, it = overpressureihmod_Ro(source_list, stat, Ro[0], v, theta, dphi, sounding_pres, sw, wind=wind, dopplershift=dopplershift)

    if mode == 'pres':
        if regime == 'ws':
            err = abs(target - dpws[-1])
        else:
            err = abs(target - dp[-1])
    else:
        if regime == 'ws':
            err = abs(target - tauws[-1])
        else:
            err = abs(target - tau[-1])


    return err

def periodSearch(p, gem_inputs, paths=False):

    ''' Uses PSO to find the Relaxation radius that returns the desired period through the Geminus program
    '''

    Ro = 10.0

    target_period = p

    
    period_ws = 0

    tol = 1e-3

    
    tau = []

    search_min = [tol]
    search_max = [30]

    Ro, f_opt = pso(overpressureErr, search_min, search_max, \
        args=gem_inputs + [p, 'period', 'ws'], swarmsize=SWARM_SIZE, maxiter=MAX_ITER, processes=1, minfunc=tol, minstep=1e-3) 

    Ro = Ro[0]
    print("Period Weak Shock: {:.2f} s Ro = {:.3f} m".format(p, Ro))

    Ro_ws = Ro

    Ro = 10.0
    period_lin = 0

    Ro, f_opt = pso(overpressureErr, search_min, search_max, \
        args=gem_inputs + [p, 'period', 'lin'], swarmsize=SWARM_SIZE, maxiter=MAX_ITER, processes=1, minfunc=tol, minstep=1e-3) 
    Ro = Ro[0]
    print("Period Linear: {:.2f} s Ro = {:.3f} m".format(p, Ro))
    Ro_lin = Ro


    if paths:
        source_list, stat, v, theta, dphi, sounding_pres, sw, wind, dopplershift = gem_inputs
        tau, tauws, Z, sR, inc, talt, dpws, dp, it = overpressureihmod_Ro(source_list, stat, Ro_ws, v, theta, dphi, sounding_pres, sw, wind=wind, dopplershift=dopplershift)
        weak_path = tau[:it] + tauws[it:]
        tau, tauws, Z, sR, inc, talt, dpws, dp, it = overpressureihmod_Ro(source_list, stat, Ro_lin, v, theta, dphi, sounding_pres, sw, wind=wind, dopplershift=dopplershift)
        lin_path = tau

        return Ro_ws, Ro_lin, weak_path, lin_path, tau, Z, it

    return Ro_ws, Ro_lin

def presSearch(p, gem_inputs, paths=False):
    Ro = 10.0

    target_pres = p
    
    pres_ws = 0

    tol = 1e-3

    source_list, stat, v, theta, dphi, sounding_pres, sw, wind, dopplershift = gem_inputs
    search_min = [tol]
    search_max = [30]
    Ro, f_opt = pso(overpressureErr, search_min, search_max, \
        args=gem_inputs + [p, 'pres', 'ws'], swarmsize=SWARM_SIZE, maxiter=MAX_ITER, processes=1, minfunc=tol, minstep=1e-3)
    Ro = Ro[0]
    print("Pressure Weak Shock: {:.2f} mPa Ro = {:.3f} m".format(p*1000, Ro))

    Ro_ws = Ro

    Ro = 10.0
    pres_lin = 0
    Ro, f_opt = pso(overpressureErr, search_min, search_max, \
        args=gem_inputs + [p, 'pres', 'lin'], swarmsize=SWARM_SIZE, maxiter=MAX_ITER, processes=1, minfunc=tol, minstep=1e-3)
    Ro = Ro[0]  
    print("Pressure Linear: {:.2f} mPa Ro = {:.3f} m".format(p*1000, Ro))


    Ro_lin = Ro



    if paths:
        source_list, stat, v, theta, dphi, sounding_pres, sw, wind, dopplershift = gem_inputs
        tau, tauws, Z, sR, inc, talt, dpws, dp, it = overpressureihmod_Ro(source_list, stat, Ro_ws, v, theta, dphi, sounding_pres, sw, wind=wind, dopplershift=dopplershift)
        weak_path = tau[:it] + tauws[it:]
        tau, tauws, Z, sR, inc, talt, dpws, dp, it = overpressureihmod_Ro(source_list, stat, Ro_lin, v, theta, dphi, sounding_pres, sw, wind=wind, dopplershift=dopplershift)
        lin_path = tau
        return Ro_ws, Ro_lin, weak_path, lin_path, tau, Z, it

    return Ro_ws, Ro_lin