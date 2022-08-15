
import numpy as np
import matplotlib.pyplot as plt

from supra.Utils.pso import pso

from supra.Atmosphere.Pressure import estPressure


class Yield:

    """ A simple yield object used for doing conversions 
    """

    def __init__(self, W):


        W_0_nuc = 4.184e12
        W_0_chem = 4.184e6 

        # yield in Joules
        self.j = W

        # Nuclear Yield in kT NE
        self.nuclear = W/W_0_nuc

        # chemcial Yield in kg HE
        self.chem = W/W_0_chem

        # yield scaling factor
        self.ysf = (self.nuclear/1)**(1/3)

    def __str__(self):

        return "Yield Object: {:.3e} J, {:.3f} kg HE, {:.3f} kT NE".format(self.j, self.chem, self.nuclear)


def r02HOB(R_0):
    """ Scaled burst height HE to Height of Burst (HOB) Yield Factor

        See Figure 10 in ANSI
    """

    hobFunc = np.poly1d(np.load("supra\\Yields\\HOB_curve.npy"))

    z = hobFunc(R_0)

    # if the fireball is high enough not to have interactions with the ground
    # Called "Free Air Burst"
    if z < 1:
        z = 1

    return z

def transmissionFactor(z, mode="mean"):
    """ Returns transmission factor given in KG85, Table XIV

    z altitude in meters
    mode = mean or alt
    """

    if mode == "mean":
        func = np.poly1d(np.load("supra\\Yields\\f_d_mean_fit.npy"))
    elif mode == "alt":
        func = np.poly1d(np.load("supra\\Yields\\f_d_alt_fit.npy"))
    else:
        print("Unrecognized Mode: {:}".format(mode))
        return None

    return func(z)

def chemOverpressureRatioInv(p_req):

    Z = (1/(p_req) * 0.82738495) + (1/(np.sqrt(p_req + (1/(p_req) * 0.30315))) * 2.496861)
    
    return Z

def chemOverpressureRatio(z):

    if z <= 500:

        r = 1/z

        func = np.poly1d(np.load("supra\\Yields\\overpressure_ratio_fit.npy"))

        return func(r)

    else:

        r = np.log10(z)

        func = np.poly1d(np.load("supra\\Yields\\overpressure_ratio_fit_extended.npy"))

        return 10**func(r)

def nucScaledDistanceFromOverpressure(p_rat):

    r = np.log(p_rat)

    func = np.poly1d(np.load("supra\\Yields\\nuc_scaled_curve_fit.npy"))

    return func(r)

def nucOverpressureFromScaledDistance(Z):

    r = 1/Z

    func = np.poly1d(np.load("supra\\Yields\\nuc_scaled_curve_fit_inv.npy"))

    return func(r)

def scaledHOB(h, W_W_0):
    """ Scaled Height of Burst (HOB) calculation from ANSI B.1.3

    If a large explosion happens at height, then it scales to a standard explosion at
    the scaled height
    """

    scaled_HOB = h/(W_W_0)**(1/3)

    return scaled_HOB

def equivAirburstYield(HOB_yf, W_W_0):

    eay = HOB_yf*W_W_0

    return HOB_yf


def yieldScalingFactor(W_W_0):

    ysf = (W_W_0)**(1/3)

    return ysf


def scaledDistance(R, W_W_0, P_P_0, f_T=1):

    #R_0 = f_T*R*(W_W_0/P_P_0)**(-1/3)

    R_0 = f_T*R*(W_W_0/P_P_0)**(1/3)

    return R_0
    
def getPP0(h):

    # pressure at ground
    P_0 = estPressure(0)

    # pressure at HOB
    P = estPressure(h)

    return P/P_0

def yieldFromScaledDistance(Z, R, P_P_0, f_T=1, mode="kgHE"):

    if mode == "kgHE":
        W_0 = 4.184e6
    elif mode == "ktNE":
        W_0 = 4.184e12
    else:
        print("Unrecognized mode: {:}".format(mode))
        return None

    W = (f_T*R/Z)**(3)*W_0*P_P_0

    return W



def scaledOverpressure(del_p, P_P_0):
    """ Scales overpressure measured at a station with atmospheric pressure differences
    """

    scaled_del_p = del_p*P_P_0

    return scaled_del_p


def kgHE2kTNE(kgHE):

    # convert from kg to kT
    kTHE = 1e-6*kgHE

    # convert from HE to NE
    kTNE = 1.1*kTHE

    return kTNE

def kTNE2kgHE(kTNE):

    # convert from NE to HE
    kTHE = kTNE/1.1

    # convert from kT to kg
    kgHE = 1e6*kTHE

    return kgHE

def overpressureDistance(R, W, P_P_0, mode="kgHE"):

    if mode ==  "kgHE":
        # W in kg HE, R in meters
        del_p = 105.93*W**(11/30)*R**(-11/10)*(P_P_0)**(19/30) # kPa
        del_p *= 1000

    elif mode == "kTNE":

        del_p = 1.349e4*W**(11/30)*R**(-11/10)*(P_P_0)**(19/30) # kPa
        del_p *= 1000

    else:
        print("Unknown mode: {:}".format(mode))
        return None

    return del_p


def distance2Overpressure(R, del_p, P_P_0, mode="kgHE"):

    del_p /= 1000 # convert to kPa

    if mode ==  "kgHE":
        # W in kg HE, R in meters
        W = (del_p/105.93*R**(11/10)*(P_P_0)**(-19/30))**(30/11)
        

    elif mode == "kTNE":
        
        W = (del_p/1.349e4*R**(11/10)*(P_P_0)**(-19/30))**(30/11)

    else:
        print("Unknown mode: {:}".format(mode))
        return None

    return W

