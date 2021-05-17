import numpy as np
import matplotlib.pyplot as plt
import scipy

from supra.Utils.Classes import Constants, Position
from supra.Utils.Formatting import *


c = Constants()

# function [dp,dpws,dpratio,tau,tauws,Z,td,talt,Ro] = overpressureihmod_Ro(meteor,stn,Ro,v,theta,dphi,atmos,sw);
def overpressureihmod_Ro(meteor, stn, Ro, v, theta, dphi, atmos, sw, wind=True, dopplershift=True):


    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # %                                                                                                      
    # %  Theoretical Acoustic Overpressure Prediction using Line Source Theory                               
    # %  (ReVelle D.O., 1976: On Meteor-Generated Infrasound, JGR, 81, pp.1217-1229.                         
    # %                                                                                                      
    # %  Usage:                                                                                              
    # %   [dp,dpws,dpratio,tau,tauws,Z,td,talt,Ro,dm] = overpressureihmod(meteor,stn,mass,rhom,v,theta,dphi,atmos,sw);  
    # %                                                                                                      
    # %  Given: meteor - [latitude,longitude,altitude] of meteor source [DD.ddd,DD.ddd,km]                   
    # %         stn    - [latitude,longitude,elevation] of observing station [DD.ddd,DD.ddd,km]              
    # %         mass   - meteoroid mass in kilograms (kg)                                                    
    # %         rhom   - meteoroid bulk density in kilograms per metre cubed (kg/m^3)             
    # %         v      - meteoroid velocity in kilometres per second (km/s)                       
    # %         theta  - entry angle of meteoroid measured from horizontal (degrees)              
    # %         dphi   - angular deviation from the meteoroid trajectory (degrees)                
    # %         atmos  - atmospheric model (nx3) - MUST encompass altitudes for direct propagation
    # %                    - altitude (m)
    # %                    - pressure (hPa/mbars)
    # %                    - temperature (oC)
    # %         sw     - switches on/off = 1/0       3 element vector [a,b,c]                       
    # %                  a - vary period to find transition (1) const. period (0)                 
    # %                  b - display figures (on/off)                                             
    # %                  c - quick integration (on/off)
    # %                                                                                           
    # %  Returns: dp      - theoretical acoustic/infrasonic overpressure (Pascals)                
    # %           dpws    - theoretical acoustic/infrasonic overpressure for completely           
    # %                     weak shock propagation
    # %           dpratio - acoustic/infrasonic overpressure ratio (dp/p) (dimensionless)         
    # %           td      - transition distance (in units of Ro)                                  
    # %           talt    - transition altitude (kilometres)                                      
    # %           tau     - signal period (seconds)
    # %           tauws   - signal period for completely weak shock propagation (seconds)
    # %           Z       - Altitude interval (km)                                                
    # %           Ro      - Blast radius in metres                                                
    # %           dm      - diameter of meteoroid (metres)                                        
    # %                                                                                           
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    try:
        meteor = validate(meteor, "Meteor Source Position") 
        stn    = validate(stn, "Station Reciever Position")
        Ro     = validate(Ro, "Blast Radius")
        v      = validate(v, "Meteor Trajectory Velocity")
        theta  = validate(theta, "Entry Angle of Meteoroid")
        dphi   = validate(dphi, "Angular Deviation from the Meteoroid Trajectory")
        atmos  = validate(atmos, "Atmospheric Data")
        sw     = validate(sw, "Switch Data")
    except TypeError as e:
        print(printMessage("Error"), " Geminus variable input error! {:}".format(e))
        raise TypeError(e)
        


    alt =  atmos[:, 0]*0.001
    if wind:
        Wm = atmos[:, 2] #% meridional wind velocity 
        Wz = atmos[:, 3] #% zonal wind velocity
    else:
        Wm = np.zeros(len(alt))
        Wz = np.zeros(len(alt))
    pres = atmos[:, 4]*100
    temp = atmos[:, 1]+273.15

    Cs0 = np.sqrt(c.GAMMA*c.R*temp/c.M_0)


    #% SWITCH SELECTIONS for wind/no wind and doppler/no doppler shift
    #disp('Include winds = 1, no winds = 2, include Doppler = 3, no Doppler = 4');
    #disp('Absorption as per ReVelle (1976) = 5');
    #disp('Absorption as per Sutherland and Bass (2004) = 6');
    #%in = input('Please make switch selections [1/2, 3/4, 5/6]:   ');
    Kw = 1
    Kd = 3
    Ka = 5

    wCs, Ceff, s, wtype = windC(meteor, stn, atmos, Cs0)
    rho = c.GAMMA*pres/wCs**2


    Csa = Cs0[0]
    Cswind = wCs[0]
    Ceffs = Ceff[0] 
    Wms = Wm[0]
    Wzs = Wz[0]
    Ps  =   pres[0]
    rhos =  rho[0]
    temps = temp[0]

    M = v/(Csa/1000)

    Ro = Ro/1000

    tau0 = 2.81*Ro/(Csa/1000)
    f0 = 1/tau0             

    theta = np.radians(theta)
    dphi = np.radians(dphi)
    epsilon = np.arctan(np.tan(theta)/(1 - 2*dphi/np.pi))

    Range = (Position(meteor[0], meteor[1], meteor[2]) - Position(stn[0], stn[1], stn[2])).mag()/1000

    dH = meteor[2] - stn[2]
    sR = np.sqrt(Range**2 + dH**2)

    inc = np.arctan(dH/Range)


    Zs = meteor[2]
    xt = sR/Ro
    Z10 = Zs - 10*Ro*np.sin(inc)
    N = 500
    dZ = (Z10 - stn[2])/N


    Z = np.linspace(Z10, stn[2], N)

    f = scipy.interpolate.interp1d(alt, pres)
    Pz = f(Z)

    f = scipy.interpolate.interp1d(alt, Cs0)
    Csz = f(Z)

    f = scipy.interpolate.interp1d(alt, rho)
    rhoz = f(Z)

    f = scipy.interpolate.interp1d(alt, temp)
    tempz = f(Z)

    f = scipy.interpolate.interp1d(alt, wCs)
    CsW = f(Z)

    CmW = meanC(meteor[2], Cswind, Z, CsW)

    # # % Find mean sound speed between source and Z
    # Cmz = meanC(meteor[2], Csa, Z, Csz)

    x = (Zs - Z)/(Ro*np.sin(inc))

    Csd = Cs0;
    Dtau = doppler(meteor, stn, tau0, Csd, atmos, Z, s, alt)
    Dtau = np.flipud(Dtau)

    tau = []
    for jj in range(len(Dtau)):
        if dopplershift:
            tau.append(0.562*Dtau[jj]*x[jj]**(0.25))
        else:
            tau.append(0.562*tau0*x[jj]**(0.25))

    # tau = 0.562*tau0*x**(1/4)
    if sw[0] == 0:
        tau = tau0
    fm = 1/np.array(tau)
    Dtogo = xt - x

    # % Calculate Linear Absorption
    # DL = lineardamping(Z, Csz, rhoz, tempz, fm, epsilon, sw[2])
    # # % Calculate Overpressure ratio via eqn 23 (eqn 93b thesis)
    # dpp = c.GAMMA/(2*(c.GAMMA + 1))*(3/8)**(-3/5)/((1 + (8/3)**(8/5)*x**2 )**(3/8) - 1)

    # dpp = dpp*(rhoz/rhos)**(1/2)*Csz/Cmz  #% Non uniform path correction (eqn 86)
    # dpp = dpp*(rhos/rhoz)*Csa**2/Csz**2;    #% Source Altitude correction for inhomogenious atmosphere

    # dpp = dpp*DL*x**(1/4)

    DL = lineardamping(Z, CsW, rhoz, tempz, fm, epsilon, sw[2])

    dpp = c.GAMMA/(2*(c.GAMMA + 1))*(3/8)**(-3/5)/((1 + (8/3)**(8/5)*x**2)**(3/8) - 1)
    dpp = dpp*(rhoz/rhos)**(1/2)*CsW/CmW
    #% Source Altitude correction for inhomogenious atmosphere
    dpp = dpp*(rhos/rhoz)*Cswind**2/CsW**2
    dpp = dpp*DL*x**(1/4)# % Linear damping & decay
                #% Linear damping & decay
    # % Calculate Distortion distance
    Dprime = Csz*tau/(34.3*dpp) 
 
    Dprime = Dprime/1000
    Dprime = Dprime/Ro


    Trans = (Dprime > Dtogo).astype(int)


    try:

        M, it = 1, np.nanargmax(Trans)
        #Transition here
    except ValueError:
        # Tested from Matlab code
        M, it = 0, 0

    if (M == 1):
        td = x[it]
        talt = Z[it]
    else:            
        td = xt
        talt = stn[2]


    Csd = Cs0

    Dtau = doppler(meteor, stn, tau0, Csd, atmos, Z, s, alt)
    Dtau = np.flipud(Dtau)

    tau = []
    for jj in range(len(Dtau)):
        if dopplershift:
            tau.append(0.562*Dtau[jj]*x[jj]**(0.25))
        else:
            tau.append(0.562*tau0*x[jj]**(0.25))


    if (sw[0] == 0):
        tau = [tau0]*len(x)
    tauws = tau.copy()
    tau[it-1:] = [tau[it-1]]*len(tau[it-1:])

    fm = 1/np.array(tau)
    fmws = 1/np.array(tauws)

    dpp = c.GAMMA/(2*(c.GAMMA + 1))*(3/8)**(-3/5)/((1 + (8/3)**(8/5)*x**2)**(3/8) - 1)

    dpp = dpp*(rhoz/rhos)**(1/2)*CsW/CmW
    dpp = dpp*(rhos/rhoz)*Cswind**2/CsW**2 
    dppws = dpp
    
    Dws = weakshockdamping(Z, CsW, rhoz, tempz, fm, epsilon, Ps, sw[2])
    dpp = dpp*Dws 
    dpp[it:] = dpp[it-1]
    Dws = weakshockdamping(Z, CsW, rhoz, tempz, fmws, epsilon, Ps, sw[2])
    dppws = dppws*Dws

    DL = lineardamping(Z[it:], CsW[it:], rhoz[it:], tempz[it:], fm[it:], epsilon, sw[2])
    pzt = (rhos/rhoz[it-1])*Cswind**2/CsW[it-1]**2
    pzg = (rhos/rhoz[it:])*Cswind**2/CsW[it:]**2
    dppl = dpp[it:]*(pzg/pzt) 
    dppl = dpp[it:]*DL*(td/x[it:])**(1/2)    


    dpratio = dpp

    for ii in range(len(dppl)):
        dpratio[it+ii-1] = dppl[ii]


    dp = dpratio*Pz
    dpws = dppws*Pz

    Ro = Ro*1000

    # if True:
    if sw[1] == 1:
        
        plt.plot(tau[0:it], Z[0:it], 'r-', label="Weak Shock Period Change")
        plt.plot(tau[it-1:], Z[it-1:], 'b-', label="Stable Period")
        plt.plot(tauws[it-1:], Z[it-1:], 'm-', label="Weak Shock: No Transition")

        plt.scatter([tau[it-1]], [Z[it-1]])
        
        plt.xlabel("Signal Period [s]")
        plt.ylabel("Geopotential Height [km]")
        plt.legend()
        plt.show()
        # txtfile = fopen('Weak Shock Output - TEST only.txt','a');
        print('FINAL OUTPUT FOR THE WEAK SHOCK')
        print('=========================================================')
        print('Period (weak shock):     {:3.4f} s'.format(tauws[-1]))
        print('  Frequency (weak shock):   {:3.3f} Hz'.format(1/tauws[-1]))
        print('Period (linear):         {:3.4f} s'.format(tau[-1]))
        print('  Frequency (linear):       {:3.3f} Hz'.format(1/tau[-1]))
        print('Slant range:             {:5.2f} km'.format(sR))
        print('Arrival (inclination):   {:3.4f} deg'.format(np.degrees(inc)))
        print('Transition height:       {:3.3f} km'.format(talt))
        print('Overpressure (weak shock):     {:3.4f} Pa'.format(dpws[-1]))
        print('Overpressure (linear):         {:3.4f} Pa'.format(dp[-1]))


    return [tau, tauws, Z, sR, np.degrees(inc), talt, dpws, dp, it]
    # # Temporary float return
    # return [tau[-1], tauws[-1], dpws[-1], dp[-1], it]

    #======================================================================================================

def meanC(Zs,Cs,Z,Cz):

    Cm = []
    for i in range(len(Z)):

        del_z = (Zs - Z[i-1])
        avgC = -intfunc([Zs, Z[i]], [Cs, Cz[i]])

        Cm.append(avgC/del_z)

    Cm = np.array(Cm)
    return Cm

# %======================================================================================================

def shearviscosity(T):

    mu0 = 1.846e-5
    T0 = 300
    Ts = 110.4
    Ta = 245.4
    Tb = 27.6

    mu = mu0*(T0 + Ts)/(np.array(T) + Ts)*(np.array(T)/T0)**(3/2)

    return mu

# %======================================================================================================

def thermcond(T):

    kap0 = 2.624e-2
    T0 = 300
    Ta = 245.4
    Tb = 27.6

    T = np.array(T)
    T_b_T = []

    for i in T:
        T_b_T.append(np.exp(-Tb/i))
    T_b_T = np.array(T_b_T)


    kap = kap0*(T/T0)**(3/2)*(T0 + Ta*T_b_T)/(T + Ta*T_b_T)
    return kap

# %======================================================================================================

def bulkviscosity(T):

    nu = 2/3*shearviscosity(T);
    return nu

# %======================================================================================================

def lineardamping(Z, Cs, rho, T, fm, epsilon, Method):

    gamma = 1.4
    Cp = 1008.56
    mu = shearviscosity(T)
    kappa = thermcond(T)
    nu = bulkviscosity(T)

    Cs = np.array(Cs)

    delta = 4*(4/3*mu + nu + kappa*(gamma-1)/Cp)
    lamda = Cs/np.array(fm)

    alpha = np.pi**2*delta/(2*rho*Cs*lamda**2)

    integ = []
    for i in range(len(Z)):
        integ.append(-intfunc(Z[:i], alpha[:i]/np.cos(epsilon)))

    DL = np.exp(-np.array(integ))
    return DL

# %======================================================================================================

def weakshockdamping(Z, Cs, rho, T, fm, epsilon, Ps, Method):

    gamma = 1.4
    Cp = 1008.56
    mu = shearviscosity(T)
    kappa = thermcond(T)
    nu = bulkviscosity(T)
    delta = 4*( 4/3*mu + nu + kappa*(gamma-1)/Cp)
    BA = 3/2*delta*fm/(gamma+1)
    l = Cs/fm
    Be = 3*delta/( 2*rho*Cs*l**2)
    integ = []
    for i in range(len(Z)):
        integ.append(-intfunc(Z[0:i], Be[0:i]/np.cos(epsilon)))
    
    dpz = 0.0575*Ps
    DWS = BA*np.exp(-np.array(integ))/(dpz*(1 - np.exp(-np.array(integ))) + BA)
    return DWS


# %======================================================================================
def intfunc(X, Y):
    # tested and works with silber version

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %  Usage: A = intfunc(X,Y)
# %
# %  Integrate an arbitrary function Y(X) defined by
# %  measured points X & Y assuming piecewise linearity
# %
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    dm = np.array(Y[1:]) - np.array(Y[0:-1])
    dx = np.array(X[1:]) - np.array(X[0:-1])
    dx2 = np.array(X[1:])**2 - np.array(X[0:-1])**2
    slope = []
    for i in range(len(dm)):
        slope.append(dm[i]/dx[i])
    slope = np.array(slope)
    inpt = Y[0:-1] - slope*X[0:-1]

    segment = slope*dx2/2 + inpt*dx

    A = np.sum(segment)
    

    return A



    

def windC(source, stn, atmos, Cs):

    sP = Position(*source)
    stP = Position(*stn)

    waveP = stP - sP
    x, y, z = waveP.x, waveP.y, waveP.z

    h = np.array(atmos[:, 0]) 
    wy = np.array(atmos[:, 2])
    wx = np.array(atmos[:, 3])
   

    Ceff =  []
    wCs =   []
    for i in range(len(h)):
        Cx, Cy, Cz = (Cs[i] + wx[i])*x, (Cs[i] + wy[i])*y, Cs[i]*z
        Ceff.append([Cx, Cy, Cz])
        wCs.append(np.sqrt(Cx*Cx + Cy*Cy + Cz*Cz))
    
    wCs = np.array(wCs)
    Ceff = np.array(wCs)

    return wCs, wx, wy, Ceff


def doppler(source, stn, tau0, Cs, atmos, Z, s, alt):
    # %==========================================================================
    # % This function is used to calculate the Doppler shift
    # %
    # % Input variables:
    # %   source: location of the source [lat, long]
    # %   stn:    location of the station [lat, long]
    # %   tau0:   fundamental period
    # %   Cs:     ambient speed of sound
    # %   atmos:  atmospheric profile (txt file)
    # %   inc:    inclination angle at observing station
    # %   Z:      altitude intervals (line 108)
    # %   wtype:  type of the wind data file (component form = 0 or
    # % magnitude/direction = 1)
    # %
    # % Last modified on: 16-Oct-2013 (E. A. Silber)
    # %==========================================================================
    # % wind components
    # % s = 0/1 switch [a,b,c], where
    # % a = homogeneous (0) or inhomogeneous (1) atmosphere
    # % b = meridional, or y-direction winds on (1) or off (0)
    # % c = zonal, or x-direction winds on (1) or off (0)
    # %
    # % Z and tau0 are top-down, therefore they need to be flipped in order to
    # % match all other variables. 
    

    wCs, wx, wy, Ceff = windC(source, stn, atmos, Cs)

    n = len(Z)
    h = alt

    Z =     np.flipud(np.array(Z))
    h =     np.flipud(np.array(h)) 
    Cs =    np.flipud(np.array(Cs)) 
    Ceff =  np.flipud(np.array(Ceff)) 
    wCs =   np.flipud(np.array(wCs))
    wx =    np.flipud(np.array(wx))
    wy =    np.flipud(np.array(wy))

    f = scipy.interpolate.interp1d(h, wx)
    iwx = f(Z)

    f = scipy.interpolate.interp1d(h, wy)
    iwy = f(Z)
 
    wz_speed = 0.01;

    iwz = np.array([wz_speed]*5000)
  
    f = scipy.interpolate.interp1d(h, Cs)
    iCs = f(Z)
   

    wind_unit = []

    wind_speed = []

    for ii in range(n):
        step = np.sqrt(iwx[ii]**2 + iwy[ii]**2 + iwz[ii]**2)
        if step == 0:
    # %         disp('Wind speed is zero in all directions (x,y,z).');
    # %         disp('Doppler shift in frequency is not applicable.');
    # %         disp('No change in frequency or period was applied.');
    # %         disp('==================================================');
            return tau0
        wind_unit.append([iwx[ii]/step, iwy[ii]/step, iwz[ii]/step])
        wind_speed.append([iwx[ii], iwy[ii], iwz[ii]])

    wind_unit = np.flipud(wind_unit)


    sP = Position(*source)
    stP = Position(*stn)

    waveP = stP - sP
    x, y, z = waveP.x, waveP.y, waveP.z*1000
    mag = np.sqrt(x**2 + y**2 + z**2)

    x, y, z = x/mag, y/mag, z/mag

    f0 = 1/tau0
    angf0 = 2*np.pi*f0

    # %
    # %==========================================================================
    # % Perform Doppler shift calculations following Morse & Ingard (1968):
    # % OMEGA = angf - kw, where

    angf = angf0
    kn = []
    kn.append(angf/iCs[0])
    k = []
    k.append(kn[0]*np.array([x, y, z]))

    kw = []
    kw.append(np.dot(np.array(k)[0, :], np.array(wind_speed)[0, :]))
    OMEGA = []
    OMEGA.append(angf0 - kw[0])

    F = []
    F.append(OMEGA[0]/(2*np.pi))
    tau = []
    tau.append(tau0)
    dOMEGA = []

    for j in range(n):
        kn.append(angf/iCs[j])
        a = kn[j]*np.array([x, y, z])
        k.append(a)
        kw.append(np.dot(a, np.array(wind_speed)[j, :]))
        dOMEGA.append(kw[j])

    for i in range(1, n):
        OMEGA.append(angf0 - (1/j)*(np.sum(dOMEGA[:i])))
        F.append(OMEGA[i]/(2*np.pi))             
        tau.append(1/F[i])


    return tau                         




if __name__ == "__main__":

    source = [51.9183000000000, -2.43020000000000, 30.0400000000000]
    stat3 = [51.2112, -0.3298, 0.080]
    Ro = 7.3
    v = 13
    theta = 41.69
    dphi3 = 24.4958
    sw = [1, 0, 0]
    data3 = np.array([[30040,   19.2830354000000,   12.0092002400000,    1.61793559900000,    27.5531885400000],
[29811.1552100000,    19.2830354000000,    12.0092002400000 ,   1.61793559900000  ,  28.3203220600000],
[29345.3559100000,    18.9468873800000,    11.1503036000000 ,   0.873037995000000 ,  29.9483897000000],
[28879.5566100000,    18.6262410900000,    10.2095693000000 ,   -0.0203514030000000, 31.6700510600000],
[28413.7573100000,    18.3262035300000,    9.21982058000000 ,   -0.967960426000000,  33.4906866200000],
[27947.9580000000,    18.0518817000000,    8.21388069600000 ,   -1.87551690200000 ,  35.4159861700000],
[27482.1587000000,    17.8083826200000,    7.22457289300000 ,   -2.64874866000000 ,  37.4519665900000],
[27016.3594000000,    17.6008132700000,    6.28472042100000 ,   -3.19338353100000 ,  39.6049906500000],
[26550.5601000000,    17.4342806700000,    5.42714652700000 ,   -3.41514934200000 ,  41.8817869300000],
[26084.7608000000,    17.3138918200000,    4.68467445900000 ,   -3.21977392400000 ,  44.2894707800000],
[25618.9615000000,    17.2442143500000,    4.08849683000000 ,   -2.52282450700000 ,  46.8355666300000],
[25153.1622000000,    17.2207883200000,    3.64251399600000 ,   -1.40455221100000 ,  49.5280314400000],
[24687.3629000000,    17.2316798900000,    3.32803113100000 ,   -0.081549479000000,  52.3752796200000],
[24221.5636000000,    17.2647300500000,    3.12567270900000 ,   1.22548382500000  ,  55.3862092900000],
[23755.7643000000,    17.3077798100000,    3.01606320100000 ,   2.29584783400000  ,  58.5702301200000],
[23289.9650000000,    17.3487023300000,    2.97981198500000 ,   2.90943988400000  ,  61.9372926900000],
[22824.1657000000,    17.3809240800000,    2.99492175400000 ,   2.94928968100000  ,  65.4979196500000],
[22358.3664000000,    17.4096920500000,    3.03384678300000 ,   2.51794798300000  ,  69.2632385400000],
[21892.5671000000,    17.4417678600000,    3.06833037700000 ,   1.74609462700000  ,  73.2450166200000],
[21426.7678000000,    17.4839131500000,    3.07011584200000 ,   0.764409448000000 ,  77.4556976000000],
[20960.9685000000,    17.5428895400000,    3.01094648800000 ,   -0.296427716000000,  81.9084405700000],
[20495.1692000000,    17.6254586500000,    2.86256561800000 ,   -1.30573703000000 ,  86.6171611000000],
[20029.3699000000,    17.7383880900000,    2.59937613600000 ,   -2.13445816200000 ,  91.5965747200000],
[19563.5706000000,    17.8886045300000,    2.26666285000000 ,   -2.69669289200000 ,  96.8622429300000],
[19097.7713000000,    18.0832068700000,    1.98644819600000 ,   -2.95327082700000 ,  102.430621800000],
[18631.9720000000,    18.3293025800000,    1.88458766000000 ,   -2.86735563300000 ,  108.319113500000],
[18166.1727000000,    18.6339991800000,    2.08693672100000 ,   -2.40211098200000 ,  114.546120500000],
[17700.3734000000,    18.9981690500000,    2.65538233000000 ,   -1.55701180700000 ,  121.131103300000],
[17234.5741000000,    19.3981237300000,    3.39983053200000 ,   -0.474568141000000,  128.094641000000],
[16768.7748000000,    19.8041279600000,    4.06815067800000 ,   0.667495312000000 ,  135.458495800000],
[16302.9755000000,    20.1864464700000,    4.40821211900000 ,   1.69145384600000  ,  143.245681000000],
[15837.1762000000,    20.5154564500000,    4.16909953200000 ,   2.41982200300000  ,  151.480532900000],
[15371.3769000000,    20.7814372300000,    3.31498619800000 ,   2.71745592000000  ,  160.188786800000],
[14905.5776000000,    21.0172786000000,    2.27055029600000 ,   2.53986513900000  ,  169.397657400000],
[14439.7783000000,    21.2611940300000,    1.51878759200000 ,   1.85443707900000  ,  179.135923900000],
[13973.9790000000,    21.5152815300000,    1.31521930100000 ,   0.666478214000000 ,  189.434020100000],
[13508.1797000000,    21.7014032700000,    1.40999834600000 ,   -0.934462214000000,  200.324129300000],
[13042.3804000000,    21.7135357700000,    1.51361931800000 ,   -2.72511417500000 ,  211.840284800000],
[12576.5811000000,    21.2844879800000,    1.60441615000000 ,   -3.32758735900000 ,  224.018476600000],
[12110.7818000000,    20.1280569500000,    1.80671965900000 ,   -1.02816098400000 ,  236.896763500000],
[11644.9825000000,    18.8857245900000,    2.22991100100000 ,   1.90497940600000  ,  250.515392400000],
[11179.1832000000,    18.6974016900000,    2.82487362100000 ,   1.25106017700000  ,  264.916923700000],
[10713.3839000000,    19.6487717600000,    2.59380956400000 ,   -1.45927153400000 ,  280.146364600000],
[10247.5846000000,    21.2946382200000,    1.48804729700000 ,   -2.93359371300000 ,  296.251309700000],
[9781.78530200000,    23.0860407700000,    1.52501285400000 ,   -1.43283009000000 ,  313.282089600000],
[9315.98600200000,    24.9832396700000,    2.49512839000000 ,   1.56148398400000  ,  331.291928300000],
[8850.18670200000,    27.1508261300000,    3.22498799700000 ,   3.38554864600000  ,  350.337109600000],
[8384.38740100000,    29.6192824800000,    3.06454747800000 ,   3.15487528100000  ,  370.477152900000],
[7918.58810100000,    32.2113985100000,    2.17559126100000 ,   2.73472409900000  ,  391.774999200000],
[7452.78880100000,    34.7896644100000,    0.998437394000000,   3.19420349500000  ,  414.297207700000],
[6986.98950100000,    37.3711217300000,    0.535159967000000,   2.82788006700000  ,  438.114164200000],
[6521.19020100000,    39.9097360800000,    1.38029961300000 ,   0.744103870000000 ,  463.300300700000],
[6055.39090100000,    42.1443019200000,    2.87548114400000 ,   -0.788679054000000,  489.934327900000],
[5589.59160100000,    44.1403914000000,    4.27619715900000 ,   -0.782552838000000,  518.099481700000],
[5123.79230100000,    46.4009971200000,    5.06369378700000 ,   -0.668129721000000,  547.883782900000],
[4657.99300100000,    48.6647982500000,    5.43402403400000 ,   -0.779143340000000,  579.380312400000],
[4192.19370100000,    50.5624311700000,    5.61620082300000 ,   -0.860071698000000,  612.687502100000],
[3726.39440100000,    52.4474632000000,    5.31358619700000 ,   -0.964503001000000,  647.909442500000],
[3260.59510100000,    54.2162936500000,    6.51736006400000 ,   -1.05726338200000 ,  685.156208200000],
[2794.79580000000,    55.9922729000000,    9.21998439200000 ,   -1.11589134500000 ,  724.544201400000],
[2328.99650000000,    58.0104813700000,    10.6576076600000 ,   -1.16115673900000 ,  766.196516200000],
[1863.19720000000,    59.6992087800000,    10.7914352100000 ,   -1.18021050200000 ,  810.243323100000],
[1397.39790000000,    60.8453442000000,    10.9822039700000 ,   -1.30434323700000 ,  856.822275700000],
[931.598600200000,    61.7477812900000,    12.5594856100000 ,   -1.44036513500000 ,  906.078940900000],
[465.799300100000,    62.1792358000000,    13.6274798500000 ,   -1.48388871000000 ,  958.167254100000],
[80,  60.7438690000000,    10.9109203900000,    -1.78912816400000,   1003.56934200000]])

    import time

    t1 = time.time()
    results = overpressureihmod_Ro(source,stat3,Ro,v,theta,dphi3,data3,sw, wind=False)
    t2 = time.time()
    print("Overpressure (base): {:.4f} s".format(t2 - t1))