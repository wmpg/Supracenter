import numpy as np
import matplotlib.pyplot as plt
import multiprocessing
import warnings
from functools import partial

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.Utils.pso import pso
from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import createLabelEditObj
from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT
from supra.Utils.TryObj import *
from supra.Utils.EigenRay import eigenConsistancy
from supra.Utils.Energy import EnergyObj

#from supra.Supracenter.cyscanIntegration import cyscan as intscan
from supra.Supracenter.anglescanYieldInt import anglescan as intscan
from supra.Supracenter.anglescan2 import anglescan
from supra.Supracenter.cyscan5 import cyscan
from supra.Utils.Classes import Position, Constants, Color

from supra.Utils.AngleConv import percDiff

from supra.Yields.YieldCalcs import *

from supra.Atmosphere.Parse import parseWeather

from supra.Supracenter.refractiveFactor import refractiveFactor

consts = Constants()
ColorGen = Color()

N_LAYERS = 150

def ReedYield(R, dp):
    """ Emperical Relation from Reed 1977
    R - metres
    pressure - Pa
    yield - given in Joules for Chemical Explosives
    """

    # reed_yield = 2100*(R/1000)**3*(dp/1000)**(5/2)*4.184e6
    reed_yield = 2.77e-7*dp**(5/2)*R**3
    return reed_yield

def ANSIYield(R, P, dp, P_0):

    # ansi_yield = (dp/105.93/1000*R**(1.1)/(P/101325)**(0.6333))**(1/0.3667)
    ansi_yield = (dp*1.35/105.93/1000)**(3/1.1)*R**(3)/(P/P_0)**(3/1.1 - 1)
    #ansi_yield = 217.32*dp**(2.727)*R**3*(P/101325)**(-1.727)
    return ansi_yield


def chemFuncMinimizer(x, *args):
    del_P = args

    Z = 10**x

    P_est = chem_func(Z)

    err = np.abs(del_P - P_est)

    return err

def nucFuncMinimizer(x, *args):
    del_P = args

    Z = 10**x

    P_est = nuc_func(Z)

    err = np.abs(del_P - P_est)

    return err

def chemFuncDurationMinimizer(x, *args):
    Jd, R, f_d = args

    Z = 10**x

    Jd_est = chemDuration(Z, R, f_d)

    err = np.abs(Jd - Jd_est)

    return err

def nucFuncDurationMinimizer(x, *args):
    Jd, R, f_d = args

    Z = 10**x

    Jd_est = nucDuration(Z, R, f_d)

    err = np.abs(Jd - Jd_est)

    return err

def chemFuncImpulseMinimizer(x, *args):
    I_A = args

    Z = 10**x

    I_A_est = chemImpulse(Z)

    err = np.abs(I_A - I_A_est)

    return err

def findScaledDistance(x, func):
    # Use magnitude for faster search
    a, b = pso(func, [-2], [7], args=x, maxiter=1000, swarmsize=1000)

    scaled_distance = 10**a[0]
    return scaled_distance

def nucDuration(Z, R, f_d):
    # Returns time in seconds
    # return (f_d*R/Z)*(4.184e12)**(1/3)*180*(1 + (Z/100)**3)**(0.5)/(1 + Z/40)**(0.5)/(1 + (Z/285)**5)**(1/6)/(1 + Z/50000)**(1/6)
    return ((f_d*R/Z)*(4.184e12)**(1/3)*180*(1 + (Z/100)**3)**(0.5)/(1 + Z/40)**(0.5)/(1 + (Z/285)**5)**(1/6)/(1 + Z/50000)**(1/6))/1000

def chemDuration(Z, R, f_d):
    # return (f_d*R/Z)*(4.184e6)**(1/3)*980*(1 + (Z/0.54)**10)/(1 + (Z/0.02)**3)/(1 + (Z/0.74)**6)/(1 + (Z/6.9)**2)**(0.5)
    return (f_d*R/Z)*(4.184e6)**(1/3)*980*(1 + (Z/0.54)**10)/(1 + (Z/0.02)**3)/(1 + (Z/0.74)**6)/(1 + (Z/6.9)**2)**(0.5)

def nuc_func(Z):
    # KG85
    return 3.2e6*Z**(-3)*(1 + (Z/87)**2)**(0.5)*(1+Z/800)

def chem_func(Z):
    ''' Chemical explosion function given in KG85 Eq 6-2

    '''
    p_Pa = 808*(1 + (Z/4.5)**2)/(1 + (Z/0.048)**2)**0.5/(1 + (Z/0.32)**2)**0.5/(1 + (Z/1.35)**2)**0.5
    return p_Pa


def chemImpulse(Z):

    return (0.067*(1 + (Z/0.23)**4)**(0.5))/Z**2/(1 + (Z/1.55)**3)**(1/3)

def integration_full(k, v, b, I, P_a):
    return np.exp(-k*v**2/b/P_a*I)

def scaledDistance(f_d, R, W, W_0):
    ''' Scaled distance from KG85
    Arguments:
    f_d [float] transmission factor
    R [float] Range in meters
    W [float] Actual blast yield
    W_0 [float] Reference blast yield

    Returns:
    Z [float] Scaled distance in m/kg^(1/3)
    '''

    #scaled distance
    Z = f_d*R/(W/W_0)**(1/3)
    
    return Z

def gunc_error(W, *args):
    p_ratio, J_m, W_0, P, P_0, c, c_m, f_d, R, P_a, k, b, I, cf, h_R, v = args
    overpressure = total_func(W[0], J_m, W_0, P, P_0, c, c_m, f_d, R, P_a, k, b, I, cf, h_R, v)
    return np.abs(overpressure - p_ratio)

def total_func(W, J_m, W_0, P, P_0, c, c_m, f_d, R, P_a, k, b, I, cf, h_R, v):

    W = 10**W
    # v = 1/2/J_m*(W_0*P/W/P_0)**(1/3)*(c/c_m)
    overpressure = chem_func(scaledDistance(f_d, R, W, W_0))*P_a*(I)*cf
    return overpressure


class Yield(QWidget):
    ''' Dialog to estimate yields from overpressures
    '''
    def __init__(self, bam, prefs, current_station):

        #################
        # Initialize GUI
        #################
        QWidget.__init__(self)
        self.setWindowTitle('Yields')
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        self.prefs = prefs
        self.bam = bam
        self.setup = bam.setup
        self.stn_list = bam.stn_list
        self.current_station = current_station
        self.iterator = 0

        theme(self)

        self.count = 0
        layout = QHBoxLayout()
        self.setLayout(layout)

        pane1 = QGridLayout()
        layout.addLayout(pane1)

        pane2 = QVBoxLayout()
        layout.addLayout(pane2)

        self.station_label = QLabel('Station: {:}'.format(self.stn_list[self.current_station].metadata.code))
        pane1.addWidget(self.station_label, 0, 1, 1, 1)

        self.station1_label = QLabel('Nominal')
        pane1.addWidget(self.station1_label, 0, 2, 1, 1)


        self.height_label, self.height_edits = createLabelEditObj('Height', pane1, 1)
        self.range_label, self.range_edits = createLabelEditObj('Range', pane1, 2)
        self.pressure_label, self.pressure_edits = createLabelEditObj('Explosion Height Pressure', pane1, 3)
        self.overpressure_label, self.overpressure_edits = createLabelEditObj('Overpressure', pane1, 4)
        self.afi_label, self.afi_edits = createLabelEditObj('Attenuation Integration Factor', pane1, 5)
        self.geo_label, self.geo_edits = createLabelEditObj('Geometric Factor', pane1, 6)
        self.p_a_label, self.p_a_edits = createLabelEditObj('Ambient Pressure', pane1, 7)
        self.Jd_label, self.Jd_edits = createLabelEditObj('Positive Phase Length [ms]', pane1, 8)
        self.fd_label, self.fd_edits = createLabelEditObj('Transmission Factor', pane1, 9)
        _, self.freq_edits = createLabelEditObj("Dominant Period", pane1, 10)
        _, self.rfangle_edits = createLabelEditObj("Refraction Angle", pane1, 11)
        # _, self.I_A_edits = createLabelEditObj("Impulse per Area (Area under curve)", pane1, 12)



        # self.fyield_button = QPushButton('Calculate Yield (Frequency)')
        # pane1.addWidget(self.fyield_button, 15, 1, 1, 4)
        # self.fyield_button.clicked.connect(self.freqYield)

        self.yield_button = QPushButton('Calculate Yield')
        pane1.addWidget(self.yield_button, 14, 1, 1, 4)
        self.yield_button.clicked.connect(self.yieldCalc)

        self.integrate_button = QPushButton('Integrate')
        pane1.addWidget(self.integrate_button, 13, 1, 1, 4)
        self.integrate_button.clicked.connect(self.intCalc)

        # Constants - Reed 1972
        self.W_0 = 4.184e6 # Standard reference explosion yield (1 kg of Chemical TNT)
        self.P_0 = 101325 # Standard pressure
        self.b = 1/consts.SCALE_HEIGHT
        self.k = consts.ABSORPTION_COEFF
        self.J_m = 0.375  # Avg positive period of reference explosion
        self.c_m = 347    # Sound speed of reference explosion

        # self.blastline_view = pg.GraphicsLayoutWidget()
        # self.blastline_canvas = self.blastline_view.addPlot()
        # pane2.addWidget(self.blastline_view)
        # self.blastline_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)


        self.blastline_plot = MatplotlibPyQT()
        self.blastline_plot.ax = self.blastline_plot.figure.add_subplot(111)
        pane2.addWidget(self.blastline_plot)

    def freqYield(self):
        print("Broken")
        # self.R = tryFloat(self.range_edits.text())
        # self.P_a = tryFloat(self.p_a_edits.text())
        # self.cf = tryFloat(self.geo_edits.text())
        # self.I = tryFloat(self.afi_edits.text())
        # self.P = tryFloat(self.pressure_edits.text())
        # self.c = tryFloat(self.c_edits.text())
        # self.f_d = tryFloat(self.fd_edits.text())
        # self.T_d = tryFloat(self.freq_edits.text())


        # #W_0 = 4.184e12 (IBM Problem)
        # W = self.W_0*self.P_a/self.P_0*(self.c*self.T_d/2/self.c_m/self.J_m)**3
        # print("Freq: {:.2f} -> {:.2E} J".format(1/self.T_d, W))

    def integration_full(self, k, v, b, I, P_a):
        return np.exp(-I*(k*v**2/b)*N_LAYERS)
        # return np.exp(-k*v**2/b/P_a*I)

  
    def inverse_gunc(self, p_ans):
        a, b = pso(gunc_error, [1], [15], args=([p_ans, self.J_m, \
                        self.W_0, self.P, self.P_0, self.Jd, self.c_m, self.f_d, self.R,\
                         self.P_a, self.k, self.b, self.I, self.cf, self.horizontal_range, self.v]), \
                        processes=1, swarmsize=1000, maxiter=1000)

        return 10**(a[0])

    def integrate(self, height, D_ANGLE=1.5, tf=1, az=1):
        ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)
        try:
            point = self.setup.trajectory.findGeo(height)
        except AttributeError:
            print("STATUS: No trajectory given, assuming lat/lon center")
            point = Position(self.setup.lat_centre, self.setup.lon_centre, height)
        point.pos_loc(ref_pos)

        stn = self.stn_list[self.current_station]
        stn.metadata.position.pos_loc(ref_pos)

        lats = [point.lat, stn.metadata.position.lat]
        lons = [point.lon, stn.metadata.position.lon]
        elevs = [point.elev, stn.metadata.position.elev]

        # make the spline lower to save time here


        sounding, perturbations = self.bam.atmos.getSounding(lats, lons, elevs, spline=N_LAYERS, ref_time=self.setup.fireball_datetime)



        trans = []
        ints = []
        ts = []
        ps = []
        rfs = []



        if perturbations is None:
            ptb_len = 1
        else:
            ptb_len = len(perturbations) + 1


        for ptb_n in range(ptb_len):

            # Temporary adjustment to remove randomness from perts
        
            if ptb_n == 0:
                zProfile = sounding
            else:
                zProfile = perturbations[ptb_n - 1]

            S = np.array([point.x, point.y, point.z])
            D = np.array([stn.metadata.position.x, stn.metadata.position.y, stn.metadata.position.z])

            _, az_n, tf_n, _ = cyscan(S, D, zProfile, wind=True,\
                h_tol=30, v_tol=30)

            self.T_d = tryFloat(self.freq_edits.text())

            self.v = 1/self.T_d

            f, g, T, P_2, P_1, path_length, pdr, reed_attenuation = intscan(S, az_n, tf_n, zProfile, self.v, wind=True)

            self.reed_attenuation = reed_attenuation

            rf = refractiveFactor(point, stn.metadata.position, zProfile, D_ANGLE=D_ANGLE)
            trans.append(f)
            ints.append(g)
            ts.append(T)
            ps = [P_2, P_1]

            rfs.append(rf)


        return trans, ints, ts, ps, rfs, path_length, pdr

    def intCalc(self):

        stn = self.stn_list[self.current_station]
        if tryFloat(self.height_edits.text()) != None:
            height = tryFloat(self.height_edits.text())
            self.rfang = tryFloat(self.rfangle_edits.text())
            trans, ints, ts, ps, rfs, path_length, pdr = self.integrate(height, D_ANGLE=self.rfang)

            f_val = np.nanmean(trans)
            g_val = np.nanmean(ints)
            t_val = np.nanmean(ts)
            r_val = np.nanmean(rfs)

            

            # h = np.linspace(20000, 34000, 56)
            # d_ang = np.array([1.5, 2.0])
            # c = ['w', 'm', 'r', 'b', 'g']
            # print('Code Started')
            # for ii, d in enumerate(d_ang):
            #     my_data = []
            #     for height in h:
            #         trans, ints, ts, ps, rfs = self.integrate(height, D_ANGLE=d)

            #         f_val = np.nanmean(trans)
            #         g_val = np.nanmean(ints)
            #         t_val = np.nanmean(ts)
            #         p_val = np.nanmean(ps)
            #         r_val = np.nanmean(rfs)
            #         my_data.append(r_val)
            #         print('RF: {:} | ANGLE: {:} | HEIGHT: {:}'.format(r_val, d, height))

            #     plt.scatter(h, my_data, label='Angle: {:} deg'.format(d), c=c[ii])
            # plt.legend()
            # plt.show()
            # print('RF - Not checking if consistant')



            self.fd_edits.setText('{:.4f}'.format(f_val))
            self.afi_edits.setText('{:.4f}'.format(g_val))
            # self.c_edits.setText('{:.4f}'.format(t_val))
            self.pressure_edits.setText('{:.4f}'.format(ps[0]))
            self.p_a_edits.setText('{:.4f}'.format(ps[1]))
            self.path_length = path_length
            self.pdr = pdr

            try:
                frag_pos = self.setup.trajectory.findGeo(height)
            except AttributeError:
                print("STATUS: No trajectory given, assuming lat/lon center")
                frag_pos = Position(self.setup.lat_centre, self.setup.lon_centre, height)
            
            self.geo_edits.setText('{:.4f}'.format(r_val))
            stn_pos = stn.metadata.position
            dist = stn_pos.pos_distance(frag_pos)
            self.range_edits.setText('{:.4f}'.format(dist))

        # if tryFloat(self.height_min_edits.text()) != None:
        #     height = tryFloat(self.height_min_edits.text())

        #     trans, ints, ts, ps, rfs = self.integrate(height)

        #     f_val = np.nanmean(trans)
        #     g_val = np.nanmean(ints)
        #     t_val = np.nanmean(ts)
        #     p_val = np.nanmean(ps)
        #     r_val = np.nanmean(rfs)

        #     self.fd_min_edits.setText('{:.4f}'.format(f_val))
        #     self.afi_min_edits.setText('{:.4f}'.format(g_val))
        #     self.c_min_edits.setText('{:.4f}'.format(t_val))
        #     self.pressure_min_edits.setText('{:.4f}'.format(p_val))
        #     self.p_a_min_edits.setText('{:.4f}'.format(estPressure(stn.metadata.position.elev)))
        #     frag_pos = self.setup.trajectory.findGeo(height)
        #     self.geo_min_edits.setText('{:.4f}'.format(r_val))
        #     stn_pos = stn.metadata.position
        #     dist = stn_pos.pos_distance(frag_pos)
        #     self.range_min_edits.setText('{:.4f}'.format(dist))

        # if tryFloat(self.height_max_edits.text()) != None:
        #     height = tryFloat(self.height_max_edits.text())

        #     trans, ints, ts, ps, rfs = self.integrate(height)

        #     f_val = np.nanmean(trans)
        #     g_val = np.nanmean(ints)
        #     t_val = np.nanmean(ts)
        #     p_val = np.nanmean(ps)
        #     r_val = np.nanmean(rfs)

        #     self.fd_max_edits.setText('{:.4f}'.format(f_val))
        #     self.afi_max_edits.setText('{:.4f}'.format(g_val))
        #     self.c_max_edits.setText('{:.4f}'.format(t_val))
        #     self.pressure_max_edits.setText('{:.4f}'.format(p_val))
        #     self.p_a_max_edits.setText('{:.4f}'.format(estPressure(stn.metadata.position.elev)))
        #     frag_pos = self.setup.trajectory.findGeo(height)
        #     self.geo_max_edits.setText('{:.4f}'.format(r_val))
        #     stn_pos = stn.metadata.position
        #     dist = stn_pos.pos_distance(frag_pos)
        #     self.range_max_edits.setText('{:.4f}'.format(dist))


    def yieldCalc(self):
        w = np.linspace(np.log10(1e1), np.log10(1e15))
        W = 10**w

        if tryFloat(self.height_edits.text()) != None:
            self.R = tryFloat(self.range_edits.text())
            self.P_a = tryFloat(self.p_a_edits.text())
            self.cf = tryFloat(self.geo_edits.text())
            self.I = tryFloat(self.afi_edits.text())
            self.P = tryFloat(self.pressure_edits.text())
            self.Jd = tryFloat(self.Jd_edits.text())
            self.f_d = tryFloat(self.fd_edits.text())
            # self.I_A = tryFloat(self.I_A_edits.text())
           


            # Kinney and Graham 1985 uses horizontal distance (see Problem 7.1 in their book)
            height = tryFloat(self.height_edits.text())
            frag_pos = self.setup.trajectory.findGeo(height)
            stn = self.stn_list[self.current_station]
            stn_pos = stn.metadata.position
            dist = stn_pos.ground_distance(frag_pos)
            true_range = stn_pos.pos_distance(frag_pos)
            self.horizontal_range = dist



            print("Height:          {:.2f} km".format(height/1000))
            print("Ground Distance: {:.2f} km".format(self.horizontal_range/1000))
            print("Slant Range:     {:.2f} km".format(true_range/1000))
            print("Path Length:     {:.2f} km".format(self.path_length/1000))
            print("")
            print("Positive Phase Duration  {:.0f} ms".format(self.Jd))

            del_P = tryFloat(self.overpressure_edits.text())

            print("Overpressure:            {:.2f} Pa".format(del_P))
            # print("Impulse / Area           {:.2f} Pa*s".format(self.I_A))
            print("Dominant Period:         {:.2f} s".format(self.T_d))

            perc_diff = percDiff(self.T_d, self.Jd/1000*2)

            print("Dominant Period and Positive Phase are different by: {:.2f}%".format(perc_diff))
            print("")

            del_P_adj = del_P

            print("### Pressure Changes")
            print("\tAttenuation Factor: {:.2f}".format(self.I))
            print("\tGeometric Factor:   {:.2f}".format(self.cf))
            print("\tCombined Factor:    {:.2f}".format(self.I*self.cf))
            print("")
            print("### Atmospheric Pressure")
            print("\tAt Source      {:8.2f} Pa".format(self.P))
            print("\tAt Reciever    {:8.2f} Pa".format(self.P_a))
            print("\tPressure Ratio {:8.2f} Pa".format(self.P/self.P_a))
            print("")
            print("Overpressure (Reciever)                      : {:.2E} Pa".format(del_P_adj))
            print("Overpressure (Reciever at Source Pressure)   : {:.2E} Pa".format(del_P_adj*(self.P_a/self.P)**(1/6)))
            print("Overpressure Ratio (Source)                  : {:.2E}".format(del_P_adj/self.P))
            print("Overpressure Ratio (Reciever)                : {:.2E}".format(del_P_adj/self.P_a))

            print("")

            #needham blast waves
            # brode
            # jones 1968

            # Reed 1972 enhancement for airborne bursts
            pres_fact = (self.P_a/self.P)**(1/6)
            pres_fact_85 = 1/self.f_d


            # f_t depends on sound speed - roughly a factor of 1
            f_t = self.f_d*(330/330)

            p_p0 = np.exp(self.reed_attenuation)

            print("Pressure Factor (Reed 1972) Attenuation: {:.2f}".format(p_p0))
            print("Pressure Factor (Reed 1972) 1/6 Power  : {:.2f}".format(pres_fact))
            print("Pressure Factor (KG   1985)            : {:.2f}".format(pres_fact_85))



            new_pres = del_P_adj*pres_fact/self.I/self.cf
            phase_duration = self.Jd*(1/f_t)


            # reed_yield = ReedYield(self.R, new_pres/pres_fact)
            # ansi_yield = ANSIYield(self.R, self.P, new_pres/pres_fact, self.P_a)
            

            print("Height Adjusted Pressure: {:.2f} Pa".format(new_pres))
            print("Height Adjusted Time    : {:.2f} ms".format(phase_duration))

            Z_chem = findScaledDistance([(new_pres)/(self.P_a)], chemFuncMinimizer)
            Z_nuc  = findScaledDistance([(new_pres)/(self.P_a)], nucFuncMinimizer)
            Z_chem_d = findScaledDistance([phase_duration, self.R, self.f_d], chemFuncDurationMinimizer)
            Z_nuc_d  = findScaledDistance([phase_duration, self.R, self.f_d], nucFuncDurationMinimizer)

            # # Z_chem_IA = findScaledDistance([self.I_A], chemFuncImpulseMinimizer)

            print("Scaled Distance from Overpressure (Chemical):    {:10.2f} km".format(Z_chem/1000))
            print("Scaled Distance from Duration (Chemical):        {:10.2f} km".format(Z_chem_d/1000))
            # print("Scaled Distance from Impulse (Chemical):         {:10.2f} km".format(Z_chem_IA/1000))
            print("Scaled Distance from Overpressure  (Nuclear):    {:10.2f} km".format(Z_nuc/1000))
            print("Scaled Distance from Duration  (Nuclear):        {:10.2f} km".format(Z_nuc_d/1000))

            #Yield_chem = overpressure2YieldKG(del_P, height, stn_pos.elev, true_range)
            Yield_chem = overpressure2YieldKG(new_pres, height, stn_pos.elev, true_range)
            Yield_nuc = (self.f_d*self.R/Z_nuc)**3*self.W_0*1e6

            Yield_chem_d = (self.f_d*self.R/Z_chem_d)**3*self.W_0
            Yield_nuc_d = (self.f_d*self.R/Z_nuc_d)**3*self.W_0*1e6
            # Yield_chem_IA = (self.f_d*self.R/Z_chem_IA)**3*self.W_0
            # print("### Previous Yields")
            # print("\tExpected Yield (ANSI 1977 Empirical):    {:.2E} J ({:.2f} kg TNT)".format(ansi_yield, ansi_yield/self.W_0))
            # print("\tExpected Yield (Reed 1977 Empirical):    {:.2E} J ({:.2f} kg TNT)".format(reed_yield, reed_yield/self.W_0))
            # print("\tExpected Yield (Chemical Overpressure):  {:.2E} J ({:.2f} kg TNT)".format(Yield_chem, Yield_chem/self.W_0))
            # print("\tExpected Yield (Chemical Duration):      {:.2E} J ({:.2f} kg TNT)".format(Yield_chem_d, Yield_chem_d/self.W_0))
            # # print("Expected Yield (Chemical Impulse):       {:.2E} J ({:.2f} kg TNT)".format(Yield_chem_IA, Yield_chem_IA/self.W_0))
            # print("\tExpected Yield (Nuclear Overpressure):   {:.2E} J ({:.2f} kT TNT)".format(Yield_nuc, Yield_nuc/self.W_0/1e6))
            # print("\tExpected Yield (Nuclear Duration):       {:.2E} J ({:.2f} kT TNT)".format(Yield_nuc_d, Yield_nuc_d/self.W_0/1e6))

            # factor = 101325/self.P # This is 101325 because that is what the KG85 equations are reference to

            # print("Yield Correction Factor: {:.2f}".format(factor))
            # print("New Factor:              {:.2f}".format((self.pdr/self.P_a/self.path_length**3)))


            # print("\tExpected Yield (Chemical Overpressure):  {:.2E} J ({:.2f} kg TNT)".format(Yield_chem, Yield_chem/self.W_0))
            # print("\tExpected Yield (Chemical Duration):      {:.2E} J ({:.2f} kg TNT)".format(factor*Yield_chem_d, factor*Yield_chem_d/self.W_0))
            # # print("Expected Yield (Chemical Impulse):       {:.2E} J ({:.2f} kg TNT)".format(factor*Yield_chem_IA, factor*Yield_chem_IA/self.W_0))
            # print("\tExpected Yield (Nuclear Overpressure):   {:.2E} J ({:.2f} kT TNT)".format(factor*Yield_nuc, factor*Yield_nuc/self.W_0/1e6))
            # print("\tExpected Yield (Nuclear Duration):       {:.2E} J ({:.2f} kT TNT)".format(factor*Yield_nuc_d, factor*Yield_nuc_d/self.W_0/1e6))

            # print("### Sach Scaling")
            # Z = np.array([Z_chem, Z_nuc, Z_chem_d, Z_nuc_d])
            # Yield = (self.W_0/self.P_a/Z**3)*self.pdr
            # print("\tExpected Yield (Chemical Overpressure):  {:.2E} J ({:.2f} kg TNT)".format(Yield[0], Yield[0]/self.W_0))
            # print("\tExpected Yield (Chemical Duration):      {:.2E} J ({:.2f} kg TNT)".format(Yield[2], Yield[2]/self.W_0))
            # # print("Expected Yield (Chemical Impulse):       {:.2E} J ({:.2f} kg TNT)".format(factor*Yield_chem_IA, factor*Yield_chem_IA/self.W_0))
            # print("\tExpected Yield (Nuclear Overpressure):   {:.2E} J ({:.2f} kT TNT)".format(Yield[1]*1e6, Yield[1]/self.W_0))
            # print("\tExpected Yield (Nuclear Duration):       {:.2E} J ({:.2f} kT TNT)".format(Yield[3]*1e6, Yield[3]/self.W_0))

            # # Plot of overpressure vs. yield
            # del_p = chem_func(scaledDistance(self.f_d, self.R, W, self.W_0))*self.P_a*\
            #             self.I*self.cf


            # self.yieldPlot(del_p, W)

            ### PLOTTING

            sample_R = np.linspace(0, 100000)
            sample_dP = new_pres/pres_fact # Pa
            sample_P = self.P

            sample_W = ReedYield(sample_R, sample_dP)
            sample_W_ansi = ANSIYield(sample_R, sample_P, sample_dP, self.P_a)

            c = ["w", "r", "m"]


            self.blastline_plot.ax.plot(sample_R/1000,  sample_W_ansi/self.W_0,         c=c[self.iterator], linestyle="--", label="(ANSI 1983) Overpressure = {:.2f} Pa".format(sample_dP))
            self.blastline_plot.ax.plot(sample_R/1000,  sample_W/self.W_0,              c=c[self.iterator], label="(Reed 1972) Overpressure = {:.2f} Pa".format(sample_dP))
            self.blastline_plot.ax.scatter(self.R/1000, Yield_chem/self.W_0,            c=c[self.iterator], marker="<", label="Chemical Overpressure")
            self.blastline_plot.ax.scatter(self.R/1000, Yield_chem_d/self.W_0,          c=c[self.iterator], marker=">", label="Chemical Duration")
            self.blastline_plot.ax.legend()

            self.iterator += 1 

            self.blastline_plot.ax.semilogy()

            self.blastline_plot.ax.set_xlabel("Range [km]")
            self.blastline_plot.ax.set_ylabel("Yield [kg TNT HE]")

            self.blastline_plot.show()


            qm = QMessageBox()
            ret = qm.question(self, '', "Save this yield?", qm.Yes | qm.No)

            if ret == qm.Yes:
                a = EnergyObj()
                a.source_type = "Fragmentation"
                a.height = height
                a.range = self.R
                a.station = self.stn_list[self.current_station]
                a.ansi_yield = ansi_yield
                a.reed_yield = reed_yield 
                a.chem_pres = Yield_chem 
                a.chem_dur = Yield_chem_d
                a.nuc_pres = Yield_nuc
                a.nuc_dur = Yield_nuc_d
                a.adj_pres = new_pres
                a.adj_dur = phase_duration

                self.bam.energy_measurements.append(a)


    def yieldPlot(self, p_ratio, W, unc='none'):
        if unc == 'none':
            self.count += 1
        
        colour = [(0, 255, 26), (3, 252, 176), (252, 3, 3), (176, 252, 3), (255, 133, 3),
                      (149, 0, 255), (76, 128, 4), (82, 27, 27), (101, 128, 125), (5, 176, 249)]
        ptb_colour = [(0, 255, 26, 150), (3, 252, 176, 150), (252, 3, 3, 150), (176, 252, 3, 150), (255, 133, 3, 150),
                      (149, 0, 255, 150), (76, 128, 4, 150), (82, 27, 27, 150), (101, 128, 125, 150), (5, 176, 249, 150)]

        self.blastline_canvas.setLogMode(False, True)
        if unc == 'none':
            self.nominal = pg.PlotCurveItem(x=p_ratio, y=np.log10(W), pen=colour[self.count-1], name='Fragmentation {:}'.format(self.count))
            self.blastline_canvas.addItem(self.nominal)
            self.p_rat = tryFloat(self.overpressure_edits.text())
        # if unc == 'min':
        #     self.min = pg.PlotCurveItem(x=p_ratio, y=np.log10(W), pen=ptb_colour[self.count-1], name='Fragmentation {:}'.format(self.count))
        #     self.blastline_canvas.addItem(self.min)
        #     self.p_rat = tryFloat(self.overpressure_min_edits.text())
        # if unc == 'max':
        #     self.max = pg.PlotCurveItem(x=p_ratio, y=np.log10(W), pen=ptb_colour[self.count-1], name='Fragmentation {:}'.format(self.count))
        #     self.blastline_canvas.addItem(self.max)
        #     self.p_rat = tryFloat(self.overpressure_max_edits.text())
        # try:
        #     pfill = pg.FillBetweenItem(self.min, self.max, brush=ptb_colour[self.count-1])
        #     self.blastline_canvas.addItem(pfill)
        #     self.min = None
        #     self.max = None
        # except:
        #     pass
            

        self.blastline_canvas.setLabel('bottom', "Overpressure", units='Pa')
        self.blastline_canvas.setLabel('left', "Yield", units='J')
        # self.blastline_canvas.addItem(pg.InfiniteLine(pos=(self.p_rat, 0), angle=90))

        W_final = self.inverse_gunc(self.p_rat)
        self.blastline_canvas.scatterPlot(x=[self.p_rat], y=[W_final], pen=(255, 255, 255), brush=(255, 255, 255), size=10, pxMode=True)

        print('Blast Estimate     {:.2E} J'.format(W_final))

        if unc == 'none':
            txt = pg.TextItem("Frag {:}: {:.2E} J".format(self.count, W_final))
            txt.setPos(self.p_rat, np.log10(W_final))
            self.blastline_canvas.addItem(txt)
        self.blastline_canvas.setTitle('Fragmentation Yield Curves for a Given Overpressure')


if __name__ == "__main__":

    W = 1e6
    W_0 = 1e6

    Z = np.logspace(-1, 2)

    R = Z*(W/W_0)**(1/3)

    plt.loglog(Z, chem_func(Z))
    plt.scatter([10, 3, 0.3], [0.1, 1, 100])
    plt.show()

    print(R)

