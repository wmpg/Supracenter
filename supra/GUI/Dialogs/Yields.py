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
from supra.Utils.TryObj import *
from supra.Utils.EigenRay import eigenConsistancy
from supra.Supracenter.cyscanIntegration import cyscan as intscan
from supra.Supracenter.anglescan import anglescan
from supra.Supracenter.cyscan2 import cyscan
from supra.Utils.Classes import Position
from supra.Supracenter.l137 import estPressure

from supra.Atmosphere.Parse import parseWeather

def nuc_func(Z):
    # KG85
    return 3.2e6*Z**(-3)*(1 + (Z/87)**2)**(0.5)*(1+Z/800)

def chem_func(Z):
    return 808*(1 + (Z/4.5)**2)/(1 + (Z/0.048)**2)**0.5/(1 + (Z/0.32)**2)**0.5/(1 + (Z/1.35)**2)**0.5

def integration_full(k, v, b, I, P_a):
    return np.exp(-k*v**2/b/P_a*I)

def phi(f_d, d, W, W_0):
    #scaled distance
    return f_d*d/(W/W_0)**(1/3)

def gunc_error(W, *args):
    p_ratio, J_m, W_0, P, P_0, c, c_m, f_d, R, P_a, k, b, I, cf = args
    return np.abs(total_func(W[0], J_m, W_0, P, P_0, c, c_m, f_d, R, P_a, k, b, I, cf) - p_ratio)

def total_func(W, J_m, W_0, P, P_0, c, c_m, f_d, R, P_a, k, b, I, cf):
    W = 10**W
    v = 1/2/J_m*(W_0*P/W/P_0)**(1/3)*(c/c_m)
    return chem_func(phi(f_d, R, W, W_0))*P_a*(integration_full(k, v, b, I, P_0))*cf


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

        self.station2_label = QLabel('Min')
        pane1.addWidget(self.station2_label, 0, 3, 1, 1)

        self.station3_label = QLabel('Max')
        pane1.addWidget(self.station3_label, 0, 4, 1, 1)

        self.height_label, self.height_edits = createLabelEditObj('Height', pane1, 1)
        self.range_label, self.range_edits = createLabelEditObj('Range', pane1, 2)
        self.pressure_label, self.pressure_edits = createLabelEditObj('Explosion Height Pressure', pane1, 3)
        self.overpressure_label, self.overpressure_edits = createLabelEditObj('Overpressure', pane1, 4)
        self.afi_label, self.afi_edits = createLabelEditObj('Attenuation Integration Factor', pane1, 5)
        self.geo_label, self.geo_edits = createLabelEditObj('Geometric Factor', pane1, 6)
        self.p_a_label, self.p_a_edits = createLabelEditObj('Ambient Pressure', pane1, 7)
        self.c_label, self.c_edits = createLabelEditObj('Speed of Sound', pane1, 8)
        self.fd_label, self.fd_edits = createLabelEditObj('Transmission Factor', pane1, 9)

        self.height_min_edits = QLineEdit('')
        pane1.addWidget(self.height_min_edits, 1, 3, 1, 1)
        self.range_min_edits = QLineEdit('')
        pane1.addWidget(self.range_min_edits, 2, 3, 1, 1)
        self.pressure_min_edits = QLineEdit('')
        pane1.addWidget(self.pressure_min_edits, 3, 3, 1, 1)
        self.overpressure_min_edits = QLineEdit('')
        pane1.addWidget(self.overpressure_min_edits, 4, 3, 1, 1)
        self.afi_min_edits = QLineEdit('')
        pane1.addWidget(self.afi_min_edits, 5, 3, 1, 1)
        self.geo_min_edits = QLineEdit('')
        pane1.addWidget(self.geo_min_edits, 6, 3, 1, 1)
        self.p_a_min_edits = QLineEdit('')
        pane1.addWidget(self.p_a_min_edits, 7, 3, 1, 1)
        self.c_min_edits = QLineEdit('')
        pane1.addWidget(self.c_min_edits, 8, 3, 1, 1)
        self.fd_min_edits = QLineEdit('')
        pane1.addWidget(self.fd_min_edits, 9, 3, 1, 1)

        self.height_max_edits = QLineEdit('')
        pane1.addWidget(self.height_max_edits, 1, 4, 1, 1)
        self.range_max_edits = QLineEdit('')
        pane1.addWidget(self.range_max_edits, 2, 4, 1, 1)
        self.pressure_max_edits = QLineEdit('')
        pane1.addWidget(self.pressure_max_edits, 3, 4, 1, 1)
        self.overpressure_max_edits = QLineEdit('')
        pane1.addWidget(self.overpressure_max_edits, 4, 4, 1, 1)
        self.afi_max_edits = QLineEdit('')
        pane1.addWidget(self.afi_max_edits, 5, 4, 1, 1)
        self.geo_max_edits = QLineEdit('')
        pane1.addWidget(self.geo_max_edits, 6, 4, 1, 1)
        self.p_a_max_edits = QLineEdit('')
        pane1.addWidget(self.p_a_max_edits, 7, 4, 1, 1)
        self.c_max_edits = QLineEdit('')
        pane1.addWidget(self.c_max_edits, 8, 4, 1, 1)
        self.fd_max_edits = QLineEdit('')
        pane1.addWidget(self.fd_max_edits, 9, 4, 1, 1)


        self.yield_button = QPushButton('Calculate Yield')
        pane1.addWidget(self.yield_button, 11, 1, 1, 4)
        self.yield_button.clicked.connect(self.yieldCalc)

        self.integrate_button = QPushButton('Integrate')
        pane1.addWidget(self.integrate_button, 10, 1, 1, 4)
        self.integrate_button.clicked.connect(self.intCalc)

        # Constants - Reed 1972
        self.W_0 = 4.2e12 # Standard reference explosion yield
        self.P_0 = 101325 # Standard pressure
        self.k = 2e-4     
        self.b = 1.19e-4  # Scale height coeff
        self.J_m = 0.375  # Avg positive period of reference explosion
        self.c_m = 347    # Sound speed of reference explosion

        self.blastline_view = pg.GraphicsLayoutWidget()
        self.blastline_canvas = self.blastline_view.addPlot()
        pane2.addWidget(self.blastline_view)
        self.blastline_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)


    def integration_full(self, k, v, b, I, P_a):
        return np.exp(-k*v**2/b/P_a*I)

    def phi(self, f_d, d, W, W_0):
        #scaled distance
        return f_d*d/(W/W_0)**(1/3)
  
    def inverse_gunc(self, p_ans):
        a, b = pso(gunc_error, [8], [12], args=([p_ans, self.J_m, self.W_0, self.P, self.P_0, self.c, self.c_m, self.f_d, self.R, self.P_a, self.k, self.b, self.I, self.cf]), processes=1, swarmsize=1000, maxiter=1000)

        return 10**a[0]

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


        sounding, perturbations = self.bam.atmos.getSounding(lats, lons, elevs, spline=50)

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
                
            f, g, T, P = intscan(S, D, zProfile, wind=True, n_theta=2000, n_phi=2000,\
                    h_tol=330, v_tol=330)
                    # h_tol=self.prefs.pso_min_ang, v_tol=self.prefs.pso_min_dist)

            rf = self.refractiveFactor(S, D, zProfile, D_ANGLE=D_ANGLE)
            trans.append(f)
            ints.append(g)
            ts.append(T)
            ps.append(P)

            rfs.append(rf)


        return trans, ints, ts, ps, rfs

    def intCalc(self):

        stn = self.stn_list[self.current_station]
        if tryFloat(self.height_edits.text()) != None:
            height = tryFloat(self.height_edits.text())
            trans, ints, ts, ps, rfs = self.integrate(height, D_ANGLE=1.5)

            f_val = np.nanmean(trans)
            g_val = np.nanmean(ints)
            t_val = np.nanmean(ts)
            p_val = np.nanmean(ps)
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
            self.c_edits.setText('{:.4f}'.format(t_val))
            self.pressure_edits.setText('{:.4f}'.format(p_val))
            self.p_a_edits.setText('{:.4f}'.format(estPressure(stn.metadata.position.elev)))

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
        w = np.linspace(np.log10(1e8), np.log10(1e12))
        W = 10**w

        if tryFloat(self.height_edits.text()) != None:
            self.R = tryFloat(self.range_edits.text())
            self.P_a = tryFloat(self.p_a_edits.text())
            self.cf = tryFloat(self.geo_edits.text())
            self.I = tryFloat(self.afi_edits.text())
            self.P = tryFloat(self.pressure_edits.text())
            self.c = tryFloat(self.c_edits.text())
            self.f_d = tryFloat(self.fd_edits.text())

            v = 1/2/self.J_m*(self.W_0*self.P/W/self.P_0)**(1/3)*(self.c/self.c_m)
            p_ratio = chem_func(self.phi(self.f_d, self.R, W, self.W_0))*self.P_a*(self.integration_full(self.k, v, self.b, self.I, self.P_0))*self.cf

            self.yieldPlot(p_ratio, W)

        # if tryFloat(self.height_min_edits.text()) != None:
        #     self.R = tryFloat(self.range_min_edits.text())
        #     self.P_a = tryFloat(self.p_a_min_edits.text())
        #     self.cf = tryFloat(self.geo_min_edits.text())
        #     self.I = tryFloat(self.afi_min_edits.text())
        #     self.P = tryFloat(self.pressure_min_edits.text())
        #     self.c = tryFloat(self.c_min_edits.text())
        #     self.f_d = tryFloat(self.fd_min_edits.text())

        #     v = 1/2/self.J_m*(self.W_0*self.P/W/self.P_0)**(1/3)*(self.c/self.c_m)
        #     p_ratio = chem_func(self.phi(self.f_d, self.R, W, self.W_0))*self.P_a*(self.integration_full(self.k, v, self.b, self.I, self.P_0))*self.cf

        #     self.yieldPlot(p_ratio, W, unc='min')

        # if tryFloat(self.height_max_edits.text()) != None:
        #     self.R = tryFloat(self.range_max_edits.text())
        #     self.P_a = tryFloat(self.p_a_max_edits.text())
        #     self.cf = tryFloat(self.geo_max_edits.text())
        #     self.I = tryFloat(self.afi_max_edits.text())
        #     self.P = tryFloat(self.pressure_max_edits.text())
        #     self.c = tryFloat(self.c_max_edits.text())
        #     self.f_d = tryFloat(self.fd_max_edits.text())

        #     v = 1/2/self.J_m*(self.W_0*self.P/W/self.P_0)**(1/3)*(self.c/self.c_m)
        #     p_ratio = chem_func(self.phi(self.f_d, self.R, W, self.W_0))*self.P_a*(self.integration_full(self.k, v, self.b, self.I, self.P_0))*self.cf

        #     self.yieldPlot(p_ratio, W, unc='max')

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
        print('Blast Estimate {:.2E} J'.format(W_final))
        if unc == 'none':
            txt = pg.TextItem("Frag {:}: {:.2E} J".format(self.count, W_final))
            txt.setPos(self.p_rat, np.log10(W_final))
            self.blastline_canvas.addItem(txt)
        self.blastline_canvas.setTitle('Fragmentation Yield Curves for a Given Overpressure')


    def addAngleComp(self, ang1, ang2, deg=False):
        ''' Adds perpendicular angles together to a combined angle
        '''

        if np.isnan(ang1) or np.isnan(ang2):
            return np.nan

        if deg:
            return np.degrees(np.arccos(np.cos(np.radians(ang1))*np.cos(np.radians(ang2))))
        else:
            return np.arccos(np.cos(ang1)*np.cos(ang2))

    def refractiveFactor(self, S, D, zProfile, D_ANGLE=1.5):

        dx, dy, dz = D - S
        tf_ideal_n = np.degrees(np.arctan2(-dz, np.sqrt((dy)**2 + (dx)**2)))
        az_ideal_n = np.degrees(np.arctan2(dx, dy))

        # make sure this angle is high or the rays won't reach the station (like 2000)
        angle = 2000
        _, az_n, tf_n, _ = cyscan(S, D, zProfile, wind=True, n_theta=angle, n_phi=angle,\
                h_tol=330, v_tol=330)
                    # h_tol=self.prefs.pso_min_ang, v_tol=self.prefs.pso_min_dist)

        # if not eigenConsistancy(S, D, az_n, tf_n, zProfile, n_angle=self.prefs.pso_theta, h_tol=self.prefs.pso_min_ang, v_tol=self.prefs.pso_min_dist):
        #     return np.nan

        d_angle_ideal = [np.nan]

        RANGE = 2
        for i in range(RANGE):
            d_az = np.degrees(np.arctan(np.cos(i*2*np.pi/RANGE)*np.tan(np.radians(D_ANGLE))))
            d_tf = np.degrees(np.arctan(np.sin(i*2*np.pi/RANGE)*np.tan(np.radians(D_ANGLE))))

            D_n = anglescan(S, az_n + d_az, tf_n + d_tf, zProfile, wind=True)
                            # h_tol=self.prefs.pso_min_ang, v_tol=self.prefs.pso_min_dist, target=D)

            dx, dy, dz = D_n[0:3] - S

            # if eigenConsistancy(S, D_n[:3], az_n + d_az, tf_n + d_tf, zProfile, h_tol=self.prefs.pso_min_ang, v_tol=self.prefs.pso_min_dist):
        
            tf = np.abs((np.degrees(np.arctan2(-dz, np.sqrt((dy)**2 + (dx)**2)))) - tf_ideal_n)
            az = np.abs((np.degrees(np.arctan2(dx, dy))) - az_ideal_n)

            d_angle_ideal.append(self.addAngleComp(tf, az, deg=True))



        if np.isnan(d_angle_ideal).all():
            return np.nan

        d_angle_ideal = np.nanmean(d_angle_ideal)

        rf = np.sqrt(D_ANGLE/d_angle_ideal)
        
        return rf
