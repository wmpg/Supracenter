import numpy as np
import matplotlib.pyplot as plt
import multiprocessing

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.Utils.pso import pso
from supra.GUI.GUITools import createLabelEditObj
from supra.Utils.TryObj import *
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.cyscanIntegration import cyscan as intscan
from supra.Fireballs.SeismicTrajectory import parseWeather
from supra.Utils.Classes import Position
from supra.Supracenter.SPPT import perturb as perturbation_method
from supra.Supracenter.Utils.l137 import estPressure

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

    def __init__(self, setup, stn_list, current_station):

        QWidget.__init__(self)
        self.setWindowTitle('Yields')
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        self.setup = setup
        self.stn_list = stn_list
        self.current_station = current_station

        stylesheet = """ 
        QTabWidget>QWidget>QWidget{background: gray;}
        QLabel{color: white;}
        QCheckBox{color: white;}
        QDockWidget{color: white; background: black;}
        QGroupBox{color: white;}
        QGroupBox{ 
        border: 2px white; 
        border-radius: 0px; }
        QMessageBox{color: white; background: black;} 
        QTableWidget{color: white; background: black;}
        """

        self.setStyleSheet(stylesheet)
        self.count = 0
        layout = QHBoxLayout()
        self.setLayout(layout)

        pane1 = QGridLayout()
        layout.addLayout(pane1)

        pane2 = QVBoxLayout()
        layout.addLayout(pane2)

        self.station_label = QLabel('Station: {:}'.format(self.stn_list[self.current_station].code))
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


        self.yield_button = QPushButton('Calculate Yield')
        pane1.addWidget(self.yield_button, 10, 1, 1, 4)
        self.yield_button.clicked.connect(self.yieldCalc)

        self.integrate_button = QPushButton('Integrate')
        pane1.addWidget(self.integrate_button, 9, 1, 1, 4)
        self.integrate_button.clicked.connect(self.intCalc)

        self.f_d = 0.573
        self.W_0 = 4.2e12
        self.P_0 = 101325
        self.k = 2e-4
        self.b = 1.19e-4
        self.J_m = 0.375
        self.c_m = 347

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
        a, b = pso(gunc_error, [8], [12], args=([p_ans, self.J_m, self.W_0, self.P, self.P_0, self.c, self.c_m, self.f_d, self.R, self.P_a, self.k, self.b, self.I, self.cf]), processes=multiprocessing.cpu_count(), swarmsize=1000, maxiter=1000)

        return 10**a[0]

    def integrate(self, height):
        ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)
        point = self.setup.trajectory.findGeo(height)
        point.pos_loc(self.setup.ref_pos)
        dataset = parseWeather(self.setup)

        ints = []
        ts = []
        ps = []

        for ptb_n in range(self.setup.perturb_times):
            stn = self.stn_list[self.current_station]
            stn.position.pos_loc(self.setup.ref_pos)
            self.sounding = self.perturbGenerate(ptb_n, dataset, self.perturbSetup())
            zProfile, _ = getWeather(np.array([point.lat, point.lon, point.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), self.setup.weather_type, \
                    [ref_pos.lat, ref_pos.lon, ref_pos.elev], self.sounding, convert=False)

            f, g, T, P = intscan(np.array([point.x, point.y, point.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True, \
                                n_theta=1000, n_phi=1000, h_tol=1e-15, v_tol=330)
            ints.append(g)
            ts.append(T)
            ps.append(P)
            # if self.setup.debug:
            #     print(zProfile)

        return ints, ts, ps

    def intCalc(self):
        stn = self.stn_list[self.current_station]
        if tryFloat(self.height_edits.text()) != None:
            height = tryFloat(self.height_edits.text())
            ints, ts, ps = self.integrate(height)

            g_val = np.nanmean(ints)
            t_val = np.nanmean(ts)
            p_val = np.nanmean(ps)

            self.afi_edits.setText('{:.4f}'.format(g_val))
            self.c_edits.setText('{:.4f}'.format(t_val))
            self.pressure_edits.setText('{:.4f}'.format(p_val))
            self.p_a_edits.setText('{:.4f}'.format(estPressure(stn.position.elev)))
            frag_pos = self.setup.trajectory.findGeo(height)
            stn_pos = stn.position
            dist = stn_pos.pos_distance(frag_pos)
            self.range_edits.setText('{:.4f}'.format(dist))

        if tryFloat(self.height_min_edits.text()) != None:
            height = tryFloat(self.height_min_edits.text())

            ints, ts, ps = self.integrate(height)

            g_val = np.nanmean(ints)
            t_val = np.nanmean(ts)
            p_val = np.nanmean(ps)

            self.afi_min_edits.setText('{:.4f}'.format(g_val))
            self.c_min_edits.setText('{:.4f}'.format(t_val))
            self.pressure_min_edits.setText('{:.4f}'.format(p_val))
            self.p_a_min_edits.setText('{:.4f}'.format(estPressure(stn.position.elev)))
            frag_pos = self.setup.trajectory.findGeo(height)
            stn_pos = stn.position
            dist = stn_pos.pos_distance(frag_pos)
            self.range_min_edits.setText('{:.4f}'.format(dist))

        if tryFloat(self.height_max_edits.text()) != None:
            height = tryFloat(self.height_max_edits.text())

            ints, ts, ps = self.integrate(height)

            g_val = np.nanmean(ints)
            t_val = np.nanmean(ts)
            p_val = np.nanmean(ps)

            self.afi_max_edits.setText('{:.4f}'.format(g_val))
            self.c_max_edits.setText('{:.4f}'.format(t_val))
            self.pressure_max_edits.setText('{:.4f}'.format(p_val))
            self.p_a_max_edits.setText('{:.4f}'.format(estPressure(stn.position.elev)))
            frag_pos = self.setup.trajectory.findGeo(height)
            stn_pos = stn.position
            dist = stn_pos.pos_distance(frag_pos)
            self.range_max_edits.setText('{:.4f}'.format(dist))


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

            v = 1/2/self.J_m*(self.W_0*self.P/W/self.P_0)**(1/3)*(self.c/self.c_m)
            p_ratio = chem_func(self.phi(self.f_d, self.R, W, self.W_0))*self.P_a*(self.integration_full(self.k, v, self.b, self.I, self.P_0))*self.cf

            self.yieldPlot(p_ratio, W)

        if tryFloat(self.height_min_edits.text()) != None:
            self.R = tryFloat(self.range_min_edits.text())
            self.P_a = tryFloat(self.p_a_min_edits.text())
            self.cf = tryFloat(self.geo_min_edits.text())
            self.I = tryFloat(self.afi_min_edits.text())
            self.P = tryFloat(self.pressure_min_edits.text())
            self.c = tryFloat(self.c_min_edits.text())

            v = 1/2/self.J_m*(self.W_0*self.P/W/self.P_0)**(1/3)*(self.c/self.c_m)
            p_ratio = chem_func(self.phi(self.f_d, self.R, W, self.W_0))*self.P_a*(self.integration_full(self.k, v, self.b, self.I, self.P_0))*self.cf

            self.yieldPlot(p_ratio, W, unc='min')

        if tryFloat(self.height_max_edits.text()) != None:
            self.R = tryFloat(self.range_max_edits.text())
            self.P_a = tryFloat(self.p_a_max_edits.text())
            self.cf = tryFloat(self.geo_max_edits.text())
            self.I = tryFloat(self.afi_max_edits.text())
            self.P = tryFloat(self.pressure_max_edits.text())
            self.c = tryFloat(self.c_max_edits.text())

            v = 1/2/self.J_m*(self.W_0*self.P/W/self.P_0)**(1/3)*(self.c/self.c_m)
            p_ratio = chem_func(self.phi(self.f_d, self.R, W, self.W_0))*self.P_a*(self.integration_full(self.k, v, self.b, self.I, self.P_0))*self.cf

            self.yieldPlot(p_ratio, W, unc='max')

    def yieldPlot(self, p_ratio, W, unc='none'):
        if unc == 'none':
            self.count += 1
        self.p_rat = tryFloat(self.overpressure_edits.text())
        colour = [(0, 255, 26), (3, 252, 176), (252, 3, 3), (176, 252, 3), (255, 133, 3),
                      (149, 0, 255), (76, 128, 4), (82, 27, 27), (101, 128, 125), (5, 176, 249)]
        ptb_colour = [(0, 255, 26, 150), (3, 252, 176, 150), (252, 3, 3, 150), (176, 252, 3, 150), (255, 133, 3, 150),
                      (149, 0, 255, 150), (76, 128, 4, 150), (82, 27, 27, 150), (101, 128, 125, 150), (5, 176, 249, 150)]

        self.blastline_canvas.setLogMode(False, True)
        if unc == 'none':
            self.nominal = pg.PlotCurveItem(x=p_ratio, y=np.log10(W), pen=colour[self.count-1], name='Fragmentation {:}'.format(self.count))
            self.blastline_canvas.addItem(self.nominal)
        if unc == 'min':
            self.min = pg.PlotCurveItem(x=p_ratio, y=np.log10(W), pen=ptb_colour[self.count-1], name='Fragmentation {:}'.format(self.count))
            self.blastline_canvas.addItem(self.min)
        if unc == 'max':
            self.max = pg.PlotCurveItem(x=p_ratio, y=np.log10(W), pen=ptb_colour[self.count-1], name='Fragmentation {:}'.format(self.count))
            self.blastline_canvas.addItem(self.max)
        try:
            pfill = pg.FillBetweenItem(self.min, self.max, brush=ptb_colour[self.count-1])
            self.blastline_canvas.addItem(pfill)
            self.min = None
            self.max = None
        except:
            pass
            

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


    def perturbSetup(self):

        if self.setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(self.setup, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(self.setup, time=-1)

        else:
            sounding_u = []
            sounding_l = []

        if self.setup.perturb_method == 'ensemble':
            ensemble_file = self.setup.perturbation_spread_file
        else:
            ensemble_file = ''

        if self.setup.perturb_times == 0: self.setup.perturb_times = 1

        if not self.setup.perturb:
            self.setup.perturb_times = 1

        return np.array([sounding_l, sounding_u, ensemble_file])

    def perturbGenerate(self, ptb_n, dataset, perturb_data):

        sounding_l, sounding_u, ensemble_file = perturb_data[0], perturb_data[1], perturb_data[2]

        if ptb_n > 0:
            
            # if self.setup.debug:
            #     print("STATUS: Perturbation {:}".format(ptb_n))

            # generate a perturbed sounding profile
            sounding_p = perturbation_method(self.setup, dataset, self.setup.perturb_method, \
                sounding_u=sounding_u, sounding_l=sounding_l, \
                spread_file=self.setup.perturbation_spread_file, lat=self.setup.lat_centre, lon=self.setup.lon_centre, \
                ensemble_file=ensemble_file, ensemble_no=ptb_n)

        else:
            sounding_p = dataset

        return sounding_p