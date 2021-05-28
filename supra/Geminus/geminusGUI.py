import numpy as np
import matplotlib.pyplot as plt

import os

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *

from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT

from supra.Geminus.overpressure2 import overpressureihmod_Ro
from supra.Geminus.geminusSearch import periodSearch, presSearch
from supra.Files.SaveLoad import save
from supra.Lightcurve.light_curve import processLightCurve, readLightCurve
from supra.Utils.Formatting import *

def Efunction(Ro, h):
    k = 1
    p = 10*101.325*np.exp(-0.00012*h)
    return Ro**2*p*k

class Geminus(QWidget):

    def __init__(self, bam, prefs):

        QWidget.__init__(self)
       
        self.bam = bam
        self.prefs = prefs

        if not hasattr(bam, "infra_curve"):
            self.bam.infra_curve = []

        self.buildGUI()

        self.current_Ro, self.current_height = None, None


    def buildGUI(self):
        self.setWindowTitle('Geminus (Silber 2014)')

        app_icon = QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)
        
        layout = QHBoxLayout()
        self.setLayout(layout)

        graph_layout = QVBoxLayout()
        layout.addLayout(graph_layout)

        right_panels = QVBoxLayout()
        layout.addLayout(right_panels)

        input_panels = QGridLayout()
        right_panels.addLayout(input_panels)

        output_panels = QGridLayout()
        right_panels.addLayout(output_panels)

        stn_name_list = []
        for stn in self.bam.stn_list:
            stn_name_list.append("{:}-{:}".format(stn.metadata.network, stn.metadata.code))

        control_panel = QGridLayout()
        input_panels.addLayout(control_panel, 6, 1, 3, 3)

        _, self.source_height = createLabelEditObj('Source Height Along Trajectory [m]', input_panels, 1, width=1, h_shift=0, tool_tip='')
        _, self.station_combo = createComboBoxObj('Station', input_panels, 2, items=stn_name_list, width=1, h_shift=0, tool_tip='')
     
        _, self.blast_radius = createLabelEditObj('Blast Radius [m]', input_panels, 3, width=1, h_shift=0, tool_tip='')
        _, self.dom_period = createLabelEditObj('Dominant Period [s]', input_panels, 4, width=1, h_shift=0, tool_tip='')
        _, self.over_pres = createLabelEditObj('Overpressure [Pa]', input_panels, 5, width=1, h_shift=0, tool_tip='')  
        self.vary_period = createToggle('Vary Period', control_panel, 1, width=1, h_shift=1, tool_tip='')
        self.add_winds = createToggle('Include Winds', control_panel, 2, width=1, h_shift=1, tool_tip='')
        self.doppler = createToggle('Doppler Shift', control_panel, 3, width=1, h_shift=1, tool_tip='')         
        self.overpressure_run = createButton("Run Blast Radius Simulation", control_panel, 4, 1, self.overpressure, args=["normal"])
        self.overpressure_period_finder = createButton("Run Period Search", control_panel, 4, 2, self.overpressure, args=["period"])
        self.overpressure_pres_finder = createButton("Run Overpressure Search", control_panel, 4, 3, self.overpressure, args=["pres"])
        self.pro_sim = createButton("Period vs. Blast Radius", control_panel, 5, 1, self.overpressure, args=["pro"])
        self.proE_sim = createButton("Period vs. Energy", control_panel, 5, 2, self.overpressure, args=["proE"])
        self.infra_curve = createButton("Infrasound Curve", control_panel, 5, 3, self.infraCurve)
        self.clear_infra = createButton("Clear Curve", control_panel, 6, 1, self.clearInfra)

        self.vary_period.setChecked(True)
        self.add_winds.setChecked(True)
        self.doppler.setChecked(True)

        self.overpressure_plot = MatplotlibPyQT()
        self.overpressure_plot.ax = self.overpressure_plot.figure.add_subplot(111)
        graph_layout.addWidget(self.overpressure_plot)

        theme(self)

    def clearInfra(self):
        self.bam.infra_curve = []

        print(printMessage("status"), "Cleared infrasound curve data!")

    def infraCurve(self):
        
        if not hasattr(self, "current_lin_Ro"):
            if self.prefs.debug:
                print(printMessage("warning"), " Run a search for overpressure or period first before plotting a point!")
            errorMessage("No points to plot!", 1, detail='Please use one of the searches to find Ro for both weak-shock and linear!')
            return None   

        if not self.current_lin_Ro is None and self.current_ws_Ro and not self.current_height is None:
            self.bam.infra_curve.append([self.current_lin_Ro, self.current_ws_Ro, self.current_height])
        
        if len(self.bam.infra_curve) == 0:
            errorMessage("No points to plot!", 1, detail='Please use one of the searches to find Ro!')    
            return None

        ax1 = plt.subplot(2, 1, 1)
        E_lin = []
        E_ws = []

        for point in self.bam.infra_curve:

            h = point[2]
            E_lin.append(Efunction(point[0], h))
            E_ws.append(Efunction(point[1], h))

        ax1.scatter(h/1000, E_lin, label='Linear')
        ax1.scatter(h/1000, E_ws, label="Weak Shock")

        ax1.set_xlabel("Height [km]")
        ax1.set_ylabel("Energy per Unit Length [J/m]")

        ax2 = plt.subplot(2, 1, 2, sharex=ax1)


        light_curve = readLightCurve(self.bam.setup.light_curve_file)

        light_curve_list = processLightCurve(light_curve)

        for L in light_curve_list:
            ax2.scatter(L.h, L.M, label=L.station)
            # light_curve_curve = pg.ScatterPlotItem(x=L.M, y=L.t)
            # self.light_curve_canvas.addItem(light_curve_curve)

        ax2.set_xlabel("Height [km]")
        ax2.set_ylabel("Absolute Magnitude")
        plt.gca().invert_yaxis()
        plt.legend()

        # ax3 = plt.subplot(3, 1, 3, sharex=ax2)
        # v = self.bam.setup.trajectory.v
        # ax3.scatter(np.array(h)/1000, np.array(E_lin)/v, label='Linear')
        # ax3.scatter(np.array(h)/1000, np.array(E_ws)/v, label="Weak Shock")
        # for L in light_curve_list:
        #     ax3.scatter(L.h, 10**(-0.4*np.array(L.M)), label=L.station)

        # ax3.set_xlabel("Height [km]")
        # ax3.set_ylabel("?? Max Intensity ??")

        # plt.legend()
        plt.show()




    def overpressure(self, mode):

        wind = self.add_winds.isChecked()
        dopplershift = self.doppler.isChecked()

        if self.prefs.debug:
            print(printMessage("debug"), " Running Geminus on mode '{:}'".format(mode))

        self.overpressure_plot.clear()

        traj = self.bam.setup.trajectory

        try:
            source = traj.findGeo(float(self.source_height.text()))
        except ValueError as e:
            if self.prefs.debug:
                print(printMessage("Error"), " No source height given!")
            errorMessage("Cannot read source height!", 2, detail='{:}'.format(e))

            return None


        source_list = [source.lat, source.lon, source.elev/1000]

        stat_idx = self.station_combo.currentIndex()
        stat = self.bam.stn_list[stat_idx]
        stat_pos = stat.metadata.position
        stat = [stat_pos.lat, stat_pos.lon, stat_pos.elev/1000]

        try:
            Ro = float(self.blast_radius.text())
        except:
            Ro = None

        v = traj.v/1000

        theta = 90 - traj.zenith.deg

        dphi = np.degrees(np.arctan2(stat_pos.lon - source.lon, stat_pos.lat - source.lat)) - traj.azimuth.deg

        # Switch 3 doesn't do anything in this version of overpressure.py
        sw = [self.vary_period.isChecked(), 0, 1]
        
        lat =   [source.lat, stat_pos.lat]
        lon =   [source.lon, stat_pos.lon]
        elev =  [source.elev, stat_pos.elev]

        sounding, _ = self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev)

        pres = 10*101.325*np.exp(-0.00012*sounding[:, 0])
        
        sounding_pres = np.zeros((sounding.shape[0], sounding.shape[1]+1))
        sounding_pres[:, :-1] = sounding
        sounding_pres[:, -1] = pres
        sounding_pres[:, 1] -= 273.15

        sounding_pres = np.flip(sounding_pres, axis=0)
        self.sounding_pres = sounding_pres

        gem_inputs = [source_list, stat, v, theta, dphi, sounding_pres, sw, wind, dopplershift]

        if mode == "normal":
            try:
                tau, tauws, Z, sR, inc, talt, dpws, dp, it = \
                            overpressureihmod_Ro(source_list, stat, Ro, v, theta, dphi, sounding_pres, sw, wind=wind, dopplershift=dopplershift)
            except TypeError as e:
                errorMessage("Error in running Geminus!", 2, detail='{:}'.format(e))
                return None

            self.overpressure_plot.ax.plot(tau[0:it], Z[0:it], 'r-', label="Weak Shock Period Change")
            self.overpressure_plot.ax.plot(tau[it-1:], Z[it-1:], 'b-', label="Stable Period")
            self.overpressure_plot.ax.plot(tauws[it-1:], Z[it-1:], 'm-', label="Weak Shock: No Transition")

            self.overpressure_plot.ax.scatter([tau[it-1]], [Z[it-1]])
            
            self.overpressure_plot.ax.set_xlabel("Signal Period [s]")
            self.overpressure_plot.ax.set_ylabel("Geopotential Height [km]")
            self.overpressure_plot.ax.legend()
            self.overpressure_plot.show()


            print('Geminus Output')
            print('=========================================================')
            print('Period (weak shock):     {:3.4f} s'.format(tauws[-1]))
            print('  Frequency (weak shock):   {:3.3f} Hz'.format(1/tauws[-1]))
            print('Period (linear):         {:3.4f} s'.format(tau[-1]))
            print('  Frequency (linear):       {:3.3f} Hz'.format(1/tau[-1]))
            print('Slant range:             {:5.2f} km'.format(sR))
            print('Arrival (inclination):   {:3.4f} deg'.format(np.degrees(inc)%360))
            print('Transition height:       {:3.3f} km'.format(talt))
            print('Overpressure (weak shock):     {:3.4f} Pa'.format(dpws[-1]))
            print('Overpressure (linear):         {:3.4f} Pa'.format(dp[-1]))

        elif mode == "period":

            p = float(self.dom_period.text())

            Ro_ws, Ro_lin, weak_path, lin_path, tau, Z, it = periodSearch(p, gem_inputs, paths=True)

            self.overpressure_plot.ax.plot(weak_path, Z, 'r-', label="Weak Shock")
            self.overpressure_plot.ax.plot(lin_path, Z, 'b-', label="Linear")

            self.overpressure_plot.ax.scatter([tau[it-1]], [Z[it-1]])
            
            self.overpressure_plot.ax.set_xlabel("Signal Period [s]")
            self.overpressure_plot.ax.set_ylabel("Geopotential Height [km]")
            self.overpressure_plot.ax.legend()
            self.overpressure_plot.show()


            print('Geminus Output')
            print('=========================================================')
            print("Blast Radius (Weak-Shock): {:.2f} m".format(Ro_ws))
            print("Blast Radius (Linear): {:.2f} m".format(Ro_lin))

        elif mode == "pres":
            
            p = float(self.over_pres.text())

            Ro_ws, Ro_lin, weak_path, lin_path, tau, Z, it = presSearch(p, gem_inputs, paths=True)

            self.overpressure_plot.ax.plot(weak_path, Z, 'r-', label="Weak Shock")
            self.overpressure_plot.ax.plot(lin_path, Z, 'b-', label="Linear")

            self.overpressure_plot.ax.scatter([tau[it-1]], [Z[it-1]])
            
            self.overpressure_plot.ax.set_xlabel("Signal Period [s]")
            self.overpressure_plot.ax.set_ylabel("Geopotential Height [km]")
            self.overpressure_plot.ax.legend()
            self.overpressure_plot.show()


            print('Geminus Output')
            print('=========================================================')
            print("Blast Radius (Weak-Shock): {:.2f} m".format(Ro_ws))
            print("Blast Radius (Linear): {:.2f} m".format(Ro_lin))

        elif mode == "pro" or mode == "proE":

            Ro = np.linspace(0.01, 30.00, 30)

            tau_list = []
            tau_ws_list = []

            for R in Ro:

                try:
                    tau, tauws, Z, sR, inc, talt, dpws, dp, it = \
                            overpressureihmod_Ro(source_list, stat, R, v, theta, dphi, sounding_pres, sw, wind=wind, dopplershift=dopplershift)
                except TypeError as e:
                    errorMessage("Error in running Geminus!", 2, detail='{:}'.format(e))
                    return None

                tau_list.append(tau[-1])
                tau_ws_list.append(tauws[-1])



            if mode == "pro":
                self.overpressure_plot.ax.plot(Ro, np.array(tau_list), 'b-', label="Linear")
                self.overpressure_plot.ax.plot(Ro, np.array(tau_ws_list), 'r-', label="Weak Shock")
                
                self.overpressure_plot.ax.set_xlabel("Blast Radius [m]")
                self.overpressure_plot.ax.set_ylabel("Period [s]")

            elif mode == "proE":

                self.overpressure_plot.ax.plot(Efunction(Ro, source.elev), np.array(tau_list), 'b-', label="Linear")
                self.overpressure_plot.ax.plot(Efunction(Ro, source.elev), np.array(tau_ws_list), 'r-', label="Weak Shock")

                
                self.overpressure_plot.ax.set_xlabel("Energy per Unit Length [J/m]")
                self.overpressure_plot.ax.set_ylabel("Period [s]")

            self.overpressure_plot.ax.legend()
            self.overpressure_plot.show()


        self.current_height = float(self.source_height.text())
