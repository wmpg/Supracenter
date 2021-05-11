import numpy as np

import os

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *

from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT

from supra.Geminus.overpressure import overpressureihmod_Ro

class Geminus(QWidget):

    def __init__(self, bam, prefs):

        QWidget.__init__(self)
       
        self.bam = bam

        self.buildGUI()


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


        _, self.source_height = createLabelEditObj('Source Height Along Trajectory [m]', input_panels, 1, width=1, h_shift=0, tool_tip='')
        _, self.station_combo = createComboBoxObj('Station', input_panels, 2, items=stn_name_list, width=1, h_shift=0, tool_tip='')
        self.vary_period = createToggle('Vary Period', input_panels, 3, width=1, h_shift=1, tool_tip='')
        _, self.blast_radius = createLabelEditObj('Blast Radius [m]', input_panels, 4, width=1, h_shift=0, tool_tip='')
        self.overpressure_run = createButton("Run", input_panels, 5, 2, self.overpressure, args=[])

        self.overpressure_plot = MatplotlibPyQT()
        graph_layout.addWidget(self.overpressure_plot)

        self.weak_shock_period = createLabel("Period (Weak-Shock):              s", output_panels, 1)
        self.weak_shock_freq = createLabel("Frequency (Weak-Shock):           Hz", output_panels, 2)
        self.linear_period = createLabel("Period (Linear):                  s", output_panels, 3)
        self.linear_freq = createLabel("Frequency (Linear):               Hz", output_panels, 4)
        self.slant_range = createLabel("Slant Range:                      km", output_panels, 5)
        self.arrival_inclination = createLabel("Arrival Inclination:              deg", output_panels, 6)
        self.tansition_height = createLabel("Transition Height:                km", output_panels, 7)
        self.weak_shock_pres = createLabel("Overpressure (Weak-Shock):        Pa", output_panels, 8)
        self.linear_pres = createLabel("Overpressure (Linear):            Pa", output_panels, 9)




        theme(self)

    def overpressure(self):

        traj = self.bam.setup.trajectory

        source = traj.findGeo(float(self.source_height.text()))

        source_list = [source.lat, source.lon, source.elev/1000]

        stat_idx = self.station_combo.currentIndex()
        stat = self.bam.stn_list[stat_idx]
        stat_pos = stat.metadata.position
        stat = [stat_pos.lat, stat_pos.lon, stat_pos.elev/1000]

        Ro = float(self.blast_radius.text())

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

        tau, tauws, Z, sR, inc, talt, dpws, dp, it = \
                    overpressureihmod_Ro(source_list, stat, Ro, v, theta, dphi, sounding_pres, sw)

        self.overpressure_plot.ax.plot(tau[0:it], Z[0:it], 'r-', label="Weak Shock Period Change")
        self.overpressure_plot.ax.plot(tau[it-1:], Z[it-1:], 'b-', label="Stable Period")
        self.overpressure_plot.ax.plot(tauws[it-1:], Z[it-1:], 'm-', label="Weak Shock: No Transition")

        self.overpressure_plot.ax.scatter([tau[it-1]], [Z[it-1]])
        
        self.overpressure_plot.ax.set_xlabel("Signal Period [s]")
        self.overpressure_plot.ax.set_ylabel("Geopotential Height [km]")
        self.overpressure_plot.ax.legend()
        self.overpressure_plot.show()

        self.weak_shock_period.setText("Period (Weak-Shock):{:14.4f} s".format(tauws[-1]))
        self.weak_shock_freq.setText("Frequency (Weak-Shock):{:11.3f} Hz".format(1/tauws[-1]))
        self.linear_period.setText("Period (Linear):{:18.4f} s".format(tau[-1]))
        self.linear_freq.setText("Frequency (Linear):{:15.4f} Hz".format(1/tau[-1]))
        self.slant_range.setText("Slant Range:{:22.2f} km".format(sR))
        self.arrival_inclination.setText("Arrival Inclination:{:14.4f} deg".format(inc))
        self.tansition_height.setText("Transition Height:{:17.3f} km".format(talt))
        self.weak_shock_pres.setText("Overpressure (Weak-Shock):{:8.4f} Pa".format(dpws[-1]))
        self.linear_pres.setText("Overpressure (Linear):{:12.4f} Pa".format(dp[-1]))

        print('FINAL OUTPUT FOR THE WEAK SHOCK')
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