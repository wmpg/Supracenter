import os
import numpy as np
import pickle
import datetime

import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from functools import partial

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *
from supra.Utils.Classes import *
from supra.GUI.Tools.CustomWidgets import TauEx

from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT
from supra.Utils.Formatting import *

from supra.Lightcurve.light_curve import processLightCurve, readLightCurve

def findArea(height, energy, data_h, data_m, magnitude=None):

    if magnitude is None:

        data_I = 1500*10**(data_m/-2.5)

        #find closest height in data
        h_indx = np.nanargmin(np.abs(data_h - height))


        total_energy = 0
        offset = -1


        while total_energy < energy:

            offset += 1
            if offset == 0:
                total_energy += data_I[h_indx]*np.abs(data_h[h_indx + 1] - data_h[h_indx - 1])/2

            else:
                total_energy += data_I[h_indx + offset]*np.abs(data_h[h_indx + offset + 1] - data_h[h_indx + offset - 1])/2
                total_energy += data_I[h_indx - offset]*np.abs(data_h[h_indx - offset + 1] - data_h[h_indx - offset - 1])/2


        final_data = []
        for i in range(h_indx - offset, h_indx + offset + 1):
            final_data.append([data_m[i], data_h[i]])
        return np.array(final_data)

    else:
        I = 1500*10**(magnitude/-2.5)

        dh = energy/I

        final_data = np.array([[magnitude, height - dh/2],
                               [magnitude, height + dh/2]])
        return final_data


class lumEffDialog(QWidget):

    def __init__(self, bam):

        QWidget.__init__(self)
        
        self.bam = bam
        self.setup = self.bam.setup

        self.v = self.bam.setup.trajectory.v
        self.tau = [5.00]*len(self.bam.energy_measurements)
        self.height_list = []
        for i in range(len(self.bam.energy_measurements)):
            self.height_list.append([None, None])

        self.buildGUI()
        self.processEnergy()
        

    def buildGUI(self):
        self.setWindowTitle('Luminous Efficiency')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        main_layout = QGridLayout()
        self.light_curve = MatplotlibPyQT()
        self.light_curve.ax = self.light_curve.figure.add_subplot(111)
        main_layout.addWidget(self.light_curve, 1, 101, 1, 100)

        self.lum_curve = MatplotlibPyQT()
        self.lum_curve.ax = self.lum_curve.figure.add_subplot(111)
        main_layout.addWidget(self.lum_curve, 1, 202, 1, 100)
        
        self.lightCurve()

        self.sources_table = QScrollArea()
        # self.sources_layout.addWidget(self.sources_table)
        self.sources_table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.sources_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.sources_table.setWidgetResizable(True)

        container = QWidget()
        container.setStyleSheet("""
            QWidget {
                background-color: rgb(0, 0, 0);
                }
            """)
        self.sources_table.setWidget(container)
        self.sources_table_layout = QVBoxLayout(container)
        self.sources_table_layout.setSpacing(10)

        main_layout.addWidget(self.sources_table, 1, 1, 1, 100)
        l, self.tau_edits, b = createFileSearchObj("Tau Curve File", main_layout, 3, width=1, h_shift=0, tool_tip='')
        b.clicked.connect(partial(fileSearch, ['CSV (*.csv)'], self.tau_edits))

        self.redraw_button = createButton("Plot", main_layout, 2, 1, self.redraw, args=[])
        self.cla_button = createButton("Clear All", main_layout, 2, 2, self.clearEnergy, args=[])

        self.setLayout(main_layout)

    def clearEnergy(self):

        self.bam.energy_measurements = []
        self.redraw()

    def plotTaus(self):


        with open(self.tau_edits.text(), 'r+') as f:
            light_curve = []
            for line in f:
                line = line.split(',')
                temp_line = []
        
                for item in line:
                    temp_line.append(float(item.strip()))

        
                light_curve.append(temp_line)

        light_curve = [x for x in light_curve if x != []]

        light_curve = np.array(light_curve)

        h = light_curve[:, 1]
        t = light_curve[:, 0]

        self.lum_curve.ax.plot(t, h)
        self.lum_curve.show()

    def redraw(self):

        self.height_list = []

        # Get Taus
        for ii in range(len(self.source_widget)):
            self.tau[ii] = self.source_widget[ii].getTau()

            min_h, max_h = self.source_widget[ii].getHeights()
            self.height_list.append([min_h, max_h])

        # Remove all Widgets
        for i in reversed(range(self.sources_table_layout.count())): 
            self.sources_table_layout.itemAt(i).widget().setParent(None)

        # Clear Graph
        self.light_curve.ax.clear()
        self.lum_curve.ax.clear()

        # Light Curve
        self.lightCurve()

        # Plot tau
        if self.tau_edits.text() != "":
            self.plotTaus()

        # Data
        self.processEnergy()

    def lightCurve(self):

        if len(self.setup.light_curve_file) > 0 or not hasattr(self.setup, "light_curve_file"):

            light_curve = readLightCurve(self.setup.light_curve_file)

            self.light_curve_list = processLightCurve(light_curve)


            for L in self.light_curve_list:
                h, M = L.interpCurve(dh=10000)

                self.light_curve.ax.plot(h, M, label=L.station)

            self.light_curve.ax.invert_yaxis()
            self.light_curve.ax.legend()


    def processEnergy(self):

        def ballEnergy2Mag(ball_E, v, tau):
            return -2.5*np.log10((ball_E*v*tau)/1500)

        def fragEnergy2Mag(frag_E, h):

            for L in self.light_curve_list:
                M_list, h_list = L.interpCurve(dh=10000)
                area = findArea(h/1000, frag_E, M_list, h_list)

            return area

        self.source_widget = [None]*len(self.bam.energy_measurements)

        for ee, energy in enumerate(self.bam.energy_measurements):

            self.source_widget[ee] = TauEx(energy, self.tau[ee], self.height_list[ee])

            if self.source_widget[ee].mark_for_deletion == True:
                self.bam.energy_measurements[ee] = None
                continue

            self.sources_table_layout.addWidget(self.source_widget[ee])

            tau = self.tau[ee]/100
            h = energy.height

            ### BALLISTIC
            if energy.source_type.lower() == "ballistic":

                lin_e = energy.linear_E
                ws_e = energy.ws_E
                
                lin_mag = ballEnergy2Mag(lin_e, self.v, tau)
                print("Magnitude {:.2f}".format(lin_mag))
                self.light_curve.ax.scatter(h/1000, lin_mag, label="Ballistic Measurement - Linear")

                ws_mag = ballEnergy2Mag(ws_e, self.v, tau)
                print("Magnitude {:.2f}".format(ws_mag))
                self.light_curve.ax.scatter(h/1000, ws_mag, label="Ballistic Measurement - Weak Shock")


            ### FRAGMENTATION
            elif energy.source_type.lower() == "fragmentation":


                chem_pres_yield = energy.chem_pres

                energy_area = fragEnergy2Mag(chem_pres_yield*tau*self.v, h)

                if energy_area is not None:
                    self.light_curve.ax.fill_between(energy_area[:, 1], energy_area[:, 0], color="w", alpha=0.3, label="Fragmentation: {:.1f} km".format(h/1000))

                h_min = float(self.source_widget[ee].min_h_edits.text())
                h_max = float(self.source_widget[ee].max_h_edits.text())

                self.light_curve.ax.axvline(x=h_min, color='w', linestyle='--', alpha=0.3)
                self.light_curve.ax.axvline(x=h_max, color='w', linestyle='--', alpha=0.3)

            self.lum_curve.ax.scatter(tau*100, h/1000, label=energy.source_type)
            self.lum_curve.ax.set_xlabel("Luminous Efficiency [%]")
            self.lum_curve.ax.set_ylabel("Height [km]")
            self.lum_curve.ax.legend()
            self.lum_curve.show()

        self.light_curve.ax.set_xlabel("Height [km]")
        self.light_curve.ax.set_ylabel("Magnitude")
        self.light_curve.ax.legend()
        self.light_curve.show()
