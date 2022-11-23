
###############
### IMPORTS
###############
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
from supra.Utils.Energy import EnergyObj

from supra.Lightcurve.light_curve import *

from supra.Files.SaveLoad import save

I_0 = 3030
mode = "cneos"

def findArea(height, energy, data_h, data_m, magnitude=None):
    """ A fragmentation is shown on the light curve as the 
    area under the intensity-height plot with the height 
    boundaries of the fragmentation.

    Given the height of largest amplitude, this function will
    symmetrically extend the heights until the total area under
    the light curve is just larger than the energy given.

    Inputs:
    height [float] - height of the fragmentation in km
    energy [float] - energy of the fragmentation in J
    data_h [list] - list of height values in the light curve in km
    data_m [list] - list of magnitude values in the light curve
    magnitude [float] - use if you want the area to be under a constant
        magnitude instead of the data_m light curve values. Set to None to
        ignore

    Outputs:
    final_data [ndarray] - array of points [magnitude, height] following the
    light curve but limited to the area calculated
    """

    # HARD CODED PER EVENT - should change this
    v = 16400

    # Use light curve
    if magnitude is None:

        # From CNEOS light curves
        if mode == "cneos":
            data_I = 10**((6 - data_m)/2.5)
        
        # From cameras using an I_0, ie/ Borovicka uses 1500 W
        else:
            data_I = I_0*10**((data_m)/-2.5)

        #find closest height in data
        h_indx = np.nanargmin(np.abs(data_h - height))

        total_energy = 0
        offset = -1

        # Main loop
        while total_energy < energy:

            offset += 1

            if offset == 0:

                dI = data_I[h_indx]
                try:
                    dh = np.abs(data_h[h_indx + 1] - data_h[h_indx - 1])/2*1000
                except IndexError:

                    # Catch if area exceeds light curve data
                    print("Index Error: Energy needs to be {:.2E} J, but curve only contains {:.2E} J".format(energy, total_energy))
                    return []

                total_energy += dI*dh/v*4*np.pi

            else:
                try:
                    # Search both sides symmetrically - not sure if this is reasonable, but works for now
                    total_energy += data_I[h_indx + offset]*np.abs(data_h[h_indx + offset + 1] - data_h[h_indx + offset - 1])*1000/2/v*4*np.pi
                    total_energy += data_I[h_indx - offset]*np.abs(data_h[h_indx - offset + 1] - data_h[h_indx - offset - 1])*1000/2/v*4*np.pi
                except IndexError:

                    # Catch if area exceeds light curve data
                    print("Index Error: Energy needs to be {:.2E} J, but curve only contains {:.2E} J".format(energy, total_energy))
                    return []

        final_data = []
        for i in range(h_indx - offset, h_indx + offset + 1):
            final_data.append([data_m[i], data_h[i]])
        return np.array(final_data)

    # Use constant magnitude
    else:
        if mode == "cneos":
            I = 10**((6 - magnitude)/-2.5)
        else:
            I = I_0*10**((magnitude)/-2.5)

        dh = energy/I

        final_data = np.array([[magnitude, height - dh/2],
                               [magnitude, height + dh/2]])
        return final_data

def findAreaUnderCurve(h_min, h_max, L):

    h_list, M_list = L.interpCurve(dh=10000)
    I_list = 10**((6 - M_list)/2.5)

    # HARD CODED PER EVENT - should change this
    v = 16400


    #find closest height in data
    min_indx = np.nanargmin(np.abs(h_list - h_min))
    max_indx = np.nanargmin(np.abs(h_list - h_max))

    h_list = h_list[max_indx:min_indx]
    I_list = I_list[max_indx:min_indx]

    try:
        dh = np.abs(h_list[1] - h_list[0])
    except:
        return 0

    total_energy = 0

    for ii in range(len(I_list) - 1):
        total_energy += (I_list[ii] + I_list[ii+1])*dh*1000/2/v*4*np.pi

    return total_energy

class lumEffDialog(QWidget):

    def __init__(self, bam):

        QWidget.__init__(self)
        
        self.bam = bam
        self.setup = self.bam.setup

        self.v = self.bam.setup.trajectory.v

        # Assume they are all 5% now, and change later
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
        self.add_button = createButton("Add Frag", main_layout, 2, 3, self.addFrag, args=[])

        self.setLayout(main_layout)

    def clearEnergy(self):

        # TODO - whenever the number of fragmentations is changed the program needs to be restarted
        self.bam.energy_measurements = []
        self.redraw()

    def addFrag(self):

        ### Currently a terminal input, will put to GUI later...
        a = EnergyObj()
        a.source_type = "Fragmentation"
        a.height = float(input("Height of Fragmentation in m: "))
        energy = float(input("Energy of Fragmentation in kT: "))
        energy_max = float(input("Max Energy [kT]: "))
        energy_min = float(input("Min Energy [kT]: "))

        a.chem_pres = energy*4.184e12
        a.chem_pres_max = energy_max*4.184e12
        a.chem_pres_min = energy_min*4.184e12

        self.bam.energy_measurements.append(a)

        save(self.bam, file_check=False)

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

            try:
                print("Acoustic Energy {:.2E} J".format(self.bam.energy_measurements[ii].chem_pres))

            except:
                print(self.bam.energy_measurements[ii].linear_E)
                continue
            
            try:
                # Luminous Energy between heights
                for L in self.light_curve_list:
                    h_list, M_list = L.interpCurve(dh=10000)

                h_list = np.array(h_list)
                idx = np.where((h_list > min_h)*(h_list < max_h))
                dh = np.abs((h_list[1] - h_list[0])*1000)

                mags = M_list[idx]
                intensity = 0
                for m in mags:
                    
                    dI = I_0*10**((m)/-2.5)
                    #dI = 10**(m - 6)/(-2.5)
                    intensity += dI*dh
                v = self.v
            except TypeError as e:
                print(e)

            #print("Luminous Energy {:.2E} J".format(intensity))
            #print("Tau {:.2f} %".format(intensity/self.bam.energy_measurements[ii].chem_pres/v*100))


            

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


            try:
                # OTHER
                light_curve = readLightCurve(self.setup.light_curve_file)
                self.light_curve_list = processLightCurve(light_curve)

            except:
            
                # CNEOS USG DATA
                self.light_curve_list = readCNEOSlc(self.setup.light_curve_file)
                self.light_curve_list[0].estimateHeight(self.bam.setup.trajectory)

            if light_curve is None:
                # CNEOS USG DATA
                self.light_curve_list = readCNEOSlc(self.setup.light_curve_file)
                self.light_curve_list[0].estimateHeight(self.bam.setup.trajectory)

            for L in self.light_curve_list:

                ## Need a way to process light curve with trajectory here
                h, M = L.interpCurve(dh=10000)

                func = np.poly1d(np.polyfit(h, M, 10))
                xp = np.linspace(h[0], h[-1], 100)
                self.light_curve.ax.plot(h, M)#, label=L.station)
                self.light_curve.ax.plot(xp, func(xp))

            self.light_curve.ax.invert_yaxis()
            # self.light_curve.ax.legend()


    def processEnergy(self):

        def ballEnergy2Mag(ball_E, v, tau):
            return -2.5*np.log10((ball_E*v*tau)/1500)

        def fragEnergy2Mag(frag_E, h, tau, v):

            """
            Total acoustic energy: frag_E
            Total luminous energy: Area/v

            given a test value of tau, what is Area?

            tau = (Area/v) / (frag_E)
            Area = frag_E*v*tau

            Display this area under curve at fragmentation height for visual comparison
            """
            print("Acoustic Energy: {:.2E} J".format(frag_E))
            print("Velocity: {:.2} m/s".format(v))
            print("Tau: {:.2f} %".format(tau*100))

            A = frag_E*tau

            print("Needed Area: {:.2E} J".format(A))

            for L in self.light_curve_list:
                h_list, M_list = L.interpCurve(dh=10000)
                area = findArea(h/1000, A, h_list, M_list)



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

                tau_max = tau
                tau_min = tau
            ### FRAGMENTATION
            elif energy.source_type.lower() == "fragmentation":


                chem_pres_yield = energy.chem_pres

                energy_area = fragEnergy2Mag(chem_pres_yield, h, tau, self.v)
                if energy_area is not None and len(energy_area) > 0:
                    
                    self.light_curve.ax.fill_between(energy_area[:, 1], energy_area[:, 0], color="w", alpha=0.1, label="Fragmentation: {:.1f} km".format(h/1000))

                h_min = float(self.source_widget[ee].min_h_edits.text())
                h_max = float(self.source_widget[ee].max_h_edits.text())

                self.light_curve.ax.axvline(x=h_min, color='w', linestyle='--')
                self.light_curve.ax.axvline(x=h_max, color='w', linestyle='--')

                for L in self.light_curve_list:
                    Height, Magnitude = L.interpCurve(dh=10000)

                    max_idx = np.argmin(abs(Height - h_min))
                    min_idx = np.argmin(abs(Height - h_max))

                    h1 = [Height[min_idx], Height[max_idx]]
                    m1 = [Magnitude[min_idx], Magnitude[max_idx]]

                    E = findAreaUnderCurve(h_min, h_max, L)

                print("#########")
                print("Fragmentation from {:.2f} - {:.2f} km".format(h_min, h_max))
                print("Energy under curve:       {:.2E}".format(E))
                print("Expected Acoustic Energy: {:.2E}".format(chem_pres_yield))
                print("Tau Estimate:       {:.2f} %".format(E/chem_pres_yield*100))
                print("#########")

                self.light_curve.ax.fill_between(h1, m1, color="r", alpha=0.3)


                try:
                    tau_max = energy.chem_pres*tau/energy.chem_pres_max
                    tau_min = energy.chem_pres*tau/energy.chem_pres_min 
                except:
                    tau_max = tau
                    tau_min = tau
            
            self.lum_curve.ax.scatter(tau*100, h/1000, label=energy.source_type)
            try:
                self.lum_curve.ax.errorbar(tau*100, h/1000, xerr=[[np.abs(tau*100 - tau_min*100)], [np.abs(tau*100 - tau_max*100)]],\
                                 yerr=[[np.abs(h/1000 - energy_area[0, 1])], [np.abs(h/1000 - energy_area[-1, 1])]],\
                                 fmt="o", capsize=5)
            except TypeError:
                pass
            except UnboundLocalError:
                pass

        self.lum_curve.ax.set_xlabel("Luminous Efficiency [%]")
        self.lum_curve.ax.set_ylabel("Height [km]")
        # self.lum_curve.ax.legend()
        self.lum_curve.show()

        self.light_curve.ax.set_xlabel("Height [km]")
        self.light_curve.ax.set_ylabel("Magnitude")
        # self.light_curve.ax.legend()
        self.light_curve.show()
