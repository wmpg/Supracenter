
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

I_0 = 1200
mode = "notcneos"
C_list = ["y", "m", "c", "r", "g", "b"]
M_list = [".", "v", "^", "s", "*", "D"]
VEL = 19000

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
    v = VEL

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

def findAreaUnderCurve(h_min, h_max, L=None, avg_curve=None):

    if L is not None:
        h_list, M_list = L.interpCurve(dh=10000)
    else:
        h_list, M_list = avg_curve
        h_list = h_list[::-1]
        M_list = M_list[::-1]

    I_list = 10**((6 - M_list)/2.5)

    # HARD CODED PER EVENT - should change this
    v = VEL

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

        if hasattr(self.bam.setup, "trajectory"):
            self.v = self.bam.setup.trajectory.v
        else:
            self.v = VEL
            print("No trajectory found, using default speed of {:.2f} km/s".format(self.v/1000))
            

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
        
        control_layout = QGridLayout()
        main_layout.addLayout(control_layout, 2, 0, 10, 100)

        lc_control_layout = QGridLayout()
        main_layout.addLayout(lc_control_layout, 2, 101, 10, 100)

        tau_control_layout =QGridLayout()
        main_layout.addLayout(tau_control_layout, 2, 202, 10, 100)

        self.viewbox = QComboBox()
        lc_control_layout.addWidget(self.viewbox, 1, 1, 1, 100)

        self.viewbox.addItem("All")

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
        l, self.tau_edits, b = createFileSearchObj("Tau Curve File", tau_control_layout, 1, width=1, h_shift=0, tool_tip='')
        b.clicked.connect(partial(fileSearch, ['CSV (*.csv)'], self.tau_edits))

        self.redraw_button = createButton("Plot", control_layout, 1, 1, self.redraw, args=[])
        self.cla_button = createButton("Clear All", control_layout, 2, 2, self.clearEnergy, args=[])
        self.add_button = createButton("Add Frag", control_layout, 1, 2, self.addFrag, args=[])
        self.FWHM_button = createButton("Fragmentation Finder", control_layout, 2, 1, self.redraw, args=[True])


        self.prop_table = QScrollArea()
        # self.sources_layout.addWidget(self.sources_table)
        self.prop_table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
        self.prop_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

        self.prop_table.setWidgetResizable(True)

        container = QWidget()
        container.setStyleSheet("""
            QWidget {
                background-color: rgb(0, 0, 0);
                }
            """)
        self.prop_table.setWidget(container)
        self.prop_table_layout = QGridLayout(container)
        self.prop_table_layout.setSpacing(10)

        self.prop_edits = [None]*len(self.light_curve_list)
        
        lc_control_layout.addWidget(self.prop_table, 2, 1, 1, 100)

        for ll, L in enumerate(self.light_curve_list):
            self.viewbox.addItem("{:}".format(L.station))
            l = QLabel(L.station)
            self.prop_edits[ll] = QLineEdit("0.5")

            self.prop_table_layout.addWidget(l, ll+2, 1)
            self.prop_table_layout.addWidget(self.prop_edits[ll], ll+2, 2)

        self.viewbox.currentIndexChanged.connect(partial(self.redraw, True))

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

    def redraw(self, fwhm=False):

        if fwhm:
            print("Fragmentation Finding Mode")
        else:
            print("Plotting Mode")

        self.height_list = []
        self.stat_height_list = []

        # Get Taus

        # Clear Graph
        self.light_curve.ax.clear()
        self.lum_curve.ax.clear()

        for ii in range(len(self.source_widget)):
            # self.tau[ii] = self.source_widget[ii].getTau()

            min_h, max_h = self.source_widget[ii].getHeights()



            try:
                print("Acoustic Energy {:.2E} J".format(self.bam.energy_measurements[ii].chem_pres))

            except:
                print(self.bam.energy_measurements[ii].linear_E)
                continue
            
            try:
                
                # Luminous Energy between heights
                for ll, L in enumerate(self.light_curve_list):

                    self.proportion[ll] = float(self.prop_edits[ll].text())

                    if not (self.viewbox.currentText() == "All" or self.viewbox.currentText() == L.station):
                        continue

                    # Returns dh number of points in list
                    DH = 10000
                    FRAC = 0.05
                    h_list, M_list = L.interpCurve(dh=DH)

                    h_list = np.array(h_list)
                    idx = np.where((h_list > min_h)*(h_list < max_h))
                    dh = np.abs((h_list[1] - h_list[0])*1000)


                    if fwhm:
                        approx_h = (min_h + max_h)/2

                        
                        res = min(enumerate(h_list), key=lambda x: abs(approx_h - x[1]))
                        idx = res[0]

                        M_val = M_list[idx]

                        # Find local maximum around approx_H

                        f = int(DH*FRAC)

                        for m in range(idx - f, idx + f):
                            if M_list[m] < M_val:
                                M_val = M_list[m]
                                idx = m

                        # print("M_VAL = {:}".format(M_val))

                        
                        # print(M_val, np.log10(PROPORTION), M_val - np.log10(PROPORTION))
                        self.light_curve.ax.axhline(y=M_val - np.log10(self.proportion[ll]), color=C_list[ll], \
                            linestyle='dotted')
                        self.light_curve.ax.scatter(h_list[idx], M_val, color=C_list[ll])
                        for ii in range(idx, 0, -1):

                            
                            if M_list[ii] > M_val - np.log10(self.proportion[ll]):
                                # print("M_LEFT = {:}".format(M_list[ii]))
                                M_left = ii
                                break

                        for ii in range(idx, len(M_list)):

                            if M_list[ii] > M_val - np.log10(self.proportion[ll]):
                                # print("M_RIGHT = {:}".format(M_list[ii]))
                                M_right = ii
                                break

                        self.height_list.append([h_list[M_right], h_list[M_left]])
                        self.stat_height_list.append([h_list[M_right], h_list[M_left]])

                    else:
                        self.height_list.append([min_h, max_h])
                        self.stat_height_list.append([min_h, max_h])

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

            self.proportion = [None]*len(self.light_curve_list)
            H1 = 17000/1000
            H2 = 100000/1000
            N = 1600

            h_avg_list = np.linspace(H1, H2, num=N)
            M_avg = []

            for ll, L in enumerate(self.light_curve_list):

                if not (self.viewbox.currentText() == "All" or self.viewbox.currentText() == L.station):
                    continue

                ## Need a way to process light curve with trajectory here
                h, M = L.interpCurve(dh=10000)

                M_list = []
                for h_val in h_avg_list:
                    a = L.getMatH(h_val)

                    M_list.append(a)

                M_avg.append(M_list)

                func = np.poly1d(np.polyfit(h, M, 10))
                xp = np.linspace(h[0], h[-1], 100)
                if self.proportion[ll] is None:
                    self.light_curve.ax.plot(h, M, label="{:}".format(L.station), alpha=1.0, color=C_list[ll])
                else:  
                    self.light_curve.ax.plot(h, M, label="{:} ({:.2f}% of Local Maximum)".format(L.station, self.proportion[ll]*100), alpha=1.0, color=C_list[ll])
                # self.light_curve.ax.plot(xp, func(xp))

            self.light_curve.ax.invert_yaxis()
            # self.light_curve.ax.legend()

            def colAvg(lst):
                #Take the column average of a list of lists
                new_list = []
                for l in lst:
                    new_list.append([10**x for x in l])


                lin_lst = np.nanmean(new_list, axis=0)

                lin_lst = [np.log10(x) for x in lin_lst]
                
                return np.array(lin_lst)

            h_avg = h_avg_list
            M_avg = colAvg(M_avg)

            self.avg_curve = [h_avg, M_avg]

            # self.light_curve.ax.plot(h_avg, M_avg, label="Average Curve")


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
                if not (self.viewbox.currentText() == "All" or self.viewbox.currentText() == L.station):
                    continue
                h_list, M_list = L.interpCurve(dh=10000)
                area = findArea(h/1000, A, h_list, M_list)



            return area

        self.source_widget = [None]*len(self.bam.energy_measurements)

        for ee, energy in enumerate(self.bam.energy_measurements):

            # I don't care which one of these show up they are meaningless unless you set them yourself
            try:
                self.source_widget[ee] = TauEx(energy, self.height_list[ee*len(self.light_curve_list)])
            except:
                self.source_widget[ee] = TauEx(energy, self.height_list[ee])

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

                # energy_area = fragEnergy2Mag(chem_pres_yield, h, tau, self.v)
                # if energy_area is not None and len(energy_area) > 0:
                    
                #     self.light_curve.ax.fill_between(energy_area[:, 1], energy_area[:, 0], color="w", alpha=0.1, label="Fragmentation: {:.1f} km".format(h/1000))

                total_h_min = 9999999
                total_h_max = 0

                print("#######")
                tau_list = []
                for ll, L in enumerate(self.light_curve_list):
                    if not (self.viewbox.currentText() == "All" or self.viewbox.currentText() == L.station):
                        continue
                    try:
                        height_idx = ll + ee*(len(self.bam.energy_measurements) + 1)

                        h_min = self.stat_height_list[height_idx][0]
                        h_max = self.stat_height_list[height_idx][1]
                    except:
                        h_min = float(self.source_widget[ee].min_h_edits.text())
                        h_max = float(self.source_widget[ee].max_h_edits.text())

                    if h_min < total_h_min:
                        total_h_min = h_min

                    if h_max > total_h_max:
                        total_h_max = h_max

                    self.light_curve.ax.axvline(x=h_min, color=C_list[ll], linestyle='--')
                    self.light_curve.ax.axvline(x=h_max, color=C_list[ll], linestyle='--')

                    Height, Magnitude = L.interpCurve(dh=10000)

                    max_idx = np.argmin(abs(Height - h_min))
                    min_idx = np.argmin(abs(Height - h_max))

                    h1 = [Height[min_idx], Height[max_idx]]
                    m1 = [Magnitude[min_idx], Magnitude[max_idx]]

                    E = findAreaUnderCurve(h_min, h_max, L)

                    print("### Station: {:} ###".format(L.station))
                    print("Fragmentation from {:.2f} - {:.2f} km".format(h_min, h_max))
                    print("Energy under curve:       {:.2E}".format(E))
                    print("Expected Acoustic Energy: {:.2E}".format(chem_pres_yield))
                    print("Tau Estimate:       {:.2f} %".format(E/chem_pres_yield*100))
                    print("#########")
                    tau_list.append(E/chem_pres_yield*100)
                    
                    try:
                        station_name = "{:}-{:}".format(energy.station.metadata.network, energy.station.metadata.code)
                    except AttributeError:
                        station_name = "Unknown Station"

                    if ee == 0 or ll == 0:
                        if ee == 0:
                            self.lum_curve.ax.scatter(E/chem_pres_yield*100, h/1000, color=C_list[ll], marker=M_list[ee], \
                                        label="Optical Station: {:}".format(L.station))
                        if ll == 0:
                            self.lum_curve.ax.scatter(E/chem_pres_yield*100, h/1000, color=C_list[ll], marker=M_list[ee], \
                                        label="Infrasound Station: {:}".format(station_name))    
                    else:
                        self.lum_curve.ax.scatter(E/chem_pres_yield*100, h/1000, color=C_list[ll], marker=M_list[ee])

                    try:
                        self.lum_curve.ax.errorbar(E/chem_pres_yield*100, h/1000,\
                            yerr=[[np.abs(h/1000 - h_min)], [np.abs(h/1000 - h_max)]],\
                            fmt=M_list[ee], color=C_list[ll], capsize=5)
                    except:
                        pass
                # # Average
                # Height, Magnitude = self.avg_curve

                # max_idx = np.argmin(abs(Height - h_min))
                # min_idx = np.argmin(abs(Height - h_max))

                # h1 = [Height[min_idx], Height[max_idx]]
                # m1 = [Magnitude[min_idx], Magnitude[max_idx]]

                # E = findAreaUnderCurve(h_min, h_max, avg_curve=self.avg_curve)

                # print("### Station: AVERAGE ###")
                # print("Fragmentation from {:.2f} - {:.2f} km".format(h_min, h_max))
                # print("Energy under curve:       {:.2E}".format(E))
                # print("Expected Acoustic Energy: {:.2E}".format(chem_pres_yield))
                # print("Tau Estimate:       {:.2f} %".format(E/chem_pres_yield*100))
                # print("#########")
                nom_tau = np.mean(tau_list)

                # self.light_curve.ax.fill_between(h1, m1, color="r", alpha=0.3)


                # try:
                #     tau_max = energy.chem_pres*tau/energy.chem_pres_max
                #     tau_min = energy.chem_pres*tau/energy.chem_pres_min 
                # except:
                #     tau_max = tau
                #     tau_min = tau
            
            tau_min = np.nanmin(tau_list)
            tau_max = np.nanmax(tau_list)

            if ee == 0:
                self.lum_curve.ax.scatter(nom_tau, h/1000, color="w", alpha=0.5, marker=M_list[ee], label="Total Per Station")
            else:
                self.lum_curve.ax.scatter(nom_tau, h/1000, color="w", alpha=0.5, marker=M_list[ee])
            try:
                self.lum_curve.ax.errorbar(nom_tau, h/1000, xerr=[[np.abs(nom_tau - tau_min)], [np.abs(nom_tau - tau_max)]],\
                                 yerr=[[np.abs(h/1000 - total_h_min)], [np.abs(h/1000 - total_h_max)]],\
                                 fmt=M_list[ee], capsize=5, color="w", alpha=0.5)
            except TypeError:
                pass
            except UnboundLocalError:
                pass

        self.lum_curve.ax.set_xlabel("Luminous Efficiency [%]")
        self.lum_curve.ax.set_ylabel("Height [km]")
        self.lum_curve.ax.grid(alpha=0.2)
        self.lum_curve.ax.legend()
        self.lum_curve.show()

        self.light_curve.ax.set_xlabel("Height [km]")
        self.light_curve.ax.set_ylabel("Magnitude")
        self.light_curve.ax.grid(alpha=0.2)
        self.light_curve.ax.legend()
        self.light_curve.show()
