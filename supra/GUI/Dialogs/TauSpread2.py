import os
import numpy as np
import pickle

from datetime import datetime, timedelta
import matplotlib.pyplot as plt


from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *
from supra.Utils.Classes import *
from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT
from supra.GUI.Dialogs.Yields import *

from supra.GUI.Tabs.SupracenterSearch import getStationData

from mpl_toolkits.basemap import Basemap

from supra.Supracenter.cyscan5 import cyscan

from supra.Supracenter.anglescanYieldInt import anglescan as intscan

from supra.Yields.YieldCalcs import *
from supra.Yields.YieldFuncs import yieldFromScaledDistance

from supra.Utils.pso import pso

import multiprocessing

FILENAME = "F:\\Desktop\\Nord_09_27_22.csv"
#INPUT_FILE = "F:\\Documents\\Meteor_Research\\Event\\Romania\\Romania_Yield_Search_Picks.csv"
#INPUT_FILE = "F:\\Documents\\Meteor_Research\\Event\\alaska_picks_2.csv"
INPUT_FILE = "F:\\Documents\\Meteor_Research\\Event\\meteoroids2020_picks_3.csv"

search_min = [59.0, 16.0, 42000]
search_max = [60.0, 17.5, 50000]


class Pick():

    def __init__(self, stat):

        self.group = int(stat[0])
        self.network = stat[1]
        self.code = stat[2]
        self.lat = float(stat[3])
        self.lon = float(stat[4])
        self.elev = float(stat[5])
        self.time = float(stat[7])
        self.delp = float(stat[9])
        self.period = float(stat[10])
        self.infra = True

        if self.delp < 0:
            self.delp = np.nan
            self.infra = False

        if self.period < 0:
            self.period = np.nan
            self.infra = False

        self.position = Position(self.lat, self.lon, self.elev)

    def __str__(self):

        if self.infra:
            A = "Infrasound Station: "
        else:
            A = "Seismic Station: "

        A += "{:}-{:}".format(self.network, self.code)

        return A

def multiLoop(x, *multi_args):


    pick_list, ref_pos, FILENAME, bam = multi_args

    source = Position(x[0], x[1], x[2])

    total_resid = 0
    total_stats = 0
    yield_list = []

    for ss, stat in enumerate(pick_list):

        print("Working...")



        stat_time = stat.time
        stat_delp = stat.delp
        stat_period = stat.period


        stat_pos = stat.position

        latrange =   [source.lat, stat_pos.lat]
        lonrange =   [source.lon, stat_pos.lon]
        elevrange =  [source.elev, stat_pos.elev]

        sounding, perturbations = bam.atmos.getSounding(lat=latrange, lon=lonrange, heights=elevrange, spline=100, ref_time=bam.setup.fireball_datetime)

        P_a = sounding[0, 4]
        P = sounding[-1, 4]

        stat_pos.pos_loc(source)
        source.pos_loc(source)

        R = source.pos_distance(stat_pos)

        source.z = source.elev
        stat_pos.z = stat_pos.elev

        S = np.array([source.x, source.y, source.z])
        D = np.array([stat_pos.x, stat_pos.y, stat_pos.z])

        r, tr, f_particle = cyscan(S, D, \
            sounding, trace=True, plot=False, particle_output=True, debug=False, \
            wind=True, print_times=True, processes=1)

        T, az, tf, err = r

        if np.isnan(T):


            continue


        T = T - stat_time

        total_resid += np.abs(T)
        total_stats += 1

        v = 1/stat_period

        if np.isnan(stat_delp):
            continue

        Yield_chem = overpressure2Yield(stat_delp, R, x[2])

        print("Solution to: {:.2f} N {:.2f} {:.2f} km - E = {:.2f} s {:.2E} J".format(x[0], x[1], x[2]/1000, T, Yield_chem))

        yield_list.append(Yield_chem)

    if total_stats == 0:
        return np.inf

    total_resid = total_resid/total_stats

    with open(FILENAME, "a+") as f:

        for y in yield_list:
            f.write("{:}, {:}, {:}, {:}, {:}, {:}, {:}\n".format(x[0], x[1], x[2], total_resid, y, total_stats, len(pick_list)))


    return total_resid

class TauSpreadGUI(QWidget):

    def __init__(self, bam, prefs):

        QWidget.__init__(self)
        
        self.bam = bam
        self.prefs = prefs

        self.W_0 = 4.184e6

        self.buildGUI()
        self.initGUI()

    def buildGUI(self):

        self.setWindowTitle('Yield Spread')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)

        # Plots
        self.timing_plot = MatplotlibPyQT()
        self.timing_plot.ax1 = self.timing_plot.figure.add_subplot(121)
        self.timing_plot.ax2 = self.timing_plot.figure.add_subplot(122)
        layout.addWidget(self.timing_plot, 1, 1, 1, 1)


        self.run_button = createButton("Run", layout, 5, 3, self.runYieldSearch)



    def initGUI(self):

        BASEMAP_SCALE = 0.5

        # Pull station information from picks file
        self.s_info, self.s_name, weights = getStationData(self.bam.setup.station_picks_file, Position(0, 0, 0), expectedcols=8)


    def runYieldSearch(self):
        
        pick_list = []
        list_of_stat_names = []

        with open(INPUT_FILE, "r+") as f:
            stat_lines = f.readlines()[1:]
        

        for stat in stat_lines:

            stat = stat.split(',')

            p = Pick(stat)

            pick_list.append(p)

            list_of_stat_names.append(p.code)
        
        number_of_pts = max([list_of_stat_names.count(i) for i in list_of_stat_names])

        print("Number of points to search for: {:}".format(number_of_pts))

 
        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

        multi_args = [pick_list, ref_pos, FILENAME, self.bam]

        with open(FILENAME, "w+") as f:
            f.write("Latitude [N], Longitude [E], Height [m], Time Error [s], Yield [J], Number of Successful Stations, Total Number of Stations \n")

        f_opt, x_opt = pso(multiLoop, search_min, search_max, \
            args=multi_args, processes=6, particle_output=False, swarmsize=20,\
                 maxiter=20)



        print("###########")
        print("DONE!!!")
        print("########-###")
        try:
            import winsound
            duration = 1000  # milliseconds
            freq = 440  # Hz
            winsound.Beep(freq, duration)
        except:
            pass
        # self.timing_plot.ax1.tricontour(lon1_list, lat1_list, t_list, linewidths=0.5, colors='w')
        # cntr = self.timing_plot.ax1.tricontourf(lon1_list, lat1_list, t_list, cmap="RdBu_r")
        # # self.timing_plot.ax1.colorbar(cntr)

        # self.timing_plot.ax2.tricontour(lon2_list, lat2_list, W_list, linewidths=0.5, colors='w')
        # cntr = self.timing_plot.ax2.tricontourf(lon2_list, lat2_list, W_list, cmap="RdBu_r")
        # # self.timing_plot.ax2.colorbar(cntr)
        # self.timing_plot.show()

if __name__ == "__main__":
    pass