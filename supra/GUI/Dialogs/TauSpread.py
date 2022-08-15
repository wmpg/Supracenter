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

import multiprocessing

FILENAME = "F:\\Desktop\\Nord_08_08_22.csv"
#INPUT_FILE = "F:\\Documents\\Meteor_Research\\Event\\Romania\\Romania_Yield_Search_Picks.csv"
#INPUT_FILE = "F:\\Documents\\Meteor_Research\\Event\\alaska_picks_2.csv"
INPUT_FILE = "F:\\Documents\\Meteor_Research\\Event\\meteoroids2020_fireball_picks.csv"


def multiLoop(args):


    lat, lon, H, stat_lines, ref_pos, bam = args

    source = Position(lat, lon, H)

    total_resid = 0
    total_stats = 0
    yield_list = []

    for ss, stat in enumerate(stat_lines):

        print("Working...")

        stat = stat.split(',')

        stat_lat = float(stat[3])
        stat_lon = float(stat[4])
        stat_elev = float(stat[5])
        stat_time = float(stat[7])
        stat_delp = float(stat[9])
        stat_period = float(stat[10])

        stat_pos = Position(stat_lat, stat_lon, stat_elev)

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


            #QUICK FIX FOR ALASKA
            break


            #continue


        T = T - stat_time

        total_resid += np.abs(T)
        total_stats += 1

        v = 1/stat_period

        if stat_delp < 0:
            continue

        ### CALCULATE YIELDS

        
        # f_d, g, _, _, path_length, pdr, reed_attenuation = intscan(S, az, tf, sounding, v, wind=True)

        # rf = refractiveFactor(source, stat_pos, sounding, D_ANGLE=1.5)

        # f_d = np.nanmean(f_d)
        # g = np.nanmean(g)
        # rf = np.nanmean(rf)


        #stat_delp = stat_delp/g/rf
        # adj_pres = pres*(P_a/P)**(1/6)

        # yields
        # Z_chem = findScaledDistance([adj_pres/P_a], chemFuncMinimizer)
        
        #Yield_chem = overpressure2YieldKG(stat_delp, H, stat_pos.elev, R)
        # Yield_chem = overpressure2YieldKG(stat_delp, H, 0, R)
        Yield_chem = overpressure2Yield(stat_delp, R, H)
        # Yield_chem = yieldFromScaledDistance(Z_chem, R, getPP0(H), f_T=f_d, mode="kgHE")

        print("Solution to: {:.2f} N {:.2f} {:.2f} km - E = {:.2f} s {:.2E} J".format(lat, lon, H/1000, T, Yield_chem))

        yield_list.append(Yield_chem)

    if total_stats == 0:
        return None

    total_resid = total_resid/total_stats

    with open(FILENAME, "a+") as f:

        for y in yield_list:
            f.write("{:}, {:}, {:}, {:}, {:}, {:}, {:}\n".format(lat, lon, H, total_resid, y, total_stats, len(stat_lines)))

        f.close()

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

        # Station Combo Box
        # Fill station combo box
        # stn_name_list = []
        # for stn in self.bam.stn_list:
        #     stn_name_list.append("{:}-{:}".format(stn.metadata.network, stn.metadata.code))

        # _, self.station_combo = createComboBoxObj('Station', layout, 2, items=stn_name_list, width=1, h_shift=1, tool_tip='')

        # self.overpressure_label, self.overpressure_edits = createLabelEditObj('Overpressure', layout, 3)
        # self.freq_label, self.freq_edits = createLabelEditObj('Dominant Period', layout, 4)

        self.run_button = createButton("Run", layout, 5, 3, self.runYieldSearch)



    def initGUI(self):

        BASEMAP_SCALE = 0.5

        ### Create Basemap
        # resolution c, l, i, h, f


        # self.m1 = Basemap(projection='merc', \
        #     llcrnrlat=np.ceil(self.bam.setup.lat_centre - BASEMAP_SCALE*self.bam.setup.deg_radius),\
        #     urcrnrlat=np.floor(self.bam.setup.lat_centre + BASEMAP_SCALE*self.bam.setup.deg_radius), \
        #     llcrnrlon=np.ceil(self.bam.setup.lon_centre - BASEMAP_SCALE*self.bam.setup.deg_radius), \
        #     urcrnrlon=np.floor(self.bam.setup.lon_centre + BASEMAP_SCALE*self.bam.setup.deg_radius), \
        #     lat_ts=1, \
        #     resolution='l', ax=self.timing_plot.ax1)

        # self.m1.fillcontinents(color='grey', lake_color='aqua')
        # self.m1.drawcountries(color='black')
        # self.m1.drawlsmask(ocean_color='aqua')

        # self.m1.drawparallels(np.arange(self.bam.setup.lat_centre - BASEMAP_SCALE*self.bam.setup.deg_radius, \
        #     self.bam.setup.lat_centre + BASEMAP_SCALE*self.bam.setup.deg_radius, 1), labels=[1,0,0,1], textcolor="white", fmt="%.1f")
        # self.m1.drawmeridians(np.arange(self.bam.setup.lon_centre - BASEMAP_SCALE*self.bam.setup.deg_radius, \
        #     self.bam.setup.lon_centre + BASEMAP_SCALE*self.bam.setup.deg_radius, 1), labels=[1,0,0,1], textcolor="white", rotation="horizontal", fmt="%.1f")

        # self.m2 = Basemap(projection='merc', \
        #     llcrnrlat=np.ceil(self.bam.setup.lat_centre - BASEMAP_SCALE*self.bam.setup.deg_radius),\
        #     urcrnrlat=np.floor(self.bam.setup.lat_centre + BASEMAP_SCALE*self.bam.setup.deg_radius), \
        #     llcrnrlon=np.ceil(self.bam.setup.lon_centre - BASEMAP_SCALE*self.bam.setup.deg_radius), \
        #     urcrnrlon=np.floor(self.bam.setup.lon_centre + BASEMAP_SCALE*self.bam.setup.deg_radius), \
        #     lat_ts=1, \
        #     resolution='l', ax=self.timing_plot.ax2)

        # self.m2.fillcontinents(color='grey', lake_color='aqua')
        # self.m2.drawcountries(color='black')
        # self.m2.drawlsmask(ocean_color='aqua')

        # self.m2.drawparallels(np.arange(self.bam.setup.lat_centre - BASEMAP_SCALE*self.bam.setup.deg_radius, \
        #     self.bam.setup.lat_centre + BASEMAP_SCALE*self.bam.setup.deg_radius, 1), labels=[1,0,0,1], textcolor="white", fmt="%.1f")
        # self.m2.drawmeridians(np.arange(self.bam.setup.lon_centre - BASEMAP_SCALE*self.bam.setup.deg_radius, \
        #     self.bam.setup.lon_centre + BASEMAP_SCALE*self.bam.setup.deg_radius, 1), labels=[1,0,0,1], textcolor="white", rotation="horizontal", fmt="%.1f")

        # Pull station information from picks file
        self.s_info, self.s_name, weights = getStationData(self.bam.setup.station_picks_file, Position(0, 0, 0), expectedcols=8)




    def runYieldSearch(self):
        
        with open(INPUT_FILE, "r+") as f:
            stat_lines = f.readlines()[1:]



        
        N = 15
        RAD = 0.3
        H = 22300
        dh = 5000


        # del_p = float(self.overpressure_edits.text())

        # make net of lats and lons
        # Alaska
        # lats = np.linspace(67.5, 68.1, N)
        # lons = np.linspace(-150.1, -148.9, N)
        # heights = np.linspace(30000, 36000, 61, endpoint=True)

        # Romania
        # lats = np.linspace(45.5, 46.1, 10)
        # lons = np.linspace(26.5, 27.3, 10)
        # heights = np.array([42700])

        # Meteoroids 2020
        lats = np.linspace(60.3, 60.5, 10)
        lons = np.linspace(16.8, 17.2, 15)
        heights = np.linspace(20000, 24000, 41, endpoint=True)

        #New Zealand
        # lats = np.linspace(-42.2, -41.2, 10)
        # lons = np.linspace(174.5, 175.5, 10)
        # heights = np.linspace(30000, 35000, 5, endpoint=True)


        # stat_idx = self.station_combo.currentIndex()
        # stat = self.bam.stn_list[stat_idx]
        # stat_pos = stat.metadata.position

        # stat_name = "{:}-{:}".format(stat.metadata.network, stat.metadata.code)
        # ref_time = 0

        # for ss, s in enumerate(self.s_name):
        #     if s == stat_name:
        #         ref_time = self.s_info[ss][3]


        with open(FILENAME, "w+") as f:
            f.write("Latitude [N], Longitude [E], Height [m], Time Error [s], Yield [J], Number of Successful Stations, Total Number of Stations \n")

        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

        for lat in lats:
            for lon in lons:
                for h in heights:
                    multiLoop([lat, lon, h, stat_lines, ref_pos, self.bam])

        print("###########")
        print("DONE!!!")
        print("########-###")
        import winsound
        duration = 1000  # milliseconds
        freq = 440  # Hz
        winsound.Beep(freq, duration)
        # self.timing_plot.ax1.tricontour(lon1_list, lat1_list, t_list, linewidths=0.5, colors='w')
        # cntr = self.timing_plot.ax1.tricontourf(lon1_list, lat1_list, t_list, cmap="RdBu_r")
        # # self.timing_plot.ax1.colorbar(cntr)

        # self.timing_plot.ax2.tricontour(lon2_list, lat2_list, W_list, linewidths=0.5, colors='w')
        # cntr = self.timing_plot.ax2.tricontourf(lon2_list, lat2_list, W_list, cmap="RdBu_r")
        # # self.timing_plot.ax2.colorbar(cntr)
        # self.timing_plot.show()

if __name__ == "__main__":
    pass