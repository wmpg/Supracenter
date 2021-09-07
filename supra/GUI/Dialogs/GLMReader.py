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

from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT
from supra.Utils.Formatting import *

class glmWindowDialog(QWidget):

    def __init__(self, bam):

        QWidget.__init__(self)
        
        self.bam = bam

        self.buildGUI()
        
      

    def buildGUI(self):
        self.setWindowTitle('Ray-Trace Viewer')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)

        self.glm_graph = MatplotlibPyQT()
      
        self.glm_graph.ax = self.glm_graph.figure.add_subplot(111)

        layout.addWidget(self.glm_graph, 1, 1, 15, 1)

        # self.hvt_graph = MatplotlibPyQT()
        # self.hvt_graph.ax = self.hvt_graph.figure.add_subplot(111)
        # layout.addWidget(self.hvt_graph, 16, 1, 15, 1)

        # stn_name_list = []
        # for stn in self.bam.stn_list:
        #     stn_name_list.append("{:}-{:}".format(stn.metadata.network, stn.metadata.code))

        # _, self.source_height = createLabelEditObj('Source Height Along Trajectory [m]', layout, 1, width=1, h_shift=1, tool_tip='', validate='float')
        # _, self.station_combo = createComboBoxObj('Station', layout, 2, items=stn_name_list, width=1, h_shift=1, tool_tip='')
        # self.trajmode = createToggle("Plot Trajectory?", layout, 3, width=1, h_shift=2, tool_tip='')
        # self.netmode = createToggle("Run Ray Net?", layout, 9, width=1, h_shift=2, tool_tip='')


        # self.run_trace_button = createButton("Run", layout, 4, 3, self.runRayTrace)
        # self.clear_trace_button = createButton("Clear", layout, 5, 3, self.clearRayTrace)
        # # _, self.ray_frac = createLabelEditObj('Fraction of Rays to Show', layout, 5, width=1, h_shift=1, tool_tip='', validate='int', default_txt='50')

        # _, self.horizontal_tol = createLabelEditObj('Horizontal Tolerance', layout, 6, width=1, h_shift=1, tool_tip='', validate='float', default_txt='330')
        # _, self.vertical_tol = createLabelEditObj('Vertical Tolerance', layout, 7, width=1, h_shift=1, tool_tip='', validate='float', default_txt='3000')

        # self.pertstog = createToggle("Use Pertubations", layout, 8, width=1, h_shift=2, tool_tip='')
        # _, self.source_lat = createLabelEditObj('Source Latitude', layout, 10, width=1, h_shift=1, tool_tip='', validate='float')
        # _, self.source_lon = createLabelEditObj('Source Longitude', layout, 11, width=1, h_shift=1, tool_tip='', validate='float')
        
        # self.save_ray = createButton("Export Ray", layout, 12, 3, self.exportRay)

        self.load_glm_label, self.load_glm_edits, self.load_glm_buton = createFileSearchObj('Load GLM: ', layout, 13, width=1, h_shift=1)
        self.load_glm_buton.clicked.connect(partial(fileSearch, ['CSV (*.csv)'], self.load_glm_edits))
        self.load_glm_buton.clicked.connect(self.procGLM)

        self.for_met_sim = createButton("Save For MetSim", layout, 13, 2, self.saveMetSim)


        # self.draw_stat = createButton("Draw Station", layout, 14, 3, self.drawStat)
        # self.draw_src  = createButton("Draw Source", layout, 15, 3, self.drawSrc)
        # self.draw_traj = createButton("Draw Trajectory", layout, 16, 3, self.drawTraj)

        # _, self.draw_beam = createLabelEditObj('Beam Azimuth', layout, 17, width=1, h_shift=1, tool_tip='', validate='float')
        # self.draw_beam_button = createButton("Draw", layout, 17, 4, self.drawBeam)

        # self.hvt_graph.ax.set_xlabel("Time after Source [s]")
        # self.hvt_graph.ax.set_ylabel("Height [km]")

        traj = self.bam.setup.trajectory

        # define line bottom boundary
        max_height = traj.pos_i.elev
        min_height = traj.pos_f.elev

        points = traj.trajInterp2(div=50, min_p=min_height, max_p=max_height)

        for pp, p in enumerate(points):
            if pp == 0:
                self.glm_graph.ax.scatter(p[1], p[0], c="g", label="Source Trajectory")
            else:
                self.glm_graph.ax.scatter(p[1], p[0], c="g")

        self.glm_graph.ax.legend()

    def saveMetSim(self):
        """ Save GLM station as a fake observer camera at the end of the trajectory given 
        at the center of the Earth
        """
        
        file_name = saveFile("csv", note="")

        time, lon, lat, energy, E, T, lc_list, h_list = self.readGLM()

        station_location = latLonAlt2ECEF(np.radians(lat[0]), np.radians(lon[0]), h_list[0]) 


        with open(file_name, 'w+') as f:
            for ii in range(len(time))[1:]:
                point = latLonAlt2ECEF(np.radians(lat[ii]), np.radians(lon[ii]), h_list[ii])

                dx = point[0] - station_location[0]
                dy = point[1] - station_location[1]
                dz = point[2] - station_location[2]
                dh = np.sqrt(dx**2 + dy**2)

                az = np.arctan2(dx, dy)
                ze = np.arctan2(dz, dh)

                f.write("{:}, {:}, {:}\n".format(T[ii], np.degrees(az), np.degrees(ze)))

        print(printMessage("status"), "Output Complete!")
        print(printMessage("info"), "Output as Time [s], Azimuth (North +East), Zenith - use MeasType = 2 in MetSim")
        print(printMessage("info"), "Station Coordinates: {:.4f}N {:.4f}E {:.4f} m".format(lat[0], lon[0], h_list[0]))
        print(printMessage("info"), "Reference Time: {:}".format(self.bam.setup.fireball_datetime))



    def energyConverter(self, energy):
        # See Jenniskens et al 2018


        # Source to GLM satellite distance
        R = 35780000 #m
        #R = 42170000

        # 4 pi r^2 : r - radius of the effective apperature
        r = 0.0095

        geo_f = 4*np.pi*R**2/r

        blackbody_f = 1.018e3

        E = np.array(energy)*geo_f*blackbody_f

        return E

    def timeConverter(self, time):

        # time is seconds since 1970, online docs are wrong! 
        timestamp = datetime.datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0, microsecond=0)

        time_list = []
        for t in time:
            time_list.append((timestamp + datetime.timedelta(seconds=t/1e3) - self.bam.setup.fireball_datetime).total_seconds())

        return time_list

    def readGLM(self):
        """ Returns time, lon, lat, energy given in GLM file
        and converted energy (E), time from reference point (T) and
        magnitude given by Borovicka definition I = 1500*10^(M/-2.5)
        h_list, approximate height
        """


        traj = self.bam.setup.trajectory
        time = []
        lon = []
        lat = []
        energy = []
        print(printMessage("status"), "Loaded: {:}".format(self.load_glm_edits.text()))
        with open(self.load_glm_edits.text(), 'r+') as f:
            for line in f:
                a = line.strip().split(',')
                
                try:
                    float(a[0])
                    time.append(float(a[0]))
                    lon.append(float(a[1]))
                    lat.append(float(a[2]))
                    energy.append(float(a[3]))
                except:
                    continue

        E = self.energyConverter(energy)
        T = self.timeConverter(time)

        h_list = []

        for t in T:
            hhh = traj.approxHeight(t)

            h_list.append(hhh)

        lc_list = []
        for ii in range(len(E) - 1):

            mag = -2.5*np.log(E[ii]/1500)

            lc_list.append(mag)



        return time, lon, lat, energy, E, T, lc_list, h_list

    def procGLM(self):

        time, lon, lat, energy, E, T, lc_list = self.readGLM()

        file_name = saveFile("csv", note="")

        with open(file_name, 'w+') as f:
            f.write("# Station: GLM\n")
            for ll in range(len(lc_list)):
                f.write("{:}, {:}, {:}\n".format(T, h_list/1000, lc_list))

        self.glm_graph.ax.scatter(lon, lat, label="GLM")
        self.glm_graph.ax.legend()