import os
import numpy as np
import pickle

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

        self.load_glm_label, self.load_glm_edits, self.load_glm_buton = createFileSearchObj('Load Ray Trace: ', layout, 13, width=1, h_shift=1)
        self.load_glm_buton.clicked.connect(partial(fileSearch, ['CSV (*.csv)'], self.load_glm_edits))
        self.load_glm_buton.clicked.connect(self.procGLM)

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

    def procGLM(self):

        time = []
        lon = []
        lat = []
        energy = []
        print("Loaded: {:}".format(self.load_glm_edits.text()))
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

        self.glm_graph.ax.scatter(lon, lat, label="GLM")
        self.glm_graph.ax.legend()