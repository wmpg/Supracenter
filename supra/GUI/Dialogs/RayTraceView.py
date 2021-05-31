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
from supra.Supracenter.cyscan2 import cyscan
from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT

from supra.Utils.Formatting import *



class rtvWindowDialog(QWidget):

    def __init__(self, bam, prefs):

        QWidget.__init__(self)
        
        self.bam = bam
        self.prefs = prefs
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

        self.rtv_graph = MatplotlibPyQT()
      
        self.rtv_graph.ax = self.rtv_graph.figure.add_subplot(111, projection='3d')

        layout.addWidget(self.rtv_graph, 1, 1, 100, 1)

        stn_name_list = []
        for stn in self.bam.stn_list:
            stn_name_list.append("{:}-{:}".format(stn.metadata.network, stn.metadata.code))

        _, self.source_height = createLabelEditObj('Source Height Along Trajectory [m]', layout, 1, width=1, h_shift=1, tool_tip='', validate='float')
        _, self.station_combo = createComboBoxObj('Station', layout, 2, items=stn_name_list, width=1, h_shift=1, tool_tip='')

        self.run_trace_button = createButton("Run", layout, 3, 3, self.runRayTrace)

    def runRayTrace(self):

        traj = self.bam.setup.trajectory

        try:
            source = traj.findGeo(float(self.source_height.text()))
        except ValueError as e:
            if self.prefs.debug:
                print(printMessage("Error"), " No source height given!")
            errorMessage("Cannot read source height!", 2, detail='{:}'.format(e))

            return None

        stat_idx = self.station_combo.currentIndex()
        stat = self.bam.stn_list[stat_idx]
        stat_pos = stat.metadata.position

        lat =   [source.lat, stat_pos.lat]
        lon =   [source.lon, stat_pos.lon]
        elev =  [source.elev, stat_pos.elev]

        sounding, _ = self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev)

        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

        stat_pos.pos_loc(ref_pos)
        source.pos_loc(ref_pos)

        t_arrival, azimuth, takeoff, E, trace_list = cyscan(np.array([source.x, source.y, source.z]), np.array([stat_pos.x, stat_pos.y, stat_pos.z]), sounding, \
                wind=self.prefs.wind_en, n_theta=self.prefs.pso_theta, n_phi=self.prefs.pso_phi,
                h_tol=self.prefs.pso_min_ang, v_tol=self.prefs.pso_min_dist, trace=True)

        print(trace_list)

        x_list, y_list, z_list = [], [], []
        for h in range(len(trace_list)):
            x_list.append(trace_list[h][0])
            y_list.append(trace_list[h][1])
            z_list.append(trace_list[h][2])

            self.rtv_graph.ax.plot(x_list, y_list, z_list)

        self.rtv_graph.show()