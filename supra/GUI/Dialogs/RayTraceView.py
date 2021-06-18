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
from supra.Supracenter.anglescan2 import anglescan
from supra.Supracenter.cyscan5 import cyscan
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
        self.trajmode = createToggle("Plot Trajectory?", layout, 3, width=1, h_shift=2, tool_tip='')


        self.run_trace_button = createButton("Run", layout, 4, 3, self.runRayTrace)
        self.clear_trace_button = createButton("Clear", layout, 5, 3, self.clearRayTrace)
        # _, self.ray_frac = createLabelEditObj('Fraction of Rays to Show', layout, 5, width=1, h_shift=1, tool_tip='', validate='int', default_txt='50')

        _, self.horizontal_tol = createLabelEditObj('Horizontal Tolerance', layout, 6, width=1, h_shift=1, tool_tip='', validate='float', default_txt='330')
        _, self.vertical_tol = createLabelEditObj('Vertical Tolerance', layout, 7, width=1, h_shift=1, tool_tip='', validate='float', default_txt='3000')

        self.pertstog = createToggle("Use Pertubations", layout, 8, width=1, h_shift=2, tool_tip='')


    def clearRayTrace(self):

        self.rtv_graph.ax.clear()
        self.rtv_graph.show()

    def rayTraceFromSource(self, source, clean_mode=False, debug=False):
        
        ### Set up parameters of source

        stat_idx = self.station_combo.currentIndex()
        stat = self.bam.stn_list[stat_idx]
        stat_pos = stat.metadata.position

        lat =   [source.lat, stat_pos.lat]
        lon =   [source.lon, stat_pos.lon]
        elev =  [source.elev, stat_pos.elev]

        sounding, perturbations = self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev, spline=100, ref_time=self.bam.setup.fireball_datetime)

        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

        stat_pos.pos_loc(source)
        source.pos_loc(source)

        source.z = source.elev
        stat_pos.z = stat_pos.elev

        h_tol = float(self.horizontal_tol.text())
        v_tol = float(self.vertical_tol.text())

        ### Ray Trace

        r, tr, f_particle = cyscan(np.array([source.x, source.y, source.z]), np.array([stat_pos.x, stat_pos.y, stat_pos.z]), \
                            sounding, trace=True, plot=False, particle_output=True, debug=False, \
                            wind=True, h_tol=h_tol, v_tol=v_tol, print_times=True)

        if not clean_mode:
            if not self.pertstog.isChecked():
                dx = np.abs(stat_pos.x - source.x)
                dy = np.abs(stat_pos.y - source.y)
                dz = np.abs(stat_pos.z - source.z)
                print("###### RESULTS ######")
                print("Time: {:.4f} s".format(r[0]))
                print("Azimuth: {:.2f} deg from North".format(r[1]))
                print("Takeoff: {:.2f} deg from up".format(r[2]))
                print("Error in Solution {:.2f} m".format(r[3]))
                print("Distance in x: {:.2f} m".format(dx))
                print("Distance in y: {:.2f} m".format(dy))
                print("Distance in z: {:.2f} m".format(dz))
                print("Horizontal Distance: {:.2f} m".format(np.sqrt(dx**2 + dy**2)))
                print("Total Distance: {:.2f} m".format(np.sqrt(dx**2 + dy**2 + dz**2)))
                print("No Winds Time: {:.2f} s".format(np.sqrt(dx**2 + dy**2 + dz**2)/330))
                print("Time Difference: {:.2f} s".format(r[0] - np.sqrt(dx**2 + dy**2 + dz**2)/330))
            else:
                t_array = []
                az_array = []
                tk_array = []
                err_array = []
                t_array.append(r[0])
                az_array.append(r[1])
                tk_array.append(r[2])
                err_array.append(r[3])


        ### Plot begin and end points of ray-trace

        self.rtv_graph.ax.scatter(source.lon, source.lat, source.elev, c='r', marker='*', s=200)
        self.rtv_graph.ax.scatter(stat_pos.lon, stat_pos.lat, stat_pos.elev, c='b', marker='^', s=200)
        
        positions=[]

        ### Plot trace of eigen-ray

        try:
            path_len = 0
            for i in range(len(tr[:, 0])):
                A = Position(0, 0, 0)
                A.x, A.y, A.z = tr[i, 0], tr[i, 1], tr[i, 2]
                if i > 0:
                    path_len += np.sqrt((tr[i, 0] - tr[i-1, 0])**2 + (tr[i, 1] - tr[i-1, 1])**2 + (tr[i, 2] - tr[i-1, 2])**2)

                else:
                    path_len += np.sqrt((tr[i, 0] - 0)**2 + (tr[i, 1] - 0)**2 + (tr[i, 2] - 0)**2)

                A.pos_geo(source)
                positions.append([A.lon, A.lat, A.elev])
            

            print("Total Path Length: {:.2f} m".format(path_len))
            print("Approximate Ray Time: {:.2f} s".format(path_len/330))
            positions = np.array(positions)
            if not clean_mode:

                self.rtv_graph.ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c='b', alpha=0.5)
            self.rtv_graph.ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], c='k')

            err = np.sqrt((stat_pos.x - tr[-1, 0])**2 + (stat_pos.y - tr[-1, 1])**2 + (stat_pos.z - tr[-1, 2])**2)

            h_err = np.sqrt((stat_pos.x - tr[-1, 0])**2 + (stat_pos.y - tr[-1, 1])**2)
            v_err = np.sqrt((stat_pos.z - tr[-1, 2])**2)

            if debug:
                print("Source Height: {:.2f} km - Final Error: {:.2f} m (v) {:.2f} m (h)".format(source.elev/1000, v_err, h_err))

            self.rtv_graph.ax.scatter(positions[-1, 0], positions[-1, 1], positions[-1, 2], c='g')
        except:
            pass


        if not clean_mode:
            for sol in f_particle:
                r = anglescan(np.array([source.x, source.y, source.z]), sol[0], sol[1], sounding, trace=True, debug=False, wind=True)
                tr = np.array(r[1])
                positions=[]
                for i in range(len(tr[:, 0])):
                    A = Position(0, 0, 0)
                    A.x, A.y, A.z = tr[i, 0], tr[i, 1], tr[i, 2]
                    A.pos_geo(source)
                    positions.append([A.lon, A.lat, A.elev])
                positions = np.array(positions)
                self.rtv_graph.ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], alpha=0.3)
                err = np.sqrt((stat_pos.x - tr[-1, 0])**2 + (stat_pos.y - tr[-1, 1])**2 + (stat_pos.z - tr[-1, 2])**2)
                if err <= 1000:
                    self.rtv_graph.ax.scatter(positions[-1, 0], positions[-1, 1], positions[-1, 2], c='g')
                else:
                    self.rtv_graph.ax.scatter(positions[-1, 0], positions[-1, 1], positions[-1, 2], c='r')
            
        if self.pertstog.isChecked() and len(perturbations) > 0:
            for pert in perturbations:
                r, tr, f_particle = cyscan(np.array([source.x, source.y, source.z]), np.array([stat_pos.x, stat_pos.y, stat_pos.z]), \
                                pert, trace=True, plot=False, particle_output=True, debug=False, \
                                wind=True, h_tol=h_tol, v_tol=v_tol, print_times=True)
                t_array.append(r[0])
                az_array.append(r[1])
                tk_array.append(r[2])
                err_array.append(r[3])

            print("###### NOMINAL RESULTS ######")
            print("Time: {:.4f} s".format(t_array[0]))
            print("Azimuth: {:.2f} deg from North".format(az_array[1]))
            print("Takeoff: {:.2f} deg from up".format(tk_array[2]))
            print("Error in Solution {:.2f} m".format(err_array[3]))
            print("### UNCERTAINTIES ###")
            print("Time: {:.4f} - {:.4f} s ({:.4f} s)".format(np.nanmin(t_array), np.nanmax(t_array), np.nanmax(t_array) - np.nanmin(t_array)))
            print("Azimuth: {:.2f} - {:.2f} deg from North ({:.2f} deg)".format(np.nanmin(az_array), np.nanmax(az_array), np.nanmax(az_array) - np.nanmin(az_array)))
            print("Takeoff: {:.2f} - {:.2f} deg from up ({:.2f} deg)".format(np.nanmin(tk_array), np.nanmax(tk_array), np.nanmax(tk_array) - np.nanmin(tk_array)))
            print("Error in Solution {:.2f} - {:.2f} m ({:.2f} m)".format(np.nanmin(err_array), np.nanmax(err_array), np.nanmax(err_array) - np.nanmin(err_array)))

    def runRayTrace(self):


        traj = self.bam.setup.trajectory


        if not self.trajmode.isChecked():
            try:
                source = traj.findGeo(float(self.source_height.text()))

            except ValueError as e:
                if self.prefs.debug:
                    print(printMessage("Error"), " No source height given!")
                errorMessage("Cannot read source height!", 2, detail='{:}'.format(e))

                return None
        
            self.rayTraceFromSource(source, debug=True)
         
        else:
            # define line bottom boundary
            max_height = traj.pos_i.elev
            min_height = traj.pos_f.elev

            points = traj.trajInterp2(div=50, min_p=min_height, max_p=max_height)

            loadingBar('Calculating Station Times: ', 0, len(points))
            for pp, pt in enumerate(points):
                loadingBar('Calculating Station Times: ', pp, len(points))
                source = Position(pt[0], pt[1], pt[2])
                self.rayTraceFromSource(source, clean_mode=True, debug=True)
                loadingBar('Calculating Station Times: ', pp+1, len(points))

            self.rtv_graph.ax.plot([traj.pos_i.lon, traj.pos_f.lon], [traj.pos_i.lat, traj.pos_f.lat], [traj.pos_i.elev, traj.pos_f.elev])

            self.rtv_graph.show()

if __name__ == '__main__':

    pass