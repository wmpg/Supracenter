import os
import sys
import numpy as np
import pickle
import multiprocessing

from datetime import datetime, timedelta
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
from supra.Lightcurve.light_curve import processLightCurve, readLightCurve
from supra.Supracenter.anglescanrev import anglescanrev

from wmpl.Utils.TrajConversions import jd2Date
from supra.Utils.pso import pso

def heightErr(height, *cyscan_inputs):
    t, traj, stat_pos, sounding = cyscan_inputs

    trial_source = traj.findGeo(height)

    stat_pos.pos_loc(trial_source)
    trial_source.pos_loc(trial_source)

    trial_source.z = trial_source.elev
    stat_pos.z = stat_pos.elev


    ### Ray Trace
    r, _ = cyscan(np.array([trial_source.x, trial_source.y, trial_source.z]), np.array([stat_pos.x, stat_pos.y, stat_pos.z]), \
                    sounding, trace=False, plot=False, particle_output=True, debug=False, \
                    wind=True, h_tol=300, v_tol=3000, print_times=True, processes=1)


    T = r[0] + traj.findTime(trial_source.elev)

    try:
        t_err = (T - t)**2
    except:
        t_err = 999
    print("Height {:.2f} km | Err {:.2f}".format(height[0]/1000, np.sqrt(t_err)))
    return t_err

def determineBackAz(p1, p2, w_spd, w_dir):
    

    # second last point
    x_1, y_1, z_1, t_1 = p1

    # last point
    x_2, y_2, z_2, t_2 = p2

    dt = t_2 - t_1

    # w = [u, v]
    w = np.array([w_spd*np.sin(np.radians(w_dir)), w_spd*np.cos(np.radians(w_dir))])

    # vector of distance ray is pushed by wind
    w_len = w*dt

    # no vertical winds
    D = [x_1 - w_len[0], y_1 - w_len[1], z_1 - 0]

    vec = [x_2 - D[0], y_2 - D[1], z_2 - D[2]]

    # flip around for back azimuth
    ba = (np.degrees(np.arctan2(vec[0], vec[1])) + 180)%360
    tf = (np.degrees(np.arctan2(np.abs(vec[2]), np.abs(np.sqrt(vec[0]**2 + vec[1]**2)))))%360

    return ba, tf


class rtvWindowDialog(QWidget):

    def __init__(self, bam, prefs):

        QWidget.__init__(self)
        
        self.bam = bam
        self.prefs = prefs
        self.buildGUI()
        
        self.current_eigen = None
        self.current_loaded_rays = []

        self.plot_ba_data = []

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

        layout.addWidget(self.rtv_graph, 1, 1, 15, 1)

        self.hvt_graph = MatplotlibPyQT()
        self.hvt_graph.ax = self.hvt_graph.figure.add_subplot(111)
        layout.addWidget(self.hvt_graph, 16, 1, 15, 1)

        stn_name_list = []
        for stn in self.bam.stn_list:
            stn_name_list.append("{:}-{:}".format(stn.metadata.network, stn.metadata.code))

        _, self.source_height = createLabelEditObj('Source Height Along Trajectory [m]', layout, 1, width=1, h_shift=1, tool_tip='', validate='float')
        _, self.station_combo = createComboBoxObj('Station', layout, 2, items=stn_name_list, width=1, h_shift=1, tool_tip='')
        self.trajmode = createToggle("Plot Trajectory?", layout, 3, width=1, h_shift=2, tool_tip='')
        self.netmode = createToggle("Run Ray Net?", layout, 9, width=1, h_shift=2, tool_tip='')


        self.run_trace_button = createButton("Run", layout, 4, 3, self.runRayTrace)
        self.clear_trace_button = createButton("Clear", layout, 5, 3, self.clearRayTrace)
        # _, self.ray_frac = createLabelEditObj('Fraction of Rays to Show', layout, 5, width=1, h_shift=1, tool_tip='', validate='int', default_txt='50')

        _, self.horizontal_tol = createLabelEditObj('Horizontal Tolerance', layout, 6, width=1, h_shift=1, tool_tip='', validate='float', default_txt='330')
        _, self.vertical_tol = createLabelEditObj('Vertical Tolerance', layout, 7, width=1, h_shift=1, tool_tip='', validate='float', default_txt='3000')

        self.pertstog = createToggle("Use Pertubations", layout, 8, width=1, h_shift=2, tool_tip='')
        _, self.source_lat = createLabelEditObj('Source Latitude', layout, 10, width=1, h_shift=1, tool_tip='', validate='float')
        _, self.source_lon = createLabelEditObj('Source Longitude', layout, 11, width=1, h_shift=1, tool_tip='', validate='float')
        
        self.save_ray = createButton("Export Ray", layout, 12, 3, self.exportRay)

        self.load_ray_label, self.load_ray_edits, self.load_ray_buton = createFileSearchObj('Load Ray Trace: ', layout, 13, width=1, h_shift=1)
        self.load_ray_buton.clicked.connect(partial(fileSearch, ['DAT (*.dat)'], self.load_ray_edits))
        self.load_ray_buton.clicked.connect(self.plotLoadedRay)

        self.draw_stat = createButton("Draw Station", layout, 14, 3, self.drawStat)
        self.draw_src  = createButton("Draw Source", layout, 15, 3, self.drawSrc)
        self.draw_traj = createButton("Draw Trajectory", layout, 16, 3, self.drawTraj)

        _, self.draw_beam = createLabelEditObj('Beam Azimuth', layout, 17, width=1, h_shift=1, tool_tip='', validate='float')
        self.draw_beam_button = createButton("Draw", layout, 17, 4, self.drawBeam)
        _, self.draw_exp_time = createLabelEditObj('Expected Arrival Time [s]', layout, 18, width=1, h_shift=1, tool_tip='', validate='float')
        self.trace_rev_button = createButton("Trace Reverse", layout, 18, 4, self.traceRev)    

        self.hvt_graph.ax.set_xlabel("Time after Source [s]")
        self.hvt_graph.ax.set_ylabel("Height [km]")

        self.load_glm_label, self.load_glm_edits, self.load_glm_buton = createFileSearchObj('Load GLM: ', layout, 19, width=1, h_shift=1)
        self.load_glm_buton.clicked.connect(partial(fileSearch, ['CSV (*.csv)'], self.load_glm_edits))
        self.load_glm_buton.clicked.connect(partial(self.procGLM, True))

        self.fireball_datetime_label, self.fireball_datetime_edits = createLabelDateEditObj("GLM Initial Datetime", layout, 20, h_shift=1)
        self.glm2lc = createButton("GLM to Light Curve", layout, 21, 4, self.glm2LC)
          

        self.pol_graph = MatplotlibPyQT()
        self.pol_graph.ax1 = self.pol_graph.figure.add_subplot(211)
        self.pol_graph.ax2 = self.pol_graph.figure.add_subplot(212)
        layout.addWidget(self.pol_graph, 1, 5, 25, 1)

        self.load_baz_label, self.load_baz_edits, self.load_baz_buton = createFileSearchObj('Load Backazimuth: ', layout, 22, width=1, h_shift=1)
        self.load_baz_buton.clicked.connect(partial(fileSearch, ['CSV (*.csv)'], self.load_baz_edits))
        self.load_baz_buton.clicked.connect(self.loadAngleCSV)
        self.height_unc_button = createButton("Height Uncertainty", layout, 23, 4, self.heightUnc)  

    def traceRev(self):
        ### Set up parameters of source


        traj = self.bam.setup.trajectory

        source = traj.findGeo(float(self.source_height.text()))
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

        az = float(self.draw_beam.text())
        D = []
        # Use 25 angles between 90 and 180 deg

        ### To do this more right, calculate D with bad winds, and then use D to find new winds and then recalc D
        ze_list = np.linspace(1, 89, 25)

        for zenith in ze_list:
            # D = anglescanrev(S.xyz, self.azimuth + offset, zenith, sounding, wind=True)
            # D = anglescanrev(S.xyz, (self.azimuth + offset + 180)%360, zenith, sounding, wind=True)

            # Plus 180 for backazimuth
            D.append(anglescanrev(stat_pos.xyz, (az + 180)%360, zenith, sounding, wind=True, trace=True, fix_phi=False))

        # pt, err = finalanglecheck(self.bam, self.bam.setup.trajectory, self.stn.metadata.position, self.azimuth)
        # print(D)

        for ii, trace in enumerate(D):

            tr = []

            for pt in trace:
                x, y, z, T = pt

                A = Position(0, 0, 0)
                A.x, A.y, A.z = x, y, z
                A.pos_geo(source)
                tr.append([A.lat, A.lon, A.elev])
            tr = np.array(tr)
            self.rtv_graph.ax.plot(tr[:, 1], tr[:, 0], tr[:, 2], c="g")
            print("Ray Trace at ze={:.2f} deg: {:.4f}N {:.4f}E {:.3f}km".format(ze_list[ii], tr[-1, 1], tr[-1, 0], tr[-1, 2]/1000))

    def getSource(self):
        traj = self.bam.setup.trajectory
        source = traj.findGeo(float(self.source_height.text()))
        return source

    def getStat(self):
        stat_idx = self.station_combo.currentIndex()
        stat = self.bam.stn_list[stat_idx]
        stat_pos = stat.metadata.position
        return stat_pos

    def getATM(self, source, stat_pos, perturbations=None):

        lat =   [source.lat, stat_pos.lat]
        lon =   [source.lon, stat_pos.lon]
        elev =  [source.elev, stat_pos.elev]

        sounding, perturbations = self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev, spline=100, \
            ref_time=self.bam.setup.fireball_datetime, perturbations=perturbations)

        return sounding, perturbations

    def heightUnc(self):
        # use given height
        # find nominal height within large tolerance around given height (prevent double solutions)
        # run perturbations within this large tolerance, return heights
        # return bar graph of solutions

       
        source = self.getSource()
        stat_pos = self.getStat()

        sounding, perturbations = self.getATM(source, stat_pos, perturbations=0)

        stat_pos.pos_loc(source)
        source.pos_loc(source)

        source.z = source.elev
        stat_pos.z = stat_pos.elev

        given_height = source.elev
        given_time = float(self.draw_exp_time.text())
        
        h_tol = float(self.horizontal_tol.text())
        v_tol = float(self.vertical_tol.text())


        def findHeightFromTime(sounding, t, approx_height):
            traj = self.bam.setup.trajectory

            found = False

            prev_err = 999

            TIME_TOL = 1e-3
            
            BOUNDS = 4000

            indxs = 100

            stat_pos = self.getStat()

            SWARM_SIZE = 100
            MAXITER = 25
            PHIP = 0.5
            PHIG = 0.5
            OMEGA = 0.5
            MINFUNC = 1e-3
            MINSTEP = 1e-2

            search_min = [approx_height - BOUNDS]
            search_max = [approx_height + BOUNDS]

            f_opt, x_opt = pso(heightErr, search_min, search_max, \
                                args=[t, traj, stat_pos, sounding], processes=multiprocessing.cpu_count()-1, particle_output=False, swarmsize=SWARM_SIZE,\
                             maxiter=MAXITER, phip=PHIP, phig=PHIG, \
                             debug=False, omega=OMEGA, minfunc=MINFUNC, minstep=MINSTEP)

            best_height = f_opt[0]


            print("Best Height = {:.5f} km".format(best_height/1000))

            return best_height

        print("Nominal")
        nominal_height = findHeightFromTime(sounding, given_time, given_height)
        pert_height = []

        for pp, pert in enumerate(perturbations):
            print("Perturbation {:}".format(pp + 1))
            temp_h = findHeightFromTime(pert, given_time, given_height)

            pert_height.append(temp_h)

        print("FINAL RESULTS")
        print("Nominal Height {:} km".format(nominal_height/1000))
        # print("Pert Heights {:} km".format(np.array(pert_height)/1000))

    def glm2LC(self):
        t, h, M = self.procGLM(plot=False)

        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = dlg[0]

        with open(file_name, 'w+') as f:
            for ii in range(len(M)):
                f.write("{:}, {:}, {:}\n".format(t[ii], h[ii], M[ii]))
            f.close()

        print("Wrote file to {:}".format(file_name))

    def procGLM(self, plot=True):
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

        initial_time = self.fireball_datetime_edits.dateTime().toPyDateTime()

        utc_times = []
        for tt, t in enumerate(time):
            utc_times.append(initial_time + timedelta(seconds=(t - time[0])/1000))

        rel_time = []

        traj = self.bam.setup.trajectory
        for tt in utc_times:
            rel_time.append((tt - self.bam.setup.fireball_datetime).total_seconds())

        

        MIN_SIZE = 10
        MAX_SIZE = 400


        lon_list = []
        lat_list = []
        height_list = []
        # I changed my mind and did intensity instead
        mag_list = []
        for ii in range(len(energy)):

            traj_h = traj.approxHeight(lat[ii], lon[ii], rel_time[ii])

            traj_pos = traj.findGeo(traj_h)

            lon_list.append(traj_pos.lon)
            lat_list.append(traj_pos.lat)
            height_list.append(traj_pos.elev)
            mag_list.append(energy[ii])

        min_mag = np.nanmin(mag_list)
        max_mag = np.nanmax(mag_list)

        size_scal = (mag_list-min_mag)*(MAX_SIZE-MIN_SIZE)/(max_mag-min_mag)+MIN_SIZE
        if plot:
            self.rtv_graph.ax.scatter(lon_list, lat_list, height_list, s=size_scal, c=mag_list, cmap='magma')


            self.rtv_graph.ax.plot([traj.pos_i.lon, traj.pos_f.lon], [traj.pos_i.lat, traj.pos_f.lat], [traj.pos_i.elev, traj.pos_f.elev])

        else:
            t = rel_time
            h = np.array(height_list)/1000
            M = -2.5*np.log10(mag_list)
            return t, h, M 


    def plothvt(self):
        
        self.hvt_graph.ax.clear()

        max_t = 0
        for ray in self.current_loaded_rays:
            self.hvt_graph.ax.plot(ray[1], ray[0][:, 2]/1000)
            for t in ray[1]:
                if max_t <= t:
                    max_t = t

        ### Wind contour
        if len(self.current_loaded_rays) > 0:
            ray = self.current_loaded_rays[0]
            pos_i = Position(ray[0][0, 0], ray[0][0, 1], ray[0][0, 2])
            pos_f = Position(ray[0][-1, 0], ray[0][-1, 1], ray[0][-1, 2])

            sounding, perturbations = self.bam.atmos.getSounding(lat=[pos_f.lat, pos_i.lat],\
                     lon=[pos_f.lon, pos_i.lon], heights=[pos_f.elev, pos_i.elev], spline=100, \
                     ref_time=self.bam.setup.fireball_datetime)

            levels = sounding[:, 0]/1000

            # maximum possible speed
            speed = sounding[:, 1] + sounding[:, 2]            

            xmin, xmax = self.hvt_graph.ax.get_xlim()
            ymin, ymax = self.hvt_graph.ax.get_ylim()

            h = []
            t = []
            z = []

            for i in range(len(levels)):
                h.append(levels[i])
                t.append(0)
                z.append(speed[i])

            for i in range(len(levels)):
                h.append(levels[i])
                t.append(max_t)
                z.append(speed[i])


            sc = self.hvt_graph.ax.tricontourf(t, h, z, alpha=0.3, cmap='viridis')
            if len(self.current_loaded_rays) == 1:
                


                cbar = self.hvt_graph.figure.colorbar(sc, pad=0.2)
                cbar.ax.set_xlabel("Maximum Effective Speed [m/s]")


        self.hvt_graph.ax.set_xlabel("Time after Source [s]")
        self.hvt_graph.ax.set_ylabel("Height [km]")

        self.hvt_graph.show()



    def drawBeam(self):
        # For short distances this is an ok approx
        stat_idx = self.station_combo.currentIndex()
        stat = self.bam.stn_list[stat_idx]
        stat_pos = stat.metadata.position

        start_point = [stat_pos.lon, stat_pos.lat]

        y = np.tan(np.radians(float(self.draw_beam.text())))

        # deg lat/lon to extend beam to
        SCALE = 1

        a = SCALE/np.sqrt(y**2 + 1)

        new_pos = [stat_pos.lon + a, stat_pos.lat + a*y]

        self.rtv_graph.ax.plot([start_point[0], new_pos[0]], [start_point[1], new_pos[1]], [0, 0])
        self.pol_graph.ax1.axhline(y=float(self.draw_beam.text()), linestyle='-')


    def drawStat(self):
        stat_idx = self.station_combo.currentIndex()
        stat = self.bam.stn_list[stat_idx]
        stat_pos = stat.metadata.position
        
        self.rtv_graph.ax.scatter(stat_pos.lon, stat_pos.lat, stat_pos.elev, marker="^", c="b", s=200)


    def drawSrc(self):

        traj = self.bam.setup.trajectory
        if len(self.source_lat.text()) != 0 and len(self.source_lon.text()) != 0:
            source = Position(float(self.source_lat.text()), float(self.source_lon.text()), float(self.source_height.text())) 
        else:
            try:
                source = traj.findGeo(float(self.source_height.text()))

            except ValueError as e:
                if self.prefs.debug:
                    print(printMessage("Error"), " No source height given!")
                errorMessage("Cannot read source height!", 2, detail='{:}'.format(e))

                return None
        self.rtv_graph.ax.scatter(source.lon, source.lat, source.elev, marker="*", c="r", s=200)


    def drawTraj(self):

        traj = self.bam.setup.trajectory

        MIN_SIZE = 10
        MAX_SIZE = 400

        if len(self.bam.setup.light_curve_file) > 0 or not hasattr(self.bam.setup, "light_curve_file"):


            light_curve = readLightCurve(self.bam.setup.light_curve_file)

            light_curve_list = processLightCurve(light_curve)


            lon_list = []
            lat_list = []
            height_list = []
            # I changed my mind and did intensity instead
            mag_list = []
            for L in light_curve_list:
                for ii in range(len(L.h)):
                    traj_pos = traj.findGeo(L.h[ii]*1000)

                    lon_list.append(traj_pos.lon)
                    lat_list.append(traj_pos.lat)
                    height_list.append(L.h[ii]*1000)
                    mag_list.append(L.I[ii])

            min_mag = np.nanmin(mag_list)
            max_mag = np.nanmax(mag_list)

            size_scal = (mag_list-min_mag)*(MAX_SIZE-MIN_SIZE)/(max_mag-min_mag)+MIN_SIZE
            self.rtv_graph.ax.scatter(lon_list, lat_list, height_list, s=size_scal, c=mag_list, cmap='magma')


        self.rtv_graph.ax.plot([traj.pos_i.lon, traj.pos_f.lon], [traj.pos_i.lat, traj.pos_f.lat], [traj.pos_i.elev, traj.pos_f.elev])


    def exportRay(self):

        filename = QFileDialog.getSaveFileName(self, 'Save File', '', 'Dat File (*.dat)')

        with open(str(filename[0]), 'w+') as f:
            
            f.write("# lat [deg]    lon [deg]   z [km]  geo. atten. [dB]    absorption [dB] time [s]\n")

            for i in range(len(self.current_eigen[0])):

                position = self.current_eigen[0]
                time = self.current_eigen[1]
                f.write("{:}, {:}, {:}, {:}, {:}, {:}\n".format(position[i, 1], \
                            position[i, 0], position[i, 2]/1000, -999, -999, time[i]))
        
        errorMessage("Ray Exported!", 0, title="Everybody Loves Ray-Tracing", detail='File saved to {:}'.format(filename))

    def plotLoadedRay(self):

        ray_pos = []
        ray_time = []
        with open(self.load_ray_edits.text(), 'r+') as f:
            for line in f:
                a = line.strip().split(',')
                
                if len(a) == 1:
                    a = line.strip().split("\t")

                if line[0] == "#":
                    continue
                lat = float(a[0])
                lon = float(a[1])
                elev = float(a[2])*1000
                t = float(a[5])

                ray_pos.append([lat, lon, elev])
                ray_time.append(t)

        ray_pos = np.array(ray_pos)

        self.rtv_graph.ax.plot(ray_pos[:, 1], ray_pos[:, 0], ray_pos[:, 2])
        self.current_loaded_rays.append([ray_pos, ray_time])
        self.plothvt()

    def clearRayTrace(self):

        self.rtv_graph.ax.clear()
        self.rtv_graph.show()

    def rayTraceFromSource(self, source, clean_mode=False, debug=False):

        traj = self.bam.setup.trajectory

        ### Set up parameters of source

        stat_idx = self.station_combo.currentIndex()
        stat = self.bam.stn_list[stat_idx]
        stat_pos = stat.metadata.position

        lat =   [source.lat, stat_pos.lat]
        lon =   [source.lon, stat_pos.lon]
        elev =  [source.elev, stat_pos.elev]

        if not clean_mode:
            print("Source Location")
            print(source)
            print("Station Location")
            print(stat_pos)

        sounding, perturbations = self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev, spline=1000, ref_time=self.bam.setup.fireball_datetime)

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
                            wind=True, h_tol=h_tol, v_tol=v_tol, print_times=True, processes=1)

        if not clean_mode:
            
            az = np.radians(r[1])
            tf = np.radians(180 - r[2])
            
            try:
                u = np.array([traj.vector.x,
                            traj.vector.y,
                            traj.vector.z])
            except AttributeError:
                # If no trajectory data is available
                pass

            v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])
            
            try:
                angle_off = abs(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))) - 90)
            except:
                pass

            if not self.pertstog.isChecked():
                dx = np.abs(stat_pos.x - source.x)
                dy = np.abs(stat_pos.y - source.y)
                dz = np.abs(stat_pos.z - source.z)
                try:
                    time_along_trajectory = traj.findTime(source.elev)
                except:
                    pass

                print("###### RESULTS ######")
                print("Acoustic Path Time: {:.4f} s".format(r[0]))

                try:
                    print("Time Along Trajectory (wrt Reference): {:.4f} s".format(time_along_trajectory))
                    print("Total Time from Reference: {:.4f} s".format(r[0] + time_along_trajectory))
                    print("Launch Angle {:.2f} deg".format(angle_off))
                except:
                    pass

                print("###### EXTRAS #######")
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
                
                try:
                    print("Time Along Trajectory: {:.4f} s".format(time_along_trajectory))
                    print("Total Time from Reference: {:.4f} s".format(r[0] + time_along_trajectory))
                except:
                    pass

            else:
                t_array = []
                az_array = []
                tk_array = []
                err_array = []
                angle_off_array = []
                t_array.append(r[0])
                az_array.append(r[1])
                tk_array.append(r[2])
                err_array.append(r[3])
                try:
                    angle_off_array.append(angle_off)
                except:
                    pass

        try:

            N_LAYERS = 1

            for i in range(N_LAYERS):
                try:
                    ba, tf = determineBackAz(tr[-(i+2), :], tr[-1, :], sounding[0, 2], np.degrees(sounding[0, 3]))

                    if hasattr(self, "plot_ba_data"):
                        self.plot_ba_data.append([source.elev/1000, ba, tf])
                except IndexError:
                    pass


        except TypeError:
            pass
        
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
            
            if not clean_mode:
                print("Total Path Length: {:.2f} m".format(path_len))
                print("Approximate Ray Time: {:.2f} s".format(path_len/330))
            positions = np.array(positions)
            if not clean_mode:

                self.rtv_graph.ax.scatter(positions[:, 0], positions[:, 1], positions[:, 2], c='b', alpha=0.5)
            self.rtv_graph.ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], c='k')

            self.current_eigen = [positions, tr[:, -1]]

            err = np.sqrt((stat_pos.x - tr[-1, 0])**2 + (stat_pos.y - tr[-1, 1])**2 + (stat_pos.z - tr[-1, 2])**2)

            h_err = np.sqrt((stat_pos.x - tr[-1, 0])**2 + (stat_pos.y - tr[-1, 1])**2)
            v_err = np.sqrt((stat_pos.z - tr[-1, 2])**2)

            if debug:
                print("Source Height: {:.2f} km - Final Error: {:.2f} m (v) {:.2f} m (h)".format(source.elev/1000, v_err, h_err))

            self.rtv_graph.ax.scatter(positions[-1, 0], positions[-1, 1], positions[-1, 2], c='g')

            self.current_loaded_rays.append([positions, tr[:, -1]])
            self.plothvt()
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
                # self.rtv_graph.ax.plot(positions[:, 0], positions[:, 1], positions[:, 2], alpha=0.3)
                err = np.sqrt((stat_pos.x - tr[-1, 0])**2 + (stat_pos.y - tr[-1, 1])**2 + (stat_pos.z - tr[-1, 2])**2)
                if err <= 1000:
                    self.rtv_graph.ax.scatter(positions[-1, 0], positions[-1, 1], positions[-1, 2], c='g')
                else:
                    self.rtv_graph.ax.scatter(positions[-1, 0], positions[-1, 1], positions[-1, 2], c='r')
            
        if self.pertstog.isChecked() and len(perturbations) > 0:
            if clean_mode:
                t_array = []
                az_array = []
                tk_array = []
                err_array = []
                angle_off_array = []
            for pert_idx, pert in enumerate(perturbations):
                sys.stdout.write("\r Working on Perturbation {:}/{:}".format(pert_idx+1, len(perturbations)))
                sys.stdout.flush()
                r, tr, f_particle = cyscan(np.array([source.x, source.y, source.z]), np.array([stat_pos.x, stat_pos.y, stat_pos.z]), \
                                pert, trace=True, plot=False, particle_output=True, debug=False, \
                                wind=True, h_tol=h_tol, v_tol=v_tol, print_times=True)
                

                t_array.append(r[0])
                az_array.append(r[1])
                tk_array.append(r[2])
                err_array.append(r[3])


                az = np.radians(r[1])
                tf = np.radians(180 - r[2])
                u = np.array([traj.vector.x,
                            traj.vector.y,
                            traj.vector.z])
                v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])
                angle_off = abs(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))) - 90)
                angle_off_array.append(angle_off)

                try:

                    ba = determineBackAz(tr[-2, :], tr[-1, :], pert[0, 2], np.degrees(pert[0, 3]))

                    # print("Height: {:.2f} km".format(source.elev/1000))
                    # print("Back Azimuth: {:.2f} deg".format(ba))
                    # print("Travel Time: {:.2f} s".format(r[0]))
                    # print("Approx: {:.2f} deg".format(last_azimuth))
                    # print("Pure Winds: {:.2f} deg".format(np.degrees(sounding[0, 3])%360))
                    self.plot_ba_data.append([source.elev/1000, ba, r[0]])

                except TypeError:
                    pass
        
            print()
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
            print("Angle Off {:.2f} - {:.2f} deg ({:.2f} deg)".format(np.nanmin(angle_off_array), np.nanmax(angle_off_array), np.nanmax(angle_off_array) - np.nanmin(angle_off_array)))

            print("Saving CSV of Perturbations")
            file_name = saveFile("csv", note="Perturbations")

            with open(file_name, "w+") as f:
                f.write("Height [km], Time [s], Azimuth [deg from North], Takeoff [deg from Up], Error [m], Angle Off [deg]\n")
                for ll, line in enumerate(t_array):
                    if ll == len(t_array):
                        f.write("{:}, {:}, {:}, {:}, {:}, {:}".format(source.elev/1000, t_array[ll], az_array[ll], tk_array[ll], err_array[ll], angle_off_array[ll]))
                    else:
                        f.write("{:}, {:}, {:}, {:}, {:}, {:}\n".format(source.elev/1000, t_array[ll], az_array[ll], tk_array[ll], err_array[ll], angle_off_array[ll]))

    def loadAngleCSV(self):
        
        times = []
        height = []
        trace = []
        incl = []
        baz = []

        with open(self.load_baz_edits.text(), 'r+') as f:
            for line in f:
                a = line.strip().split(',')
                
                times.append(float(a[0]))
                height.append(float(a[1]))
                trace.append(float(a[2]))
                incl.append(float(a[3]))
                baz.append(float(a[4]))

        self.pol_graph.ax1.scatter(np.array(height)/1000, np.array(baz))
        self.pol_graph.ax1.set_xlabel("Height [km]")
        self.pol_graph.ax1.set_ylabel("Backazimuth [deg]")

        self.pol_graph.ax2.scatter(np.array(height)/1000, np.array(incl))
        self.pol_graph.ax2.set_xlabel("Height [km]")
        self.pol_graph.ax2.set_ylabel("Inclination [deg]")



    def runRayTrace(self):


        traj = self.bam.setup.trajectory

        try:
            source_lat = float(self.source_lat.text())
            source_lon = float(self.source_lon.text())
            print("Overwritting source with source lat and lon")
        except:
            source_lat = None
            source_lon = None

        if not self.trajmode.isChecked():
            if self.netmode.isChecked():

                daz = 0.05
                # az_net = np.arange(0, 360 - daz, daz)
                az_net = np.arange(105, 160, daz)

                dtf = 0.05
                tf_net = np.arange(90 + dtf, 115 - dtf, dtf)

                stat_idx = self.station_combo.currentIndex()
                stat = self.bam.stn_list[stat_idx]
                stat_pos = stat.metadata.position
                source = traj.findGeo(float(self.source_height.text()))

                if source_lat is not None:
                    source.lat = source_lat
                    source.lon = source_lon

                h_tol = float(self.horizontal_tol.text())
                v_tol = float(self.vertical_tol.text())
                lat =   [source.lat, source.lat]
                lon =   [source.lon, source.lon]
                elev =  [source.elev, 0]
                sounding, perturbations = self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev, spline=1000, ref_time=self.bam.setup.fireball_datetime)

                ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

                source.pos_loc(source)

                source.z = source.elev
                stat_pos.pos_loc(source)
                stat_pos.z = stat_pos.elev

                az_list = []
                tf_list = []
                T_list = []
                h_err_list = []
                v_err_list = []

                for az in az_net:
                    for tf in tf_net:
                        r = anglescan(np.array([source.x, source.y, source.z]), az, tf, sounding, trace=False, debug=False, wind=True)
                        x, y, z, T = r

                        h_err = np.sqrt((x - stat_pos.x)**2 + (y - stat_pos.y)**2)
                        v_err = np.abs(z - stat_pos.z)

                        if h_err <= h_tol and v_err <= v_tol:
                            eigen = True
                        else:
                            eigen = False
                        print("Ray Trace - Azimuth: {:.2f} Takeoff: {:.2f} Time: {:.2f} s H_err: {:.2f} V_err: {:.2f}".format(az, tf, T, h_err, v_err))
 
                        if eigen == True:
                            az_list.append(az)
                            tf_list.append(tf)
                            T_list.append(T)
                            h_err_list.append(h_err)
                            v_err_list.append(v_err)


                print("Done")
                
                print("Arrival Range: {:.2f}-{:.2f} s".format(np.nanmin(T_list), np.nanmax(T_list)))
                print("Azimuths: {:.2f}-{:.2f}".format(az_list[np.nanargmin(T_list)], az_list[np.nanargmax(T_list)]))
                print("Takeoffs: {:.2f}-{:.2f}".format(tf_list[np.nanargmin(T_list)], tf_list[np.nanargmax(T_list)]))
            else:

                if len(self.source_lat.text()) != 0 and len(self.source_lon.text()) != 0:
                    source = Position(float(self.source_lat.text()), float(self.source_lon.text()), float(self.source_height.text())) 
                else:
                    try:
                        source = traj.findGeo(float(self.source_height.text()))

                    except ValueError as e:
                        if self.prefs.debug:
                            print(printMessage("Error"), " No source height given!")
                        errorMessage("Cannot read source height!", 2, detail='{:}'.format(e))

                        return None
                
                self.rayTraceFromSource(source, debug=True)
        
        else:
            
            self.plot_ba_data = []
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

            self.drawTraj()

            self.rtv_graph.show()

            self.plot_ba_data = np.array(self.plot_ba_data)



            self.pol_graph.ax1.scatter(self.plot_ba_data[:, 0], self.plot_ba_data[:, 1])
            self.pol_graph.ax1.set_xlabel("Height [km]")
            self.pol_graph.ax1.set_ylabel("Backazimuth [deg]")

            self.pol_graph.ax2.scatter(self.plot_ba_data[:, 0], self.plot_ba_data[:, 2])
            self.pol_graph.ax2.set_xlabel("Height [km]")
            self.pol_graph.ax2.set_ylabel("Inclination [deg]")
            

            print("Back Azimuths and Zeniths:")
            for i in range(len(self.plot_ba_data)):
                print("{:.4f} km, {:.2f} deg, {:.2f} deg".format(self.plot_ba_data[i, 0], \
                                                                self.plot_ba_data[i, 1], \
                                                                self.plot_ba_data[i, 2]))
            # self.pol_graph.ax2.axhline(y=float(self.draw_exp_time.text()), linestyle='-')
            
            # self.pol_graph.ax1.axvline(x=float(self.source_height.text()), linestyle='-')
            # self.pol_graph.ax2.axvline(x=float(self.source_height.text()), linestyle='-')

if __name__ == '__main__':

    pass