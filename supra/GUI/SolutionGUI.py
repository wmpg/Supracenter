
import os
import time
import datetime
import copy
import webbrowser
import pickle

from netCDF4 import Dataset
from PyQt5.QtWidgets import *
import sys
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
from matplotlib.colors import Normalize
import matplotlib.pyplot as plt
import numpy as np
import obspy
import scipy.signal

from functools import partial

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Fireballs.GetIRISData import readStationAndWaveformsListFile, butterworthBandpassFilter, convolutionDifferenceFilter
from supra.Fireballs.SeismicTrajectory import parseWeather, timeOfArrival, getStationList, estimateSeismicTrajectoryAzimuth, plotStationsAndTrajectory

from supra.Supracenter.slowscan import slowscan
from supra.Supracenter.stationDat import convStationDat
from supra.Supracenter.psoSearch import psoSearch
from supra.Supracenter.netCDFconv import findECMWFSound
from supra.Supracenter.SPPT import perturb as perturbation_method
from supra.Supracenter.fetchCopernicus import copernicusAPI
from supra.Supracenter.cyscan import cyscan
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.SPPT import perturb

from supra.GUI.GUITools import *

from wmpl.Utils.TrajConversions import datetime2JD, jd2Date
from wmpl.Utils.Earth import greatCircleDistance

from supra.Utils.AngleConv import loc2Geo, chauvenet
from supra.Utils.Formatting import *
from supra.Utils.Classes import Position, Station, Config, Supracenter, Trajectory, Constants, Pick
from supra.Utils.TryObj import tryBool, tryFloat, tryInt, tryPosition, trySupracenter, tryAngle, tryTrajectory, saveDefaults

global arrTimes 
global sounding

DATA_FILE = 'data.txt'
OUTPUT_CSV = 'data_picks.csv'
  
class SolutionGUI(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setup = Config()

        self._main = QWidget()
        self.setCentralWidget(self._main)
        layout = QGridLayout(self._main)

        self.setWindowTitle('WMPG Seismic and Infrasound Meteor Program')
        app_icon = QtGui.QIcon()
        app_icon.addFile('wmpl.png', QtCore.QSize(16,16))
        self.setWindowIcon(app_icon)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        self.slider_scale = 0.25
        self.bandpass_scale = 0.1

        self.tab_widget = QTabWidget()
        self.tab_widget.blockSignals(True)

        self.ini_dock = QDockWidget("Variables", self)
        self.addDockWidget(Qt.LeftDockWidgetArea, self.ini_dock)
        self.ini_dock.setFeatures(QtGui.QDockWidget.DockWidgetFloatable | QtGui.QDockWidget.DockWidgetMovable)
        

        self.addIniDockWidgets()
        self.addPicksReadWidgets()
        self.addSupraWidgets()
        self.addSupWidgets()
        self.addMakePicksWidgets()
        self.addSeisTrajWidgets()
        self.addFetchATMWidgets()
        self.addProfileWidgets()
        self.addRayTracerWidgets()

        self.var_typ = 't'

        self.tab_widget.blockSignals(False)
        layout.addWidget(self.tab_widget, 1, 1)

        menu_bar = self.menuBar() 
        layout.addWidget(menu_bar, 0, 1)
        file_menu = menu_bar.addMenu('&File')
        about_menu = menu_bar.addMenu('&About')
        view_menu = menu_bar.addMenu('&View')

        file_exit = QAction("Exit", self)
        file_exit.setShortcut('Ctrl+Q')
        file_exit.setStatusTip('Exit application')
        file_exit.triggered.connect(self.quitApp)
        file_menu.addAction(file_exit)

        file_qsave = QAction("Quick Save", self)
        file_qsave.setShortcut('Ctrl+S')
        file_qsave.setStatusTip('Saves setup file')
        file_qsave.triggered.connect(partial(self.saveINI, True, True))
        file_menu.addAction(file_qsave)

        file_save = QAction("Save", self)
        file_save.setShortcut('Ctrl+Shift+S')
        file_save.setStatusTip('Saves setup file')
        file_save.triggered.connect(partial(self.saveINI, True))
        file_menu.addAction(file_save)

        file_load = QAction("Load", self)
        file_load.setShortcut('Ctrl+L')
        file_load.setStatusTip('Loads setup file')
        file_load.triggered.connect(self.loadINI)
        file_menu.addAction(file_load)

        about_github = QAction("GitHub", self)
        about_github.triggered.connect(self.openGit)
        about_menu.addAction(about_github)

        about_docs = QAction("Documentation", self)
        about_docs.triggered.connect(self.openDocs)
        about_menu.addAction(about_docs)

        view_vartools = QAction("Show/Hide Toolbar", self)
        view_vartools.setShortcut('Ctrl+V')
        view_vartools.setStatusTip('Toggle if the variable toolbar is visible')
        view_vartools.triggered.connect(self.viewToolbar)
        view_menu.addAction(view_vartools)

        stylesheet = """ 
        QTabWidget>QWidget>QWidget{background: gray;}
        QLabel{color: white;}
        QCheckBox{color: white;}
        QDockWidget{color: white; background: black;}
        QGroupBox{color: white;}
        QGroupBox{ 
        border: 2px white; 
        border-radius: 0px; }
        QMessageBox{color: white; background: black;} 
        """

        self.setStyleSheet(stylesheet)

        pg.setConfigOptions(antialias=True)
        self.ray_pick = pg.ScatterPlotItem()
        self.ray_pick_traj = pg.ScatterPlotItem()
        self.ray_pick_point = [0, 0, 0]
        self.ctrl_pressed = False

    def viewToolbar(self):
        self.ini_dock.toggleViewAction().trigger()

    def quitApp(self):

        reply = QMessageBox.question(self, 'Quit Program', 
                 'Are you sure you want to quit?', QMessageBox.Yes, QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            qApp.quit()
        else:
            return None
        

    def openGit(self):
        webbrowser.open_new_tab("https://github.com/dvida/Supracenter")

    def openDocs(self):
        webbrowser.open_new_tab("supra/Fireballs/docs/index.html")

    def addPicksReadWidgets(self):
        picks_read_tab = QWidget()
        picks_read_tab_content = QGridLayout()
        picks_read_tab.setLayout(picks_read_tab_content)

        self.tab_widget.addTab(picks_read_tab, "Picks Read")

        self.csv_table = QTableWidget(0, 9)
        picks_read_tab_content.addWidget(self.csv_table, 1, 1, 1, 4)
        self.csv_table.setHorizontalHeaderLabels(['Pick Group', 'Network', 'Code', 'Latitude', 'Longitude', 'Elevation', 'Pick JD', 'Pick Time', 'station_number'])
        header = self.csv_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)

        self.csv_table_add = QPushButton("+")
        picks_read_tab_content.addWidget(self.csv_table_add, 2, 2, 1, 1)
        self.csv_table_add.clicked.connect(partial(self.changeRows, self.csv_table, 1))
        self.csv_table_add.setToolTip("Add row")

        self.csv_table_min = QPushButton("-")
        picks_read_tab_content.addWidget(self.csv_table_min, 2, 1, 1, 1)
        self.csv_table_min.clicked.connect(partial(self.changeRows, self.csv_table, -1))
        self.csv_table_min.setToolTip("Remove row")

        self.csv_table_load = QPushButton("Load")
        picks_read_tab_content.addWidget(self.csv_table_load, 2, 3, 1, 1)
        self.csv_table_load.clicked.connect(self.csvLoad)

        self.csv_table_save = QPushButton("Save")
        picks_read_tab_content.addWidget(self.csv_table_save, 2, 4, 1, 1)
        self.csv_table_save.clicked.connect(self.csvSave)

    def csvLoad(self):
        
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setNameFilters(['CSV File (*.csv)'])
        dlg.exec_()

        filename = dlg.selectedFiles()

        data_table = []

        with open(filename[0]) as f:
            
            next(f)
            
            for line in f:
                data_table.append(line.split(','))

        toTable(self.csv_table, data_table)


    def csvSave(self):

        dlg = QFileDialog.getSaveFileName(self, 'Save File')


        if '.csv' not in dlg[0]:
            file_name = dlg[0] + '.csv'
        else:
            file_name = dlg[0]

        data_set = fromTable(self.csv_table)
        # Open the output CSV
        with open(os.path.join(file_name), 'w') as f:

            # Write the header
            f.write('Pick group, Network, Code, Lat, Lon, Elev, Pick JD, Pick time, station_number \n')

            # Go through all picks
            for line in data_set:
                line[-1] = int(line[-1])
                # Write the CSV entry
                f.write("{:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}\n".format(*line))

        errorMessage('Output to CSV!', 0, title='Exported!')


    def supSaveChanges(self):
        pass

    def supraSaveChanges(self):

        self.sounding_file_edits.setText(self.atmospheric_file_edit.text())
        self.station_picks_edits.setText(self.picks_file_edit.text())
        self.fireball_datetime_edits.setDateTime(self.ref_edit.dateTime().toPyDateTime())
        self.lat_frag_edits.setText(str(self.lat_edit.text()))
        self.lon_frag_edits.setText(str(self.lon_edit.text()))
        self.elev_frag_edits.setText(str(self.elev_edit.text()))
        self.time_frag_edits.setText(str(self.time_edit.text()))

        errorMessage('File Copied!', 0, title='Ini Information Copied')
   
    def seisSearch(self):

        # Read station file
        station_list = getStationList(os.path.join(self.setup.working_directory, self.setup.fireball_name, self.setup.station_picks_file))

        # Extract all Julian dates
        jd_list = [entry[6] for entry in station_list]
        # Calculate the arrival times as the time in seconds from the earliest JD
        jd_ref = min(jd_list)
        ref_indx = np.argmin(jd_list)

        try:
            _, _, _, lat0, lon0, elev0, ref_time, pick_time, station_no = station_list[ref_indx]

        except:
            errorMessage("Old data_picks.csv file detected!", 2, detail='data_picks.csv files created previous to Jan. 8, 2019 are lacking a channel tag added. Redownloading the waveform files will likely fix this')

        # Date for weather
        ref_time = jd2Date(jd_ref)
        self.setup.ref_time = datetime.datetime(*(map(int, ref_time)))

        sounding = parseWeather(self.setup)

        # Set up search parameters
        p0 = [self.setup.lat_f, self.setup.lon_f, self.setup.t0, self.setup.v, self.setup.azimuth, self.setup.zenith]

        # Set up perturbed array
        if self.setup.perturb == True:

            try:

                allTimes = np.load(self.setup.arrival_times_file)
                if self.setup.debug:
                    print("Status: Loaded picks from perturbations")


            except:
                errorMessage("Unable to find perturbed times file", 2, detail='Place a file "all_pick_times.npy" generated by MakeIRISPicks into the file_name indicated in the SeismicTrajectory.ini file')
        else:
            allTimes = []

        self.setup.run_mode = self.setup.run_mode.lower()

        fig = plt.figure(figsize=plt.figaspect(0.5))
        plt.style.use('dark_background')
        fig.set_size_inches(8, 5)
        self.seis_traj_ax = fig.add_subplot(1, 1, 1, projection='3d')

        # If searching
        if self.setup.run_mode == 'search':                
            sup, errors, results = estimateSeismicTrajectoryAzimuth(station_list, self.setup, sounding, p0=p0, \
                azim_range=[self.setup.azimuth_min, self.setup.azimuth_max],
                elev_range=[self.setup.zenith_min, self.setup.zenith_max], v_fixed=self.setup.v_fixed, \
                allTimes=allTimes, ax=self.seis_traj_ax)
        
        # Replot 
        elif self.setup.run_mode == 'replot':
            dat = readPoints(self.setup.points_name, header=1)
            p0 = [dat[-1, 0], dat[-1, 1], dat[-1, 2], dat[-1, 3], dat[-1, 4], dat[-1, 5]]
            plotStationsAndTrajectory(station_list, p0, self.setup, sounding)
     
        elif self.setup.run_mode == 'manual':
            p0[3] *= 1000
            plotStationsAndTrajectory(station_list, p0, self.setup, sounding)

        else:
            errorMessage('Invalid mode! Use search, replot, or manual', 2)
            return None
        
        fig.set_facecolor("none")
        self.seis_tab_output.removeWidget(self.seis_three_canvas)
        self.seis_three_canvas = FigureCanvas(Figure(figsize=(15, 15)))
        self.seis_three_canvas = FigureCanvas(fig)
        self.seis_three_canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.seis_tab_output.addWidget(self.seis_three_canvas)    
        self.seis_three_canvas.draw()
        self.seis_traj_ax.mouse_init()
        SolutionGUI.update(self)

        self.seis_two_lat_canvas.scatterPlot(x=sup[0::10, 0], y=sup[0::10, 1], pen='r', update=True)
        self.seis_two_lat_canvas.setTitle('Position of Geometric Landing Point')
        self.seis_two_lat_canvas.setLabel('bottom', "x (+ East)", units='m')
        self.seis_two_lat_canvas.setLabel('left', "y (+ North)", units='m')

        self.seis_two_time_canvas.scatterPlot(x=sup[0::10, 2], y=sup[0::10, 3]*1000, pen='r', update=True)
        self.seis_two_time_canvas.setTitle('Velocty and Time of Fireball')
        self.seis_two_time_canvas.setLabel('bottom', "Time after Reference", units='s')
        self.seis_two_time_canvas.setLabel('left', "Velocity", units='m/s')

        self.seis_two_angle_canvas.scatterPlot(x=sup[0::10, 4], y=sup[0::10, 5], pen='r', update=True)
        self.seis_two_angle_canvas.setTitle('Angles of Trajectory')
        self.seis_two_angle_canvas.setLabel('bottom', "Azimuth Angle", units='deg')
        self.seis_two_angle_canvas.setLabel('left', "Zenith Angle", units='deg')

        X = [None]*len(sup[:, 0])
        Y = [None]*len(sup[:, 0])

        for i in range(len(sup[0::10, 0])):
            X[i], Y[i], _ = loc2Geo(self.setup.lat_centre, self.setup.lon_centre, 0, [sup[i*10, 0], sup[i*10, 1], 0])

        self.seis_two_plot_canvas.scatterPlot(x=X, y=Y, pen='r', update=True)
        self.seis_two_plot_canvas.setTitle('Position of Geometric Landing Point')
        self.seis_two_plot_canvas.setLabel('bottom', "Latitude", units='deg N')
        self.seis_two_plot_canvas.setLabel('left', "Longitiude", units='deg E')

        final_pos = results[0]
        t0 = results[2]
        v_est = results[3]
        azim = results[4]
        zangle = results[5]
        residuals = results[6]

        table_data = [[final_pos.lat, final_pos.lon, t0, v_est, azim, zangle]]
        toTable(self.seis_table, table_data)
        self.seis_table.setHorizontalHeaderLabels(['Latitude', 'Longitude', 'Time', 'Velocity', 'Azimuth', 'Zenith'])
        header = self.seis_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)

        toTable(self.seis_resids, residuals)
        self.seis_resids.setHorizontalHeaderLabels(['Station', 'Residual'])
        header2 = self.seis_resids.horizontalHeader()
        header2.setSectionResizeMode(QHeaderView.Stretch)

    def trajSolver(self):

        try:
            A = Position(self.setup.lat_i, self.setup.lon_i, self.setup.elev_i)
            B = Position(self.setup.lat_f, self.setup.lon_f, self.setup.elev_f)
        except:
            errorMessage("I didn't program this one yet!", 2, detail="SolutionGui (trajSolver)")

        A.pos_loc(B)
        B.pos_loc(B)

        v = np.array([B.x - A.x, B.y - A.y, B.z - A.z])

        n = (tryFloat(self.ray_height_edits.text()) - A.z)/v[2]

        P = n*v + np.array([A.x, A.y, A.z])

        pt = Position(0, 0, 0)
        pt.x = P[0]
        pt.y = P[1]
        pt.z = P[2]
        pt.pos_geo(B)

        if self.setup.debug:
            print("Solved Point:")
            print(pt)

        self.ray_lat_edits.setText(str(pt.lat))
        self.ray_lon_edits.setText(str(pt.lon))

        self.ray_pick_traj.setPoints(x=[pt.lon], y=[pt.lat], pen=(255, 0, 110))
        self.ray_canvas.addItem(self.ray_pick_traj, update=True)

        if self.ray_pick_point != [0, 0, 0]:
            self.rayTrace()

    def rayTrace(self):

        A = Position(float(self.ray_lat_edits.text()), float(self.ray_lon_edits.text()), float(self.ray_height_edits.text()))
        B = Position(self.ray_pick_point[0], self.ray_pick_point[1], self.ray_pick_point[2])

        A.pos_loc(B)
        B.pos_loc(B)

        consts = Constants()
        try:
            sounding = parseWeather(self.setup, consts)
        except:
            errorMessage('Error reading weather profile in rayTrace', 2)
            return None

        if self.setup.debug:
            print("Starting and End points of Ray Trace")
            print(A)
            print(B)

        trace_data =    [None]*self.setup.perturb_times
        t_arrival =     [None]*self.setup.perturb_times
        t_arrival_cy =  [None]*self.setup.perturb_times
        err =           [None]*self.setup.perturb_times

        if self.setup.perturb_method == 'ensemble':
            ensemble_file = self.setup.perturbation_spread_file
        else:
            ensemble_file = ''

        for ptb_n in range(self.setup.perturb_times):

            if ptb_n > 0 and self.ray_enable_perts.isChecked():
                
                if self.setup.debug:
                    print("STATUS: Perturbation {:}".format(ptb_n))

                # generate a perturbed sounding profile
                sounding_p = perturb(self.setup, sounding, self.setup.perturb_method, \
                    spread_file=self.setup.perturbation_spread_file, lat=self.setup.lat_centre, lon=self.setup.lon_centre, ensemble_file=ensemble_file, ensemble_no=ptb_n)
            else:

                # if not using perturbations on this current step, then return the original sounding profile
                sounding_p = sounding


            z_profile, _ = getWeather(np.array([A.x, A.y, A.z]), np.array([B.x, B.y, B.z]), \
                    self.setup.weather_type, A, copy.copy(sounding_p))

            trace_data[ptb_n], t_arrival[ptb_n], err[ptb_n] = slowscan(A.xyz, B.xyz, z_profile,\
                                wind=True, n_theta=self.setup.n_theta, n_phi=self.setup.n_theta, precision=self.setup.angle_precision, tol=self.setup.angle_error_tol)
            t_arrival_cy[ptb_n], _, _ = cyscan(np.array([A.x, A.y, A.z]), np.array([B.x, B.y, B.z]), z_profile,\
                                wind=True, n_theta=self.setup.n_theta, n_phi=self.setup.n_theta, precision=self.setup.angle_precision, tol=self.setup.angle_error_tol)

        if self.setup.debug:
            error_list = []
            for ii in range(self.setup.perturb_times):
                error_list.append(abs(t_arrival[ptb_n] - t_arrival_cy[ptb_n]))
            avg_error = np.mean(error_list)
            print("Mean error in computation from loss in speed: {:5.2f}".format(avg_error))
        # if (np.isnan(trace_data)) or (np.isnan(t_arrival)):

        ax = plt.axes(projection='3d')
        xline = []
        yline = []
        zline = []

        try:
            for line in trace_data[0]:
                #line[0], line[1], line[2] = loc2Geo(A.lat, A.lon, A.elev, [line[0], line[1], line[2]])

                xline.append(line[0])
                yline.append(line[1])
                zline.append(line[2])
        except:            
            if not perturb:
                errorMessage('Cannot trace rays!', 2)
                return None


        plt.style.use('dark_background')
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(5, 5)
        ax = fig.add_subplot(1, 1, 1, projection='3d')

        ax.plot3D(xline, yline, zline, 'white')
        ax.scatter(xline, yline, zline, 'black')

        x_pts = [None]*len(xline)
        y_pts = [None]*len(xline)
        for i in range(len(xline)):
            
            x_pts[i], y_pts[i], _ = loc2Geo(B.lat, B.lon, B.elev, [xline[i], yline[i], zline[i]])

        self.ray_canvas.plot(y_pts, x_pts, pen=(255, 255, 255), update=True)

        if self.ray_enable_perts.isChecked():
            for ptb_n in range(self.setup.perturb_times):
                xline = []
                yline = []
                zline = []
                if ptb_n > 0:
                    try:
                        for line in trace_data[ptb_n]:
                            #line[0], line[1], line[2] = loc2Geo(A.lat, A.lon, A.elev, [line[0], line[1], line[2]])

                            xline.append(line[0])
                            yline.append(line[1])
                            zline.append(line[2])
                    except:            
                        pass

                    ax.plot3D(xline, yline, zline, '#15ff00')
                    x_pts = [None]*len(xline)
                    y_pts = [None]*len(xline)
                    for i in range(len(xline)):
                        
                        x_pts[i], y_pts[i], _ = loc2Geo(B.lat, B.lon, B.elev, [xline[i], yline[i], zline[i]])
                    self.ray_canvas.plot(y_pts, x_pts, pen=(21, 255, 0), update=True)
                    #ax.scatter(xline, yline, zline, 'black')

        if self.ray_enable_windfield.isChecked():
            c = (z_profile[:, 1])
            mags = (z_profile[:, 2])
            dirs = (z_profile[:, 3])

            # Init the constants
            consts = Constants()

            #convert speed of sound to temp
            t = np.square(c)*consts.M_0/consts.GAMMA/consts.R

            #convert to EDN
            #dirs = np.radians(angle2NDE(np.degrees(dirs)))
            norm = Normalize()
            norm.autoscale(t)

            #convert mags and dirs to u and v
            u = mags*np.sin(dirs)
            v = mags*np.cos(dirs)

            c = t[:-1]
            c = (c.ravel() - c.min()) / c.ptp()
            c = np.concatenate((c, np.repeat(c, 2)))
            c = plt.cm.seismic(c)

            xline = []
            yline = []
            zline = []

            max_mag = np.nanmax(mags)/2500

            try:
                for line in trace_data[0]:
                    #line[0], line[1], line[2] = loc2Geo(A.lat, A.lon, A.elev, [line[0], line[1], line[2]])

                    xline.append(line[0])
                    yline.append(line[1])
                    zline.append(line[2])

                ax.quiver(xline, yline, zline, u[:-1]/max_mag, v[:-1]/max_mag, 0, color=c)
            except:
                errorMessage("Cannot trace rays!", 1)

        try:

            print('Final point error {:4.2f}: {:4.2f} m x {:4.2f} m y at {:6.2f} s'.format(err[0], xline[-1], yline[-1], t_arrival[0]))

            if setup.debug:
                F = Position(0, 0, 0)
                
                F.x = xline[-1]
                F.y = yline[-1]
                F.z = zline[-1]

                F.pos_geo(B)
            
                print('Final point:')
                print(F)

        except:
            print('Final point error {:4.2f}: {:4.2f} m x {:4.2f} m y'.format(err[0], np.nan, np.nan))


        if self.setup.perturb and self.ray_enable_perts.isChecked():
            ptb = np.nanargmin(err)

            xline = []
            yline = []
            zline = []
            try:
                for line in trace_data[ptb]:
                    #line[0], line[1], line[2] = loc2Geo(A.lat, A.lon, A.elev, [line[0], line[1], line[2]])

                    xline.append(line[0])
                    yline.append(line[1])
                    zline.append(line[2])
            except:            
                if not perturb:
                    errorMessage('Cannot trace rays!', 2)
                    return None
            print('Best perturbation {:}: error {:4.2f}: {:4.2f} m x {:4.2f} m y at {:6.2f} s'.format(ptb, err[ptb], xline[-1], yline[-1], t_arrival[ptb]))



        # try:
        #     self.ray_graphs.removeWidget(self.three_ray)
        # except:
        #     pass
        self.ray_graphs.removeWidget(self.ray_line_canvas)
        self.ray_line_canvas = FigureCanvas(Figure(figsize=(3, 3)))
        self.ray_line_canvas = FigureCanvas(fig)
        self.ray_line_canvas.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        # self.three_ray = NavigationToolbar(self.ray_line_canvas, self)
        self.ray_graphs.addWidget(self.ray_line_canvas)
        # self.ray_graphs.addWidget(self.three_ray)
        self.ray_line_canvas.draw()

        ax.mouse_init()
        SolutionGUI.update(self)

    def loadRayGraph(self):

        try:
            if not os.path.exists(self.setup.working_directory):
                os.makedirs(self.setup.working_directory)
        except FileNotFoundError:
            errorMessage("No such file or directory: '{:}'".format(self.setup.working_directory), 2)
            return None

            #Build seismic data path
        dir_path = os.path.join(self.setup.working_directory, self.setup.fireball_name)

        # Load the station and waveform files list
        data_file_path = os.path.join(dir_path, DATA_FILE)
        if os.path.isfile(data_file_path):
            
            stn_list = readStationAndWaveformsListFile(data_file_path, rm_stat=self.setup.rm_stat, debug=self.setup.debug)

        else:
            errorMessage('Station and waveform data file not found! Download the waveform files first!', 2)
            return None

        for stn in stn_list:
            text = pg.TextItem(text='{:}-{:}'.format(stn.network, stn.code),\
             border='w', color=(255, 255, 255), fill=(255, 255, 255, 100))
            
            text.setPos(stn.position.lon, stn.position.lat)
            self.ray_canvas.addItem(text)

        end_point = pg.ScatterPlotItem()
        end_point.addPoints(x=[self.setup.lon_f], y=[self.setup.lat_f], pen=(66, 232, 244), symbol='+')
        self.ray_canvas.addItem(end_point, update=True)

        x=[self.setup.lon_i, self.setup.lon_f]
        y=[self.setup.lat_i, self.setup.lat_f]
        self.ray_canvas.plot(x, y, pen=(66, 232, 244))

        SolutionGUI.update(self)

    def rayMouseClicked(self, evt):

        if self.ctrl_pressed:

            #self.ray_canvas.removeItem(self.ray_pick, update=True)

            mousePoint = self.ray_canvas.vb.mapToView(evt.pos())
            
            self.ray_pick.setPoints(x=[mousePoint.x()], y=[mousePoint.y()], pen=(255, 0, 110))
            self.ray_canvas.addItem(self.ray_pick, update=True)
            self.ray_pick_point = [mousePoint.y(), mousePoint.x(), 0]
            self.ray_pick_label.setText("Lat: {:10.4f} Lon: {:10.4f} Elev {:10.2f}".format(*self.ray_pick_point))


    def addRayTracerWidgets(self):
        ray_tab = QWidget()

        self.master_ray = QVBoxLayout()
        self.ray_graphs = QHBoxLayout()
        self.ray_control = QGridLayout()

        self.master_ray.addLayout(self.ray_graphs)
        self.master_ray.addLayout(self.ray_control)

        ray_tab.setLayout(self.master_ray)
        self.tab_widget.addTab(ray_tab, "Ray Tracer")

        self.ray_view = pg.GraphicsLayoutWidget()
        self.ray_canvas = self.ray_view.addPlot()
        self.ray_graphs.addWidget(self.ray_view)
        self.ray_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)
        self.ray_canvas.setLabel('bottom', "Longitude", units='deg E')
        self.ray_canvas.setLabel('left', "Latitude", units='deg N')

        self.ray_line_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.ray_line_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.ray_graphs.addWidget(self.ray_line_canvas)

        self.ray_label = QLabel('Starting Point')
        self.ray_control.addWidget(self.ray_label, 1, 0)
        
        self.ray_label2 = QLabel('Ending Point')
        self.ray_control.addWidget(self.ray_label2, 3, 0, 1, 1)

        self.ray_height_label, self.ray_height_edits = createLabelEditObj("Height", self.ray_control, 1, width=3)
        
        self.ray_button = QPushButton('Solve for Lat/Lon')
        self.ray_control.addWidget(self.ray_button, 7, 3, 4, 1)
        self.ray_button.clicked.connect(self.trajSolver)

        self.ray_button = QPushButton('Load Data')
        self.ray_control.addWidget(self.ray_button, 7, 1, 1, 1)
        self.ray_button.clicked.connect(self.loadRayGraph)

        self.ray_lat_label, self.ray_lat_edits = createLabelEditObj("Lat", self.ray_control, 2)
        self.ray_lon_label, self.ray_lon_edits = createLabelEditObj("Lon", self.ray_control, 2, h_shift=2)

        self.ray_pick_label = QLabel('')
        self.ray_control.addWidget(self.ray_pick_label, 3, 1, 1, 5)

        self.ray_enable_windfield = QCheckBox('Enable Wind Field')
        self.ray_control.addWidget(self.ray_enable_windfield, 1, 5)
        self.ray_enable_windfield.stateChanged.connect(self.trajSolver)

        self.ray_enable_perts = QCheckBox('Enable Perturbations')
        self.ray_control.addWidget(self.ray_enable_perts, 2, 5)
        self.ray_enable_perts.stateChanged.connect(self.trajSolver)

        self.ray_canvas.scene().sigMouseClicked.connect(self.rayMouseClicked)

    def addSeisTrajWidgets(self):

        seis_tab = QWidget()
        self.master_seis = QHBoxLayout()
        self.seis_tab_input = QVBoxLayout()
        self.seis_tab_output = QVBoxLayout()

        seis_tab.setLayout(self.master_seis)
        self.tab_widget.addTab(seis_tab, "Seismic Trajectory")

        self.master_seis.addLayout(self.seis_tab_input)
        self.master_seis.addLayout(self.seis_tab_output)

        table_group = QGridLayout()
        self.seis_tab_input.addLayout(table_group)

        tab_layout = QGridLayout()
        self.seis_tab_input.addLayout(tab_layout)

        self.seis_search = QPushButton('Search')
        tab_layout.addWidget(self.seis_search, 1, 1, 1, 4)
        self.seis_search.clicked.connect(self.seisSearch)
        
        self.seis_table = QTableWidget()
        table_group.addWidget(self.seis_table, 1, 1, 1, 4)

        self.seis_resids = QTableWidget()
        table_group.addWidget(self.seis_resids, 2, 1, 1, 4)

        self.seis_three_canvas = FigureCanvas(Figure(figsize=(5, 5)))
        self.seis_three_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.seis_tab_output.addWidget(self.seis_three_canvas)

        two_graphs = QGridLayout()
        self.seis_tab_output.addLayout(two_graphs)

        self.seis_two_lat_view = pg.GraphicsLayoutWidget()
        self.seis_two_lat_canvas = self.seis_two_lat_view.addPlot()
        two_graphs.addWidget(self.seis_two_lat_view, 1, 1)
        self.seis_two_lat_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)
        
        self.seis_two_time_view = pg.GraphicsLayoutWidget()
        self.seis_two_time_canvas = self.seis_two_time_view.addPlot()
        two_graphs.addWidget(self.seis_two_time_view, 1, 2)
        self.seis_two_time_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)
        
        self.seis_two_angle_view = pg.GraphicsLayoutWidget()
        self.seis_two_angle_canvas = self.seis_two_angle_view.addPlot()
        two_graphs.addWidget(self.seis_two_angle_view, 2, 1)
        self.seis_two_angle_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)
        
        self.seis_two_plot_view = pg.GraphicsLayoutWidget()
        self.seis_two_plot_canvas = self.seis_two_plot_view.addPlot()
        two_graphs.addWidget(self.seis_two_plot_view, 2, 2)
        self.seis_two_plot_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

    
    def addSupWidgets(self):

        sup_tab = QWidget()
        self.master_sup = QHBoxLayout()
        self.sup_tab_content = QGridLayout()
        self.sup_plots = QVBoxLayout()

        self.sup_two_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.sup_two_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.sup_plots.addWidget(self.sup_two_canvas)
        
        self.sup_three_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.sup_three_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.sup_plots.addWidget(self.sup_three_canvas)

        self.sup_search_button = QPushButton('Search')
        self.sup_tab_content.addWidget(self.sup_search_button, 5, 1, 1, 75)
        self.sup_search_button.clicked.connect(self.supSearch)

        self.sup_results_label = QLabel("Results: ")
        self.sup_tab_content.addWidget(self.sup_results_label, 3, 1, 1, 1)

        self.sup_results_table = QTableWidget(0, 0)
        self.sup_tab_content.addWidget(self.sup_results_table, 4, 1, 1, 75)

        self.master_sup.addLayout(self.sup_tab_content)
        self.master_sup.addLayout(self.sup_plots)

        sup_tab.setLayout(self.master_sup)
        self.tab_widget.addTab(sup_tab, "Supracenter PSO Search")
  
    def addSupraWidgets(self):

        supra_tab = QWidget()
        self.master_supra = QGridLayout()
        self.supra_tab_content = QVBoxLayout()
        self.plots = QVBoxLayout()

        self.two_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.two_canvas.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        self.plots.addWidget(self.two_canvas)

        self.three_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.three_canvas.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        self.plots.addWidget(self.three_canvas)

        self.results_label = QLabel("Results: ")
        self.supra_tab_content.addWidget(self.results_label)

        self.tableWidget = QTableWidget(0, 0)
        self.supra_tab_content.addWidget(self.tableWidget)

        self.search_button = QPushButton('Search')
        self.supra_tab_content.addWidget(self.search_button)
        self.search_button.clicked.connect(self.supraSearch)

        self.master_supra.addLayout(self.supra_tab_content, 1, 1, 1, 100)
        self.master_supra.addLayout(self.plots, 1, 101)

        supra_tab.setLayout(self.master_supra)
        self.tab_widget.addTab(supra_tab, "Supracenter Manual Search")

    def fatmPlot(self):

        consts = Constants()

        self.setup.lat_centre = self.fatm_lat_slide.value()*self.slider_scale
        self.setup.lon_centre = self.fatm_lon_slide.value()*self.slider_scale
        self.setup.sounding_file = self.fatm_name_edits.text()
        self.setup.weather_type = 'ecmwf'

        try:
            dataset = parseWeather(self.setup, consts)
            sounding = findECMWFSound(self.setup.lat_centre, self.setup.lon_centre, dataset)
        except:
            errorMessage('Error reading weather profile in fatmPlotProfile', 2)
            return None

        self.atm_canvas.setLabel('left', "Height", units='m')
        if self.fatm_variable_combo.currentText() == 'Temperature':
            X = sounding[:, 1]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Temperature", units='K')
        elif self.fatm_variable_combo.currentText() == 'Wind Magnitude':
            X = sounding[:, 2]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Wind Magnitude", units='m/s')
        elif self.fatm_variable_combo.currentText() == 'Wind Direction':
            X = sounding[:, 3]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Wind Direction", units='deg from N')
        elif self.fatm_variable_combo.currentText() == 'U-Component of Wind':
            dirs = angle2NDE(np.degrees(sounding[:, 3]))
            mags = sounding[:, 2]
            X = mags*np.cos(np.radians(dirs))
            Y = sounding[:, 0]
        elif self.fatm_variable_combo.currentText() == 'V-Component of Wind':
            dirs = angle2NDE(np.degrees(sounding[:, 3]))
            mags = sounding[:, 2]
            X = mags*np.sin(np.radians(dirs))
            Y = sounding[:, 0]
        else:
            errorMessage('Error reading fatmPlotProfile combo box', 2, detail=self.fatm_variable_combo.currentText())
            X = []
            Y = []

        self.fatm_canvas.clear()
        self.fatm_canvas.plot(x=X, y=Y, pen='w')
        SolutionGUI.update(self)

    def fatmFetch(self, download):

        variables = []
        self.fatm_variable_combo.clear() 

        if self.fatm_temp.isChecked():
            variables.append('temperature')
            self.fatm_variable_combo.addItem('Temperature')
        if self.fatm_u_wind.isChecked():
            variables.append('u_component_of_wind')
            self.fatm_variable_combo.addItem('U-Component of Wind')
        if self.fatm_v_wind.isChecked():
            variables.append('v_component_of_wind')
            self.fatm_variable_combo.addItem('V-Component of Wind')
        if self.fatm_u_wind.isChecked() and self.fatm_v_wind.isChecked():
            self.fatm_variable_combo.addItem('Wind Magnitude')
            self.fatm_variable_combo.addItem('Wind Direction')

        year = str(self.fatm_datetime_edits.dateTime().date().year())
        month = str(self.fatm_datetime_edits.dateTime().date().month())
        day = str(self.fatm_datetime_edits.dateTime().date().day())

        time_of = str("{:02d}".format(self.fatm_datetime_edits.dateTime().time().hour())) + \
        ':' + str("{:02d}".format(self.fatm_datetime_edits.dateTime().time().minute()))

        loc = self.fatm_name_edits.text()

        if download:
            errorMessage('Please wait until your file has downloaded...', 0)
            copernicusAPI(variables, year, month, day, time_of, loc)
            return None
        else:
            self.fatmPlot()

    def fatmValueChange(self, obj, slider):
        
        if obj == self.fatm_lat_label:
            obj.setText('Latitude: {:8.2f}'.format(slider.value()*self.slider_scale))
        elif obj == self.fatm_lon_label:
            obj.setText('Longitude: {:8.2f}'.format(slider.value()*self.slider_scale))
        else:
            errorMessage('Bad atm slider pass in fatmValueChange', 2)

        #self.atmPlotProfile(self.atm_lat_slide.value(), self.atm_lon_slide.value(), self.var_typ)
    
    def parseGeneralECMWF(self, file_name, lat, lon, time_of, variables):

        try:
        # Read the file
            dataset = Dataset(file_name, "r+", format="NETCDF4")
        except:
            errorMessage("Unable to read weather file: {:}".format(file_name), 2)
            return None

        lon_index = int(np.around((lon%360)*4))
        lat_index = int(np.around(-(lat+90)*4))
        time_index = int(time_of)

        sounding = []
        pressures = np.array(dataset.variables['level'])

        sounding.append(pressures)

        for var in variables:
            if (set([var]) < set(dataset.variables.keys())):
                sounding.append(np.array(dataset.variables[var][time_index, :, lat_index, lon_index]))
            else:
                errorMessage("WARNING: Variable {:} not found in dataset!".format(var), 1)

        sounding = np.array(sounding).transpose()

        dataset.close()
        return sounding

    def fatmPrint(self):
        
        filename = QFileDialog.getSaveFileName(self, 'Save File', '', 'Text File (*.txt)')

        setup = self.saveINI(False)

        setup.lat_centre = self.fatm_lat_slide.value()*self.slider_scale
        setup.lon_centre = self.fatm_lon_slide.value()*self.slider_scale
        setup.sounding_file = self.fatm_name_edits.text()
        atm_time = self.fatm_datetime_edits.dateTime().time().hour()
        setup.weather_type = 'ecmwf'

        variables = []

        if self.fatm_temp.isChecked():
            variables.append('t')
        if self.fatm_u_wind.isChecked():
            variables.append('u')
        if self.fatm_v_wind.isChecked():
            variables.append('v')
            
        sounding = self.parseGeneralECMWF(setup.sounding_file, setup.lat_centre, setup.lon_centre, atm_time, variables)

        if '.txt' not in filename[0]:
            filename[0] = filename[0] + '.txt'


        header = 'Pressures (hPa)'

        for element in variables:

            if element == variables[-1]:
                header = header + ', ' + element + '\n'
            else:
                header = header + ', ' + element



        with open(str(filename[0]), 'w') as f:
            
            f.write(header)

            for line in sounding:
                info = ''
                for element in line:
                    if element == line[-1]:
                        info = info + str(element) + '\n'
                    else:
                        info = info + str(element) + ','
                f.write(info)

        errorMessage('Printed out sounding data', 0, title="Print Done")

    def addFetchATMWidgets(self):
        fetch = QWidget()
        fetch_master = QHBoxLayout()
        fetch_content = QGridLayout()
        fetch_plots = QVBoxLayout()
        fetch_master.addLayout(fetch_plots)
        fetch_master.addLayout(fetch_content)
        fetch.setLayout(fetch_master)
        self.tab_widget.addTab(fetch, "Fetch Atmosphere")

        self.fatm_view = pg.GraphicsLayoutWidget()
        self.fatm_canvas = self.fatm_view.addPlot()
        fetch_plots.addWidget(self.fatm_view)
        self.fatm_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        self.fatm_variable_combo = QComboBox()
        fetch_plots.addWidget(self.fatm_variable_combo)
        self.fatm_variable_combo.currentTextChanged.connect(self.fatmPlot)

        self.fatm_name_label, self.fatm_name_edits = createLabelEditObj('Name:', fetch_content, 1)

        self.fatm_button = QPushButton('Browse')
        fetch_content.addWidget(self.fatm_button, 1, 3)
        self.fatm_button.clicked.connect(partial(self.folderSearch, self.fatm_name_edits))

        self.fatm_datetime_label = QLabel("Time of Profile:")
        self.fatm_datetime_edits = QDateTimeEdit()
        fetch_content.addWidget(self.fatm_datetime_label, 2, 1)
        fetch_content.addWidget(self.fatm_datetime_edits, 2, 2)
        self.fatm_datetime_edits.setCalendarPopup(True)

        self.fatm_variable_group = QGroupBox("Variables")
        fetch_content.addWidget(self.fatm_variable_group, 4, 1, 1, 2)

        #############################
        group_box = QGridLayout()
        self.fatm_variable_group.setLayout(group_box)

        self.fatm_temp = QCheckBox('Temperature')
        group_box.addWidget(self.fatm_temp, 1, 1)

        self.fatm_u_wind = QCheckBox('U Wind')
        group_box.addWidget(self.fatm_u_wind, 1, 2)

        self.fatm_v_wind = QCheckBox('V Wind')
        group_box.addWidget(self.fatm_v_wind, 1, 3)
        ##############################

        self.fatm_fetch = QPushButton("Download")
        fetch_content.addWidget(self.fatm_fetch, 5, 1, 1, 2)
        self.fatm_fetch.clicked.connect(partial(self.fatmFetch, True))

        self.fatm_open = QPushButton("Open")
        fetch_content.addWidget(self.fatm_open, 6, 1, 1, 2)
        self.fatm_open.clicked.connect(partial(self.fatmFetch, False))

        self.fatm_print = QPushButton("Print")
        fetch_content.addWidget(self.fatm_print, 9, 1, 1, 2)
        self.fatm_print.clicked.connect(self.fatmPrint)

        self.fatm_lat_label = QLabel("Latitude: 0")
        self.fatm_lat_slide = QSlider(Qt.Horizontal)
        fetch_content.addWidget(self.fatm_lat_label, 7, 1)
        fetch_content.addWidget(self.fatm_lat_slide, 7, 2)
        self.fatm_lat_slide.setMinimum(-90/self.slider_scale)
        self.fatm_lat_slide.setMaximum(90/self.slider_scale)
        self.fatm_lat_slide.setValue(0)
        self.fatm_lat_slide.setTickInterval(0.5)
        self.fatm_lat_slide.valueChanged.connect(partial(self.fatmValueChange, self.fatm_lat_label, self.fatm_lat_slide))

        self.fatm_lon_label = QLabel("Longitude: 0")
        self.fatm_lon_slide = QSlider(Qt.Horizontal)
        fetch_content.addWidget(self.fatm_lon_label, 8, 1)
        fetch_content.addWidget(self.fatm_lon_slide, 8, 2)
        self.fatm_lon_slide.setMinimum(-180/self.slider_scale)
        self.fatm_lon_slide.setMaximum(180/self.slider_scale)
        self.fatm_lon_slide.setValue(0)
        self.fatm_lon_slide.setTickInterval(0.5)
        self.fatm_lon_slide.valueChanged.connect(partial(self.fatmValueChange, self.fatm_lon_label, self.fatm_lon_slide))

    def addIniDockWidgets(self):

        all_vars = QGroupBox()
        self.ini_dock.setWidget(all_vars)

        dock_layout = QVBoxLayout()
        all_vars.setLayout(dock_layout)

        ini_tabs = QTabWidget()
        dock_layout.addWidget(ini_tabs)

        self.load_button = QPushButton('Load')
        dock_layout.addWidget(self.load_button)
        self.load_button.clicked.connect(self.loadINI)

        self.save_button = QPushButton('Save')
        dock_layout.addWidget(self.save_button)
        self.save_button.clicked.connect(partial(self.saveINI, True))


        tab1 = QWidget()
        tab1_content = QGridLayout()
        tab1.setLayout(tab1_content)
        ini_tabs.addTab(tab1, "General")

        self.fireball_name_label, self.fireball_name_edits = createLabelEditObj('Fireball Name:', tab1_content, 1, tool_tip='fireball_name')
        self.get_data_label, self.get_data_edits = createComboBoxObj('Get Data: ', tab1_content, 2, items=['True', 'False'], tool_tip='get_data')
        self.run_mode_label, self.run_mode_edits = createComboBoxObj('Run Mode: ', tab1_content, 3, items=['Search', 'Replot', 'Manual'], tool_tip='run_mode')
        self.debug_label, self.debug_edits = createComboBoxObj('Debug: ', tab1_content, 4, items=['True', 'False'], tool_tip='debug')
        self.working_directory_label, self.working_directory_edits, self.working_directory_buton = createFileSearchObj('Working Directory: ', tab1_content, 5, width=1, h_shift=0, tool_tip='working_directory')
        self.working_directory_buton.clicked.connect(partial(self.folderSearch, self.working_directory_edits))
        self.arrival_times_label, self.arrival_times_edits, self.arrival_times_buton = createFileSearchObj('Arrival Times:', tab1_content, 6, width=1, h_shift=0, tool_tip='arrival_times_file')
        self.arrival_times_buton.clicked.connect(partial(self.fileSearch, ['Numpy Array (*.npy)'], self.arrival_times_edits))
        self.sounding_file_label, self.sounding_file_edits, self.sounding_file_buton = createFileSearchObj('Sounding File:', tab1_content, 7, width=1, h_shift=0, tool_tip='sounding_file')
        self.sounding_file_buton.clicked.connect(partial(self.fileSearch, ['NetCDF (*.nc)', 'HDF (*.HDF)'], self.sounding_file_edits))
        self.perturbation_file_label, self.perturbation_file_edits, self.perturbation_file_buton = createFileSearchObj('Perturbation', tab1_content, 8, width=1, h_shift=0, tool_tip='perturbation_spread_file')
        self.perturbation_file_buton.clicked.connect(partial(self.fileSearch, ['NetCDF (*.nc)'], self.perturbation_file_edits))
        self.station_picks_label, self.station_picks_edits, self.station_picks_buton = createFileSearchObj('Station Picks File: ', tab1_content, 9, width=1, h_shift=0, tool_tip='station_picks_file')
        self.station_picks_buton.clicked.connect(partial(self.fileSearch, ['CSV (*.csv)', 'Text File (*.txt)'], self.station_picks_edits))
        self.points_name_label, self.points_name_edits, self.points_name_buton = createFileSearchObj('Replot Points File: ', tab1_content, 10, width=1, h_shift=0, tool_tip='points_name')
        self.points_name_buton.clicked.connect(partial(self.fileSearch, ['CSV (*.csv)'], self.points_name_edits))
        self.lat_centre_label, self.lat_centre_edits = createLabelEditObj('Latitude Center:', tab1_content, 11, tool_tip='lat_centre')
        self.lon_centre_label, self.lon_centre_edits = createLabelEditObj('Longitude Center:', tab1_content, 12, tool_tip='lon_centre')
        self.deg_radius_label, self.deg_radius_edits = createLabelEditObj('Degrees in Search Radius:', tab1_content, 13, tool_tip='deg_radius')
        self.fireball_datetime_label, self.fireball_datetime_edits = createLabelDateEditObj("Fireball Datetime", tab1_content, 14, tool_tip='fireball_datetime')
        self.v_sound_label, self.v_sound_edits = createLabelEditObj('Average Speed of Sound:', tab1_content, 15, tool_tip='v_sound')
        self.t0_label, self.t0_edits = createLabelEditObj('t0:', tab1_content, 16, width=3, tool_tip='t0')
        self.v_label, self.v_edits = createLabelEditObj('v:', tab1_content, 17, width=3, tool_tip='v')
        self.azim_label, self.azim_edits = createLabelEditObj('azim:', tab1_content, 18, width=3, tool_tip='azim')
        self.zangle_label, self.zangle_edits = createLabelEditObj('zangle:', tab1_content, 19, width=3, tool_tip='zangle')
        self.lat_i_label, self.lat_i_edits = createLabelEditObj('lat_i:', tab1_content, 20, tool_tip='lat_i')
        self.lon_i_label, self.lon_i_edits = createLabelEditObj('lon_i:', tab1_content, 21, tool_tip='lon_i')
        self.elev_i_label, self.elev_i_edits = createLabelEditObj('elev_i:', tab1_content, 22, tool_tip='elev_i')
        self.lat_f_label, self.lat_f_edits = createLabelEditObj('lat_f:', tab1_content, 23, tool_tip='lat_f')
        self.lon_f_label, self.lon_f_edits = createLabelEditObj('lon_f:', tab1_content, 24, tool_tip='lon_f')
        self.elev_f_label, self.elev_f_edits = createLabelEditObj('elev_f:', tab1_content, 25, tool_tip='elev_f')
        self.show_ballistic_waveform_label, self.show_ballistic_waveform_edits = createComboBoxObj('Show Ballistic Waveform', tab1_content, 26, items=['True', 'False'], width=2, tool_tip='show_ballistic_waveform')
        
        tab2 = QWidget()
        tab2_content = QGridLayout()
        tab2.setLayout(tab2_content)
        ini_tabs.addTab(tab2, "Sources")

        self.fragmentation_point_label = QLabel("Fragmentation Point(s)")
        self.fragmentation_point = QTableWidget(0, 4)
        tab2_content.addWidget(self.fragmentation_point, 1, 1, 1, 4)
        tab2_content.addWidget(self.fragmentation_point_label, 0, 1, 1, 4)
        self.fragmentation_point.setHorizontalHeaderLabels(['Latitude', 'Longitude', 'Elevation', 'Time'])
        header = self.fragmentation_point.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        self.fragmentation_point_label.setToolTip(toolTime('fragmentation_point'))
        self.fragmentation_point_add = QPushButton("+")
        tab2_content.addWidget(self.fragmentation_point_add, 2, 3, 1, 2)
        self.fragmentation_point_add.clicked.connect(partial(self.changeRows, self.fragmentation_point, 1))
        self.fragmentation_point_add.setToolTip("Add row")
        self.fragmentation_point_min = QPushButton("-")
        tab2_content.addWidget(self.fragmentation_point_min, 2, 1, 1, 2)
        self.fragmentation_point_min.clicked.connect(partial(self.changeRows, self.fragmentation_point, -1))
        self.fragmentation_point_min.setToolTip("Remove row")
        self.show_fragmentation_waveform_label, self.show_fragmentation_waveform_edits = createComboBoxObj("Show Fragmentation Waveform: ", tab2_content, 3, width=2, items=['True', 'False'], tool_tip='show_fragmentation_waveform')
        self.manual_label = QLabel("Manual Fragmentation Search:")
        tab2_content.addWidget(self.manual_label, 4, 1, 1, 4)
        self.manual_label.setToolTip("manual_fragmentation_search")
        self.lat_frag_label = QLabel("Latitude:")
        self.lat_frag_edits = QLineEdit("")
        tab2_content.addWidget(self.lat_frag_label, 5, 1)
        tab2_content.addWidget(self.lat_frag_edits, 6, 1)
        self.lon_frag_label = QLabel("Longitude:")
        self.lon_frag_edits = QLineEdit("")
        tab2_content.addWidget(self.lon_frag_label, 5, 2)
        tab2_content.addWidget(self.lon_frag_edits, 6, 2)
        self.elev_frag_label = QLabel("Elevation:")
        self.elev_frag_edits = QLineEdit("")
        tab2_content.addWidget(self.elev_frag_label, 5, 3)
        tab2_content.addWidget(self.elev_frag_edits, 6, 3)
        self.time_frag_label = QLabel("Time:")
        self.time_frag_edits = QLineEdit("")
        tab2_content.addWidget(self.time_frag_label, 5, 4)
        tab2_content.addWidget(self.time_frag_edits, 6, 4)
        self.v_fixed_label, self.v_fixed_edits = createLabelEditObj('v_fixed:', tab2_content, 7, width=3, tool_tip='v_fixed')
        self.restricted_time_check = QCheckBox("Enable Restricted Time: ")
        tab2_content.addWidget(self.restricted_time_check, 9, 4, 1, 1)
        self.restricted_time_label, self.restricted_time_edits = createLabelDateEditObj("Restricted Time: ", tab2_content, 9, width=2, tool_tip='restricted_time')
        self.azimuth_min_label, self.azimuth_min_edits = createLabelEditObj('azimuth_min:', tab2_content, 10, tool_tip='azimuth_min')
        self.azimuth_max_label, self.azimuth_max_edits = createLabelEditObj('azimuth_max:', tab2_content, 10, h_shift=2, tool_tip='azimuth_max')
        self.zangle_min_label, self.zangle_min_edits = createLabelEditObj('zangle_min:', tab2_content, 11, tool_tip='zenith_min')
        self.zangle_max_label, self.zangle_max_edits = createLabelEditObj('zangle_max:', tab2_content, 11, h_shift=2, tool_tip='zenith_max')
        self.lat_min_label, self.lat_min_edits = createLabelEditObj('lat_min', tab2_content, 12, tool_tip='x_min')
        self.lat_max_label, self.lat_max_edits = createLabelEditObj('lat_max', tab2_content, 12, h_shift=2, tool_tip='x_max')
        self.lon_min_label, self.lon_min_edits = createLabelEditObj('lon_min', tab2_content, 13, tool_tip='y_min')
        self.lon_max_label, self.lon_max_edits = createLabelEditObj('lon_max', tab2_content, 13, h_shift=2, tool_tip='y_max')
        self.elev_min_label, self.elev_min_edits = createLabelEditObj('elev_min:', tab2_content, 14, tool_tip='z_min')
        self.elev_max_label, self.elev_max_edits = createLabelEditObj('elev_max:', tab2_content, 14, h_shift=2, tool_tip='z_max')
        self.t_min_label, self.t_min_edits = createLabelEditObj('t_min:', tab2_content, 15, tool_tip='t_min')
        self.t_max_label, self.t_max_edits = createLabelEditObj('t_max:', tab2_content, 15, h_shift=2, tool_tip='t_max')
        self.v_min_label, self.v_min_edits = createLabelEditObj('v_min:', tab2_content, 16, tool_tip='v_min')
        self.v_max_label, self.v_max_edits = createLabelEditObj('v_max:', tab2_content, 16, h_shift=2, tool_tip='v_max')
        self.weight_distance_min_label, self.weight_distance_min_edits = createLabelEditObj('weight_distance_min', tab2_content, 17, tool_tip='weight_distance_min')
        self.weight_distance_max_label, self.weight_distance_max_edits = createLabelEditObj('weight_distance_max', tab2_content, 17, h_shift=2, tool_tip='weight_distance_max')

        tab3 = QWidget()
        tab3_content = QGridLayout()
        tab3.setLayout(tab3_content)
        ini_tabs.addTab(tab3, "Atmosphere")

        self.enable_winds_label, self.enable_winds_edits = createComboBoxObj('Enable Winds: ', tab3_content, 1, items=['True', 'False'], tool_tip='enable_winds')
        self.weather_type_label, self.weather_type_edits = createComboBoxObj('Weather Type: ', tab3_content, 2, items=['none', 'ecmwf', 'ukmo', 'merra', 'custom'], tool_tip='weather_type')
        self.perturb_times_label, self.perturb_times_edits = createLabelEditObj('Perturbation Times', tab3_content, 3, tool_tip='perturb_times')
        self.frag_no_label, self.frag_no_edits = createLabelEditObj('Fragmentation Number', tab3_content, 4, tool_tip='fragno')
        self.perturb_label, self.perturb_edits = createComboBoxObj('Perturb: ', tab3_content, 5, items=['True', 'False'], tool_tip='perturb')
        self.perturb_method_label, self.perturb_method_edits = createComboBoxObj('Perturb Method', tab3_content, 6, items=['none', 'bmp', 'sppt', 'temporal', 'spread', 'spread_r', 'ensemble'], tool_tip='perturb_method')

        tab4 = QWidget()
        tab4_content = QGridLayout()
        tab4.setLayout(tab4_content)
        ini_tabs.addTab(tab4, "Tweaks")

        self.n_theta_label, self.n_theta_edits = createLabelEditObj('Theta Resolution', tab4_content, 1, tool_tip='n_theta')
        self.n_phi_label, self.n_phi_edits = createLabelEditObj('Phi Resolution', tab4_content, 2, tool_tip='n_phi')
        self.angle_precision_label, self.angle_precision_edits = createLabelEditObj('Angle Precision', tab4_content, 3, tool_tip='angle_precision')
        self.angle_error_tol_label, self.angle_error_tol_edits = createLabelEditObj('Angle Error Tolerance', tab4_content, 4, tool_tip='angle_error_tol')
        self.maxiter_label, self.maxiter_edits = createLabelEditObj('Max Iterations: ', tab4_content, 5, tool_tip='maxiter')
        self.swarmsize_label, self.swarmsize_edits = createLabelEditObj('Swarm Size: ', tab4_content, 6, tool_tip='swarmsize')
        self.run_times_label, self.run_times_edits = createLabelEditObj('Run Times:', tab4_content, 7, tool_tip='run_times')
        self.minfunc_label, self.minfunc_edits = createLabelEditObj('minfunc:', tab4_content, 8, tool_tip='minfunc')
        self.minstep_label, self.minstep_edits = createLabelEditObj('minstep:', tab4_content, 9, tool_tip='minstep')
        self.phip_label, self.phip_edits = createLabelEditObj('phip:', tab4_content, 10, tool_tip='phip')
        self.phig_label, self.phig_edits = createLabelEditObj('phig:', tab4_content, 11, tool_tip='phig')
        self.omega_label, self.omega_edits = createLabelEditObj('omega:', tab4_content, 12, tool_tip='omega')
        self.pso_debug_label, self.pso_debug_edits = createComboBoxObj("PSO Debug: ", tab4_content, 13, items=['True', 'False'], tool_tip='pso_debug')
        self.contour_res_label, self.contour_res_edits = createLabelEditObj('Contour Resolution:', tab4_content, 14, tool_tip='contour_res')
        self.high_f_label, self.high_f_edits = createLabelEditObj('Highlight Fragmentation:', tab4_content, 15, tool_tip='high_f')
        self.high_b_label, self.high_b_edits = createLabelEditObj('Highlight Ballistic:', tab4_content, 16, tool_tip='high_b')
        self.rm_stat_label, self.rm_stat_edits = createLabelEditObj('Remove Stations:', tab4_content, 17, tool_tip='rm_stat')


        tab5 = QWidget()
        tab5_content = QGridLayout()
        tab5.setLayout(tab5_content)
        ini_tabs.addTab(tab5, "Misc")

        self.extra_label = QLabel("Manually Added Stations")
        self.extra_point = QTableWidget(0, 8)
        tab5_content.addWidget(self.extra_point, 1, 1, 1, 2)
        tab5_content.addWidget(self.extra_label, 0, 1, 1, 2)
        self.extra_point.setHorizontalHeaderLabels(\
            ['Station Network', 'Station Code', 'Latitude', 'Longitude', \
            'Elevation', 'Station Display Name', 'Channel', 'File Name'])
        header = self.extra_point.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        self.extra_label.setToolTip(toolTime('stations'))
        self.extra_point_add = QPushButton("+")
        tab5_content.addWidget(self.extra_point_add, 2, 2, 1, 1)
        self.extra_point_add.clicked.connect(partial(self.changeRows, self.extra_point, 1))
        self.extra_point_min = QPushButton("-")
        tab5_content.addWidget(self.extra_point_min, 2, 1, 1, 1)
        self.extra_point_min.clicked.connect(partial(self.changeRows, self.extra_point, -1))

    def atmPlotProfile(self, lat, lon, var_typ='t', perturb='none'):
        
        consts = Constants()

        if self.setup.weather_type == 'none':
            errorMessage('Weather type is set to "none", no weather can be displayed', 1)
            return None

        try:
            dataset = parseWeather(self.setup, consts)
            sounding = findECMWFSound(lat, lon, dataset)
        except:
            errorMessage('Error reading weather profile in atmPlotProfile', 2)
            return None

        self.var_typ = var_typ
        self.atm_canvas.setLabel('left', "Height", units='m')
        if self.var_typ == 't':
            X = sounding[:, 1]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Temperature", units='K')
        elif self.var_typ == 'm':
            X = sounding[:, 2]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Wind Magnitude", units='m/s')
        elif self.var_typ == 'd':
            X = sounding[:, 3]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Wind Direction", units='deg from N')
        else:
            errorMessage('Error reading var_typ in atmPlotProfile', 2)
            return None

        self.atm_canvas.clear()
        self.atm_canvas.plot(x=X, y=Y, pen='w')
        SolutionGUI.update(self)

        if self.setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(self.setup, consts, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(self.setup, consts, time=-1)

        else:
            sounding_u = []
            sounding_l = []

        if self.setup.perturb_method == 'ensemble':
            ensemble_file = self.setup.perturbation_spread_file
        else:
            ensemble_file = ''

        if self.setup.perturb_method != 'none':
            for ptb_n in range(self.setup.perturb_times):

                if ptb_n > 0:
                    
                    if self.setup.debug:
                        print("STATUS: Perturbation {:}".format(ptb_n))

                    # generate a perturbed sounding profile
                    sounding_p = perturbation_method(self.setup, dataset, setup.perturb_method, \
                        sounding_u=sounding_u, sounding_l=sounding_l, \
                        spread_file=self.setup.perturbation_spread_file, lat=self.setup.lat_centre, lon=self.setup.lon_centre, ensemble_file=ensemble_file, ensemble_no=ptb_n)
                    sounding_p = findECMWFSound(lat, lon, sounding_p)
                    
                    if self.var_typ == 't':
                        X = sounding_p[:, 1]
                        Y = sounding_p[:, 0]
                    elif self.var_typ == 'm':
                        X = sounding_p[:, 2]
                        Y = sounding_p[:, 0]
                    elif self.var_typ == 'd':
                        X = sounding_p[:, 3]
                        Y = sounding_p[:, 0]
                    else:
                        print('error, atmPlotProfile')

                    self.atm_canvas.plot(x=X, y=Y, pen='g')
                    SolutionGUI.update(self)

    def atmValueChange(self, obj, slider):
        
        if obj == self.atm_lat_label:
            obj.setText('Latitude: {:8.2f}'.format(slider.value()*self.slider_scale))
        elif obj == self.atm_lon_label:
            obj.setText('Longitude: {:8.2f}'.format(slider.value()*self.slider_scale))
        else:
            errorMessage(self, 'Bad atm slider pass in atmValueChange', 2)

        self.atmPlotProfile(self.atm_lat_slide.value()*self.slider_scale, self.atm_lon_slide.value()*self.slider_scale, self.var_typ)
    
    def addProfileWidgets(self):
        profile_tab = QWidget()
        profile_master = QHBoxLayout()
        profile_tab_content = QGridLayout()
        profile_tab_content_graph = QVBoxLayout()
        profile_tab.setLayout(profile_master)
        profile_master.addLayout(profile_tab_content_graph)
        profile_master.addLayout(profile_tab_content)

        self.atm_lat_label = QLabel("Latitude: 0")
        self.atm_lat_slide = QSlider(Qt.Horizontal)
        profile_tab_content.addWidget(self.atm_lat_label, 6, 1)
        profile_tab_content.addWidget(self.atm_lat_slide, 6, 2, 1, 3)
        self.atm_lat_slide.setMinimum(-90/self.slider_scale)
        self.atm_lat_slide.setMaximum(90/self.slider_scale)
        self.atm_lat_slide.setValue(0)
        self.atm_lat_slide.setTickInterval(0.5)
        self.atm_lat_slide.valueChanged.connect(partial(self.atmValueChange, self.atm_lat_label, self.atm_lat_slide))

        self.atm_lon_label = QLabel("Longitude: 0")
        self.atm_lon_slide = QSlider(Qt.Horizontal)
        profile_tab_content.addWidget(self.atm_lon_label, 7, 1)
        profile_tab_content.addWidget(self.atm_lon_slide, 7, 2, 1, 3)
        self.atm_lon_slide.setMinimum(-180/self.slider_scale)
        self.atm_lon_slide.setMaximum(180/self.slider_scale)
        self.atm_lon_slide.setValue(0)
        self.atm_lon_slide.setTickInterval(0.5)
        self.atm_lon_slide.valueChanged.connect(partial(self.atmValueChange, self.atm_lon_label, self.atm_lon_slide))

        self.atm_view = pg.GraphicsLayoutWidget()
        self.atm_canvas = self.atm_view.addPlot()

        profile_tab_content_graph.addWidget(self.atm_view)
        self.atm_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        self.atm_T_button = QPushButton('Temperature')
        profile_tab_content_graph.addWidget(self.atm_T_button)
        self.atm_T_button.clicked.connect(partial(self.atmPlotProfile, self.atm_lat_slide.value()*self.slider_scale, self.atm_lon_slide.value()*self.slider_scale, var_typ='t'))

        self.atm_mag_button = QPushButton('Wind Magnitude')
        profile_tab_content_graph.addWidget(self.atm_mag_button)
        self.atm_mag_button.clicked.connect(partial(self.atmPlotProfile, self.atm_lat_slide.value()*self.slider_scale, self.atm_lon_slide.value()*self.slider_scale, var_typ='m'))

        self.atm_dir_button = QPushButton('Wind Direction')
        profile_tab_content_graph.addWidget(self.atm_dir_button)
        self.atm_dir_button.clicked.connect(partial(self.atmPlotProfile, self.atm_lat_slide.value()*self.slider_scale, self.atm_lon_slide.value()*self.slider_scale, var_typ='d'))

        self.tab_widget.addTab(profile_tab, "Atmospheric Profile")

    def makeStationObj(self, lst):

        new_lst = []

        for line in lst:
            pos = Position(line[2], line[3], line[4])
            stn = Station(line[0], line[1], pos, line[5], line[6], line[7])
            new_lst.append(stn)

        return new_lst
    
    def makePicks(self):


        if not os.path.exists(self.setup.working_directory):
            os.makedirs(self.setup.working_directory)

            #Build seismic data path
        dir_path = os.path.join(self.setup.working_directory, self.setup.fireball_name)

        # Load the station and waveform files list
        data_file_path = os.path.join(dir_path, DATA_FILE)
        if os.path.isfile(data_file_path):
            
            stn_list = readStationAndWaveformsListFile(data_file_path, rm_stat=self.setup.rm_stat)

        else:
            errorMessage('Station and waveform data file not found! Download the waveform files first!', 2)
            sys.exit()

        # if self.setup.stations is not None:
        #     stn_list = stn_list + self.makeStationObj(self.setup.stations)

        # Init the constants
        self.setup.search_area = [0, 0, 0, 0]
        self.setup.search_area[0] = self.setup.lat_centre - self.setup.deg_radius 
        self.setup.search_area[1] = self.setup.lat_centre + self.setup.deg_radius
        self.setup.search_area[2] = self.setup.lon_centre - self.setup.deg_radius
        self.setup.search_area[3] = self.setup.lon_centre + self.setup.deg_radius

        sounding = parseWeather(self.setup)
        arrTimes = np.array([])

        if len(stn_list) == 0:
            errorMessage('No Stations to load', 2)
            return None 

        try:
            #turn coordinates into position objects
            self.setup.traj_i = Position(self.setup.lat_i, self.setup.lon_i, self.setup.elev_i)
            self.setup.traj_f = Position(self.setup.lat_f, self.setup.lon_f, self.setup.elev_f)
        except:
            self.setup.traj_i = Position(0, 0, 0)
            self.setup.traj_f = Position(0, 0, 0)
            errorMessage("Warning: Unable to build trajectory points", 1)

        self.waveformPicker(dir_path, self.setup, sounding, stn_list, stn_list=stn_list)
    
    # def trajCalc(setup):
    #     """ Creates trajectory between point A and the ground (B) based off of the initial position and the angle of travel

    #     Arguments:
    #     setup: [Object] ini file parameters

    #     Returns: 
    #     A [list] lat/lon/elev of the tail of the trajectory
    #     B [list] lat/lon/elev of the head of the trajectory
    #     """
    #     print("Trajectory Calculation...")

    #     ref = Position(setup.lat_centre, setup.lon_centre, 0)

    #     # Tail of the trajectory
    #     A = geo2Loc(ref.lat, ref.lon, ref.elev, setup.lat_i, setup.lon_i, setup.elev_i*1000)

    #     # convert angles to radians
    #     ze = np.radians(setup.zangle)
    #     az = np.radians(setup.azim)

    #     # Create trajectory vector
    #     traj = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

    #     # How far along the trajectory until it reaches the ground
    #     n = -A[2]/traj[2]

    #     # B is the intersection between the trajectory vector and the ground
    #     B = A + n*traj

    #     # Convert back to geo coordinates
    #     B = np.array(loc2Geo(ref.lat, ref.lon, ref.elev, B))
    #     A = np.array(loc2Geo(ref.lat, ref.lon, ref.elev, A))

    #     print("Created Trajectory between A and B:")
    #     print("     A = {:10.4f}N {:10.4f}E {:10.2f}m".format(A[0], A[1], A[2]))
    #     print("     B = {:10.4f}N {:10.4f}E {:10.2f}m".format(B[0], B[1], B[2]))

    #     return A, B


    def calcAllTimes(self, stn_list, setup, sounding):
        """ 
        Method to calculate all arrival times and place them into an array, so that they are not recalculated every
        time a new waveform is opened
                
        Arguments:
        data_list: [ndarray] array containing all station names and locations
        setup: [object] ini file parameters
        sounding: [ndarray] atmospheric profile

        Returns:
        allTimes: [ndarray] a ndarray that stores the arrival times for all stations over every perturbation
        [perturbation, station, 0 - ballistic/ 1 - fragmentation, frag number (0 for ballistic)]
        """

        #All perturbation happens here
        allTimes = [0]*setup.perturb_times

        # Ballistic Prediction
        ref_pos = Position(setup.lat_centre, setup.lon_centre, 0)

        no_of_frags = len(setup.fragmentation_point)

        # array of frags and ballistic arrivals have to be the same size. So, minimum can be 1
        if no_of_frags == 0:
            no_of_frags = 1

        # Initialize variables
        b_time = 0

        consts = Constants()

        # For temporal perturbations, fetch the soudning data for the hour before and after the event
        if setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(setup, consts, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(setup, consts, time=-1)

        else:
            sounding_u = []
            sounding_l = []

        if setup.perturb_method == 'ensemble':
            ensemble_file = setup.perturbation_spread_file
        else:
            ensemble_file = ''

        if setup.show_fragmentation_waveform and setup.show_ballistic_waveform:
            d_time = 2*(setup.perturb_times*len(stn_list)*no_of_frags)
        else:
            d_time = (setup.perturb_times*len(stn_list)*no_of_frags)
        count = 0

        #number of perturbations
        for ptb_n in range(setup.perturb_times):

            if ptb_n > 0:
                
                if setup.debug:
                    print("STATUS: Perturbation {:}".format(ptb_n))

                # generate a perturbed sounding profile
                sounding_p = perturb(setup, sounding, setup.perturb_method, \
                    sounding_u=sounding_u, sounding_l=sounding_l, \
                    spread_file=setup.perturbation_spread_file, lat=setup.lat_centre, lon=setup.lon_centre, ensemble_file=ensemble_file, ensemble_no=ptb_n)
            else:

                # if not using perturbations on this current step, then return the original sounding profile
                sounding_p = sounding

            # Initialize station times array
            stnTimes = [0]*len(stn_list)

            #number of stations
            for n, stn in enumerate(stn_list):

                # For ballistic arrivals
                if setup.show_ballistic_waveform:
                    bTimes = [0]*no_of_frags
                    for i in range(no_of_frags):
                        count += 1
                        loadingBar('Calculating all times:', count, d_time)
                        # sys.stdout.write("\rCalculating all times: {:5.2f} % ".format(count/d_time * 100))
                        # sys.stdout.flush()
                        #need filler values to make this a numpy array with fragmentation
                        if i == 0:
                            
                            stn.position.pos_loc(ref_pos)
                            setup.traj_f.pos_loc(ref_pos)

                            # Time to travel from trajectory to station
                            b_time = timeOfArrival(stn.position.xyz, self.setup.traj_f.x/1000, setup.traj_f.y/1000, self.setup.t0, 1000*self.setup.v, \
                                                        self.setup.azimuth.rad, self.setup.zenith.rad, self.setup, sounding=sounding_p, travel=False, fast=False, ref_loc=ref_pos)# + setup.t 

                            bTimes[i] = b_time
                        else:
                            bTimes[i] = np.nan

                else:
                    bTimes = [np.nan]*no_of_frags

                # Fragmentation Prediction
                f_time = np.array([0]*no_of_frags)

                # If manual fragmentation search is on
                if setup.show_fragmentation_waveform:
                    fTimes = [0]*no_of_frags
                    for i, line in enumerate(setup.fragmentation_point):
                        count += 1
                        loadingBar('Calculating all times:', count, d_time)

                        # location of supracenter
                        supra = Position(float(line[0]), float(line[1]), float(line[2]))
                        
                        # convert to local coordinates based off of the ref_pos
                        supra.pos_loc(ref_pos)

                        # convert station coordinates to local coordinates based on the ref_pos
                        stn.position.pos_loc(ref_pos)

                        # Cut down atmospheric profile to the correct heights, and interp
                        zProfile, _ = getWeather(np.array([supra.lat, supra.lon, supra.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), setup.weather_type, \
                                [ref_pos.lat, ref_pos.lon, ref_pos.elev], copy.copy(sounding_p), convert=False)

                        # Travel time of the fragmentation wave
                        f_time, _, _ = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True, \
                            n_theta=setup.n_theta, n_phi=setup.n_phi, precision=setup.angle_precision, tol=setup.angle_error_tol)

                        fTimes[i] = f_time + line[3]

                else:

                    # Repack all arrays into allTimes array
                    fTimes = [np.nan]*no_of_frags

                stnTimes[n] = ([np.array(bTimes), np.array(fTimes)])

            allTimes[ptb_n] = np.array(stnTimes)

        allTimes = np.array(allTimes)

        # Save as .npy file to be reused in SeismicTrajectory and Supracenter
        np.save(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times'), allTimes)
        errorMessage("All Times File saved as {:}".format(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times.npy')), 0)

        return allTimes


    def waveformPicker(self, dir_path, setup, sounding, data_list, waveform_window=600, \
        difference_filter_all=False, stn_list=[]):
        """

        Arguments:
            data_list: [list]

        Keyword arguments:
            waveform_window: [int] Number of seconds for the waveform window.
            difference_filter_all: [bool] If True, the Kalenda et al. (2014) difference filter will be applied
                on the data plotted in the overview plot of all waveforms.
        """
        self.dir_path = dir_path

        self.v_sound = setup.v_sound
        self.t0 = setup.t0

        self.stn_list = stn_list

        # Filter out all stations for which the mseed file does not exist
        filtered_stn_list = []

        names = []
        lats = []
        lons = []

        for stn in self.stn_list:
            
            mseed_file = stn.file_name
            mseed_file_path = os.path.join(self.dir_path, mseed_file)

            if os.path.isfile(mseed_file_path):
                filtered_stn_list.append(stn)

            else:
                if self.setup.debug:
                    print('mseed file does not exist:', mseed_file_path)

        self.stn_list = filtered_stn_list


        self.lat_centre = setup.lat_centre
        self.lon_centre = setup.lon_centre

        self.waveform_window = waveform_window


        self.current_station = 0
        self.current_waveform_raw = None
        self.current_waveform_delta = None
        self.current_waveform_processed = None

        # List of picks
        self.pick_list = []

        self.pick_group = 0

        # Define a list of colors for groups
        self.pick_group_colors = ['r', 'g', 'm', 'w', 'y']

        # Current station map handle
        self.current_station_scat = None

        # Station waveform marker handle
        self.current_station_all_markers = None

        # Picks on all waveform plot handle
        self.all_waves_picks_handle = None

        # Handle for pick text
        self.pick_text_handle = None
        self.pick_markers_handles = []

        # handle for pick marker on the waveform
        self.pick_waveform_handle = None


        # Default bandpass values
        self.bandpass_low_default = 2.0
        self.bandpass_high_default = 8.0

        # Flag indicating whether CTRL is pressed or not
        self.ctrl_pressed = False
        self.shift_pressed = False

        ### Sort stations by distance from source ###

        # Calculate distances of station from source
        self.source_dists = []


        for stn in self.stn_list:

            stat_name, stat_lat, stat_lon = stn.code, stn.position.lat, stn.position.lon

            names.append(stat_name)
            lats.append(stat_lat)
            lons.append(stat_lon)

            # Calculate the distance in kilometers
            dist = greatCircleDistance(np.radians(self.setup.lat_centre), np.radians(self.setup.lon_centre), \
                np.radians(stat_lat), np.radians(stat_lon))

            self.source_dists.append(dist)

        # Get sorted arguments
        dist_sorted_args = np.argsort(self.source_dists)

        # Sort the stations by distance
        self.stn_list = [self.stn_list[i] for i in dist_sorted_args]
        self.source_dists = [self.source_dists[i] for i in dist_sorted_args]
        #############################################
        
        # Init the plot framework
        self.initPlot(setup, sounding)

        # Extract the list of station locations
        # self.lat_list = [stn.position.lat_r for stn in stn_list]
        # self.lon_list = [stn.position.lon_r for stn in stn_list]

        # plt.style.use('dark_background')
        # fig = plt.figure(figsize=plt.figaspect(0.5))
        # fig.set_size_inches(8, 5)
        # self.map_ax = fig.add_subplot(1, 1, 1)

        # Init ground map
        #self.m = GroundMap(self.lat_list, self.lon_list, ax=self.map_ax, color_scheme='light')

        # Extract coordinates of the reference station
        ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)
        self.make_picks_map_graph_canvas.setLabel('bottom', "Longitude", units='deg E')
        self.make_picks_map_graph_canvas.setLabel('left', "Latitude", units='deg N')

        for ii, stn in enumerate(self.stn_list):

            self.station_marker[ii].setPoints(x=[stn.position.lon], y=[stn.position.lat], pen=(255, 255, 255), brush=(255, 255, 255), symbol='t')
            self.make_picks_map_graph_canvas.addItem(self.station_marker[ii], update=True)
            # # Plot stations
            # if stn.code in setup.high_f:
            #     self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='g', s=2)
            # elif stn.code in setup.high_b:
            #     self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='b', s=2)
            # else:
            #     self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='k', s=2)
        
        # Manual Supracenter search
        if self.setup.show_fragmentation_waveform:
            
            # Fragmentation plot
            for i, line in enumerate(self.setup.fragmentation_point):
                self.make_picks_map_graph_canvas.scatterPlot(x=[float(line[1])], y=[float(line[0])],\
                    pen=(0 + i*255/len(self.setup.fragmentation_point), 255 - i*255/len(self.setup.fragmentation_point), 0), symbol='+')



        # Plot source location
        self.make_picks_map_graph_canvas.scatterPlot(x=[setup.lon_centre], y=[setup.lat_centre], symbol='+', pen=(255, 255, 0))

        # Manual trajectory search
        if self.setup.show_ballistic_waveform:

            # Plot the trajectory with the bottom point known
            self.make_picks_map_graph_canvas.plot([self.setup.traj_i.lon, self.setup.traj_f.lon], [self.setup.traj_i.lat, self.setup.traj_f.lat],\
                             pen=(0, 0, 255))
            # Plot intersection with the ground
            self.make_picks_map_graph_canvas.scatterPlot(x=[self.setup.traj_f.lon], y=[self.setup.traj_f.lat], symbol='+', pen=(0, 0, 255))


        if self.setup.arrival_times_file != '':
            try:
                self.arrTimes = np.load(self.setup.arrival_times_file)
            except:
                errorMessage("WARNING: Unable to load allTimes_file {:} . Please check that file exists".format(self.setup.arrival_times_file), 1)
                self.arrTimes = self.calcAllTimes(self.stn_list, setup, sounding)
        else:  
            # Calculate all arrival times
            self.arrTimes = self.calcAllTimes(self.stn_list, setup, sounding)
    
        SolutionGUI.update(self)

        self.updatePlot(self.setup)

    def makeValueChange(self, obj, slider):
        
        if obj == self.low_bandpass_label:
            obj.setText('Low: {:8.2f} Hz'.format(slider.value()*self.bandpass_scale))
        elif obj == self.high_bandpass_label:
            obj.setText('High: {:8.2f} Hz'.format(slider.value()*self.bandpass_scale))
        else:
            errorMessage('Bad atm slider pass in makeValueChange', 2)

        self.updatePlot()
    
    def chooseFilter(self, obj):

        if obj.currentText() == 'Bandpass Filter':
            self.filterBandpass()
        elif obj.currentText() == 'Spectrogram of Raw Data':
            self.showSpectrogram()
        elif obj.currentText() == 'Difference Filter':
            self.filterConvolution()
        else:
            pass

    def navStats(self):

        self.current_station = self.make_picks_station_choice.currentIndex()
        self.updatePlot()


    def initPlot(self, setup, sounding):
        """ Initializes the plot framework. """
        
               ### Init the basic grid ###


        # Register a mouse press event on the waveform axis
        # plt.gca().figure.canvas.mpl_connect('button_press_event', self.onWaveMousePress)


        # Register window resize
        #plt.gca().figure.canvas.mpl_connect('resize_event', self.onResize)

        self.make_picks_station_choice.clear()

        self.pick_list = []

        self.prev_stat.clicked.connect(self.decrementStation)
        self.next_stat.clicked.connect(self.incrementStation)

        self.filter_combo_box.addItem('Raw Data')
        self.filter_combo_box.addItem('Bandpass Filter')
        self.filter_combo_box.addItem('Spectrogram of Raw Data')
        self.filter_combo_box.addItem('Difference Filter')

        self.filter_combo_box.currentTextChanged.connect(partial(self.chooseFilter, self.filter_combo_box))

        self.export_to_csv.clicked.connect(self.exportCSV)
        self.export_to_all_times.clicked.connect(self.exportToAllTimes)

        self.station_marker = []
        self.station_waveform = []
        for stn in self.stn_list:
            self.make_picks_station_choice.addItem("{:}-{:}".format(stn.network, stn.code))
            self.station_marker.append(pg.ScatterPlotItem())
            self.station_waveform.append(pg.PlotCurveItem())

        self.make_picks_station_choice.currentTextChanged.connect(self.navStats)

        plt.style.use('dark_background')
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(8, 5)
        self.station_ax = fig.add_subplot(1, 1, 1)
        
        self.make_picks_waveform_canvas.scene().sigMouseClicked.connect(self.mouseClicked)
        pg.QtGui.QApplication.processEvents()
        # Plot all waveforms
        ######################################################################################3

        max_wave_value = 0
        min_wave_value = np.inf
        min_time = np.inf
        max_time = 0

        lats = []
        lons = []
        for i in range(len(self.stn_list)):
            lats.append(self.stn_list[i].position.lat)
            lons.append(self.stn_list[i].position.lon)
        
        self.ballistic_idx = []
        self.fragmentation_idx = []


        # Go though all stations and waveforms
        bad_stats = []

        for idx, stn in enumerate(self.stn_list):

            sys.stdout.write('\rPlotting: {:} {:}              '.format(stn.network, stn.code))
            sys.stdout.flush()
            time.sleep(0.001)

            mseed_file_path = os.path.join(self.dir_path, stn.file_name)

            try:
                
                # Read the miniSEED file
                if os.path.isfile(mseed_file_path):

                    mseed = obspy.read(mseed_file_path)

                else:
                    bad_stats.append(idx)
                    if self.setup.debug:
                        print('File {:s} does not exist!'.format(mseed_file_path))
                    continue


            except TypeError as e:
                bad_stats.append(idx)
                if self.setup.debug:
                    print('Opening file {:} failed with error: {:}'.format(mseed_file_path, e))
                continue

            # Find channel with BHZ, HHZ, or BDF

            for i in range(len(mseed)):
                if mseed[i].stats.channel == 'BDF':
                    stn.channel = 'BDF'
                    stream = i

            for i in range(len(mseed)):
                if mseed[i].stats.channel == 'BHZ':
                    stn.channel = 'BHZ'
                    stream = i

            for i in range(len(mseed)):
                if mseed[i].stats.channel == 'HHZ':
                    stn.channel = 'HHZ'
                    stream = i

            for i in range(len(mseed)):
                if mseed[i].stats.channel == 'EHZ':
                    stn.channel = 'EHZ'
                    stream = i

            for i in range(len(mseed)):
                if mseed[i].stats.channel == 'SHZ':
                    stn.channel = 'SHZ'
                    stream = i

            # Unpack miniSEED data
            delta = mseed[stream].stats.delta
            waveform_data = mseed[stream].data

            # Extract time
            start_datetime = mseed[stream].stats.starttime.datetime
            end_datetime = mseed[stream].stats.endtime.datetime

            stn.offset = (start_datetime - setup.fireball_datetime - datetime.timedelta(minutes=5)).total_seconds()

            # Skip stations with no data
            if len(waveform_data) == 0:
                continue



            waveform_data = convolutionDifferenceFilter(waveform_data)




            # Calculate the distance from the source point to this station (kilometers)
            station_dist = greatCircleDistance(np.radians(setup.lat_centre), np.radians(setup.lon_centre), stn.position.lat_r, stn.position.lon_r)

            # Construct time array, 0 is at start_datetime
            time_data = np.arange(0, (end_datetime - start_datetime).total_seconds(), delta)

            # Cut the waveform data length to match the time data
            waveform_data = waveform_data[:len(time_data)]
            time_data = time_data[:len(waveform_data)] + stn.offset
            
            # Detrend the waveform and normalize to fixed width
            waveform_data = waveform_data - np.mean(waveform_data)

            #waveform_data = waveform_data/np.percentile(waveform_data, 99)*2
            waveform_data = waveform_data/np.max(waveform_data)*10

            # Add the distance to the waveform
            waveform_data += station_dist


            # Cut the waveforms around the time of arrival, if the window for cutting was given.
            if self.waveform_window is not None:

                # Time of arrival
                toa = station_dist/(setup.v_sound/1000) + setup.t0

                # Cut the waveform around the time of arrival
                crop_indices = (time_data >= toa - self.waveform_window/2 - 300) & (time_data <= toa + self.waveform_window/2 - 300)
                time_data = time_data[crop_indices] + 300 #HARD CODED we start 5 min before!
                waveform_data = waveform_data[crop_indices]
                

                # Skip plotting if array empty
                if len(time_data) == 0:
                    continue

            # Replace all NaNs with 0s
            waveform_data = np.nan_to_num(waveform_data, 0)
            
            max_time = np.max([max_time, np.max(time_data)])
            min_time = np.min([min_time, np.min(time_data)])

            # Keep track of minimum and maximum waveform values (used for plotting)
            max_wave_value = np.max([max_wave_value, np.max(waveform_data)])
            min_wave_value = np.min([min_wave_value, np.min(waveform_data)])
            #if data_list[idx][1].strip() not in setup.rm_stat: 
                
            # Plot the waveform on the the time vs. distance graph
            self.station_waveform[idx].setData(waveform_data*1000, time_data, pen=(255, 255, 255))
            self.make_picks_station_graph_canvas.addItem(self.station_waveform[idx])

            if stn.code in setup.high_f:
                self.fragmentation_idx.append(idx)
            if stn.code in setup.high_b:
                self.ballistic_idx.append(idx)


        toa_line_time = np.linspace(0, max_time, 10)

        # Plot the constant sound speed line (assumption is that the release happened at t = 0)
        self.make_picks_station_graph_canvas.plot((toa_line_time)*setup.v_sound, (toa_line_time + setup.t0), pen=(255, 0, 0))

        print('')
        
        self.make_picks_station_graph_canvas.setLabel('bottom', "Distance", units='m')
        self.make_picks_station_graph_canvas.setLabel('left', "Time", units='s')

        SolutionGUI.update(self)

    def keyPressEvent(self, event):
        
        if event.key() == QtCore.Qt.Key_Control:
            self.ctrl_pressed = True

        if event.key() == QtCore.Qt.Key_Shift:
            self.shift_pressed = True

        if event.key() == QtCore.Qt.Key_D:
            try:
                self.incrementStation()
            except:
                pass

        if event.key() == QtCore.Qt.Key_A:
            try:
                self.decrementStation()
            except:
                pass

        if event.key() == QtCore.Qt.Key_W:
            try:
                self.make_picks_waveform_canvas.clear()
                self.filterBandpass(event=event)
            except:
                pass

        if event.key() == QtCore.Qt.Key_S:

            try:
                self.showSpectrogram(event=event)
            except:
                pass

        if event.key() == QtCore.Qt.Key_C:
            try:
                self.make_picks_waveform_canvas.clear()
                self.filterConvolution(event=event)
            except:
                pass

    def keyReleaseEvent(self, event):

        if event.key() == QtCore.Qt.Key_Control:
            self.ctrl_pressed = False

        if event.key() == QtCore.Qt.Key_Control:
            self.ctrl_pressed = False

    def incrementStation(self, event=None):
        """ Increments the current station index. """

        self.make_picks_waveform_canvas.clear()

        self.current_station += 1

        if self.current_station >= len(self.stn_list):
            self.current_station = 0

        while self.checkExists() == False:
            self.current_station += 1
            if self.current_station >= len(self.stn_list):
                self.current_station = 0

        self.updatePlot()


    def decrementStation(self, event=None):
        """ Decrements the current station index. """

        self.make_picks_waveform_canvas.clear()

        self.current_station -= 1

        if self.current_station < 0:
            self.current_station = len(self.stn_list) - 1

        while self.checkExists() == False:
            self.current_station -= 1
            if self.current_station < 0:
                self.current_station = len(self.stn_list) - 1

        self.updatePlot()



    def markCurrentStation(self):
        """ Mark the position of the current station on the map. """


        # Init ground map
        # self.m = GroundMap(self.lat_list, self.lon_list, ax=self.map_ax, color_scheme='light')

        # for stn in self.stn_list:
        #     self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='k', s=2)

        # current_stn = self.stn_list[self.current_station]
        # self.m.scatter(current_stn.position.lat_r, current_stn.position.lon_r, c='r', s=2)

        # self.make_picks_map_graph_canvas.draw()
        SolutionGUI.update(self)


    def checkExists(self):
        """
        Checks if the current waveform is readable
        """

        # Extract current station
        stn = self.stn_list[self.current_station]

        # Get the miniSEED file path
        mseed_file_path = os.path.join(self.dir_path, stn.file_name)

        try:
            
            if os.path.isfile(mseed_file_path):
                pass
            else:
                if self.setup.debug:
                    print('File {:s} does not exist!'.format(mseed_file_path))
                return False

        except TypeError as e:
            if self.setup.debug:
                print('Opening file {:s} failed with error: {:s}'.format(mseed_file_path, str(e)))
            return False

        try:
            obspy.read(mseed_file_path)

        except TypeError:
            if self.setup.debug:
                print('mseed file could not be read:', mseed_file_path)
            return False

        return True

    def mouseClicked(self, evt):
        # class Pick:
        #     def __init__(self, time, stn, stn_no, channel):
        #         self.time = time
        #         self.stn = stn
        #         self.stn_no = stn_no
        #         self.channel = channel


        if self.ctrl_pressed:
            mousePoint = self.make_picks_waveform_canvas.vb.mapToView(evt.pos())

            self.make_picks_waveform_canvas.scatterPlot(x=[mousePoint.x()], y=[0], pen='r', update=True)

            pick = Pick(mousePoint.x(), self.stn_list[self.current_station], self.current_station, self.stn_list[self.current_station].channel)
            self.pick_list.append(pick)

            if self.setup.debug:
                print("New pick object made: {:} {:} {:}".format(mousePoint.x(), self.stn_list[self.current_station].code, self.current_station))

        elif self.shift_pressed:

            self.make_picks_waveform_canvas.clear()
            for ii, pick in enumerate(self.pick_list):
                if pick.stn_no == self.current_station:
                    self.pick_list.pop(ii)
                    if self.setup.debug:
                        print('Pick removed!')

                self.make_picks_waveform_canvas.scatterPlot(x=[pick.time], y=[0], pen='r', update=True)
            self.drawWaveform()

    def drawWaveform(self, channel_changed=0, waveform_data=None):
        """ Draws the current waveform from the current station in the waveform window. Custom waveform 
            can be given an drawn, which is used when bandpass filtering is performed. 

        """
        setup = self.setup

        # Clear waveform axis
        self.make_picks_waveform_canvas.clear()

        # Extract current station
        stn = self.stn_list[self.current_station]

        # Get the miniSEED file path
        mseed_file_path = os.path.join(self.dir_path, stn.file_name)

        # Try reading the mseed file, if it doesn't work, skip to the next frame
        try:
            mseed = obspy.read(mseed_file_path)

        except TypeError:
            if self.setup.debug:
                print('mseed file could not be read:', mseed_file_path)
            return None

        if channel_changed == 0:
            # Populate channel list
            self.make_picks_channel_choice.blockSignals(True)
            self.make_picks_channel_choice.clear()
            for i in range(len(mseed)):
                self.make_picks_channel_choice.addItem(mseed[i].stats.channel)
            self.make_picks_channel_choice.blockSignals(False)
        
        current_channel = 0
        for i in range(len(mseed)):
            if mseed[i].stats.channel == self.make_picks_channel_choice.currentText():
                self.stn_list[self.current_station].channel = mseed[i].stats.channel
                current_channel = i

        # Unpact miniSEED data
        delta = mseed[current_channel].stats.delta
        start_datetime = mseed[current_channel].stats.starttime.datetime
        end_datetime = mseed[current_channel].stats.endtime.datetime

        stn.offset = (start_datetime - setup.fireball_datetime).total_seconds()

        # Check if the waveform data is already given or not
        if waveform_data is None or channel_changed != 2:
            waveform_data = mseed[current_channel].data

            # Store raw data for bookkeeping on first open
            self.current_waveform_raw = waveform_data

        self.current_waveform_delta = delta
        self.current_waveform_time = np.arange(0, (end_datetime - start_datetime).total_seconds() + delta, \
            delta)

        # Construct time array, 0 is at start_datetime
        time_data = np.copy(self.current_waveform_time)

        # Cut the waveform data length to match the time data
        waveform_data = waveform_data[:len(time_data)]
        time_data = time_data[:len(waveform_data)] + stn.offset

        # Store currently plotted waveform
        self.current_waveform_processed = waveform_data

        # Calculate the time of arrival assuming constant propagation with the given speed of sound
        t_arrival = self.source_dists[self.current_station]/(self.v_sound/1000) + self.t0

        # Plot the waveform
        self.make_picks_waveform_canvas.plot(x=time_data, y=waveform_data, pen='w')
        self.make_picks_waveform_canvas.setXRange(t_arrival-100, t_arrival+100, padding=1)
        self.make_picks_waveform_canvas.setLabel('bottom', "Time after {:}".format(setup.fireball_datetime), units='s')
        self.make_picks_waveform_canvas.setLabel('left', "Signal Response")

        for pick in self.pick_list:
            if pick.stn_no == self.current_station:
                self.make_picks_waveform_canvas.scatterPlot(x=[pick.time], y=[0], pen='r', update=True)

        SolutionGUI.update(self)
        # Initialize variables
        b_time = 0

        # Extract coordinates of the reference station
        ref_pos = Position(setup.lat_centre, setup.lon_centre, 0)

        # Calculate ground distances
        stn.stn_ground_distance(ref_pos)

        print('####################')
        print("Current Station: {:}-{:}".format(stn.network, stn.code))
        print("Channel: {:}".format(stn.channel))
        print("Ground Distance: {:7.3f} km".format(stn.ground_distance/1000))

        # If manual ballistic search is on
        if setup.show_ballistic_waveform:

            # az = np.radians(setup.azim)
            # ze = np.radians(setup.zangle)
            # Plot Ballistic Prediction
            b_time = self.arrTimes[0, self.current_station, 0, 0]
            # sounding = parseWeather(setup, consts)
            # p = waveReleasePointWinds([stn.position.x, stn.position.y, stn.position.z], setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, az, \
            #                 ze, setup, sounding, [ref_pos.lat, ref_pos.lon, ref_pos.elev])
            # check if nan

            if b_time == b_time:
                self.make_picks_waveform_canvas.plot(x=[b_time]*2, y=[np.min(waveform_data), np.max(waveform_data)], pen=pg.mkPen(color=(0, 0, 255), width=2) , label='Ballistic')
                print("Ballistic Arrival: {:.3f} s".format(b_time))
            else:
                print("No Ballistic Arrival")

            for i in range(setup.perturb_times):
                if i >= 1:
                    if i == 1:
                        data, remove = chauvenet(self.arrTimes[:, self.current_station, 0, 0])
                        try:
                            print('Perturbation Arrival Range: {:.3f} - {:.3f}s'.format(np.nanmin(data), np.nanmax(data)))
                            print('Removed points {:}'.format(remove))
                        except ValueError:
                            print('No Perturbation Arrivals')
                    try:
                        self.make_picks_waveform_canvas.plot(x=[self.arrTimes[i, self.current_station, 0, 0]]*2, \
                         y=[np.min(waveform_data), np.max(waveform_data)], pen=pg.mkPen(color=(0, 0, 255), style=QtCore.Qt.DotLine) )
                    except:
                        pass
            # Fragmentation Prediction

            # If manual fragmentation search is on
        if setup.show_fragmentation_waveform:

            for i, line in enumerate(setup.fragmentation_point):

                f_time = self.arrTimes[0, self.current_station, 1, i]
            #     # check if nan

                print('++++++++++++++++')
                print('Fragmentation {:} ({:6.2f} km)'.format(i+1, line[2]/1000))
                if f_time == f_time:
                    # Plot Fragmentation Prediction
                    self.make_picks_waveform_canvas.plot(x=[f_time]*2, y=[np.min(waveform_data), np.max(waveform_data)], pen=pg.mkPen(color=self.pick_group_colors[(i+1)%4], width=2), label='Fragmentation')
                    stn.stn_distance(Position(line[0], line[1], line[2]))
                    print("Range: {:7.3f} km".format(stn.distance/1000))                   
                    print('Arrival: {:.3f} s'.format(f_time))

                    #if abs(f_time-395.3038374383332) <= 4:
                    #print('Fragmentation {:} Arrival: {:.3f}s / {:.3f}s'.format(i+1, f_time-282, f_time-290.5))

                else:
                    pass
                    print('No Fragmentation {:} ({:6.2f} m) Arrival'.format(i+1, line[2]))

                for j in range(setup.perturb_times):
                    if j >= 1:
                        #try:
                            if j == 1:
                                data, remove = chauvenet(self.arrTimes[:, self.current_station, 1, i])
                                try:
                                    print('Perturbation Arrival Range: {:.3f} - {:.3f}s'.format(np.nanmin(data), np.nanmax(data)))
                                    print('Removed points {:}'.format(remove))
                                except ValueError:
                                    print('No Perturbation Arrivals')
                            self.make_picks_waveform_canvas.plot(x=[self.arrTimes[j, self.current_station, 1, i]]*2, y=[np.min(waveform_data),\
                                np.max(waveform_data)], alpha=0.3,\
                                pen=pg.mkPen(color=self.pick_group_colors[(i+1)%4], style=QtCore.Qt.DotLine), zorder=3)
                        #except:
                        #    pass

    def markStationWaveform(self):
        """ Mark the currently shown waveform in the plot of all waveform. """

        pass

    def showSpectrogram(self, event=None):
        """ Show the spectrogram of the waveform in the current window. """


        # Get time limits of the shown waveform
        #x_min, x_max = self.wave_ax.get_xlim()

        # Extract the time and waveform
        #crop_window = (self.current_waveform_time >= x_min) & (self.current_waveform_time <= x_max)
        wave_arr = self.current_waveform_raw#[crop_window]


        ### Show the spectrogram ###
        
        fig = plt.figure()
        ax_spec = fig.add_subplot(111)

        ax_spec.specgram(wave_arr, Fs=1.0/self.current_waveform_delta, cmap=plt.cm.inferno)

        ax_spec.set_xlabel('Time (s)')
        ax_spec.set_ylabel('Frequency (Hz)')

        fig.show()

        ###


    def filterBandpass(self, event=None):
        """ Run bandpass filtering using values set on sliders. """

        # Get bandpass filter values
        bandpass_low = self.low_bandpass_slider.value()*self.bandpass_scale
        bandpass_high = self.high_bandpass_slider.value()*self.bandpass_scale


        # Limit the high frequency to be lower than the Nyquist frequency
        max_freq = (1.0/self.current_waveform_delta)/2

        if bandpass_high > max_freq:
            bandpass_high = max_freq - 0.1

            self.high_bandpass_slider.setValue(bandpass_high*self.bandpass_scale)
        

        # Init the butterworth bandpass filter
        butter_b, butter_a = butterworthBandpassFilter(bandpass_low, bandpass_high, \
            1.0/self.current_waveform_delta, order=6)

        # Filter the data
        waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(self.current_waveform_raw))


        # Plot the updated waveform
        self.drawWaveform(channel_changed=2, waveform_data=waveform_data)


    def filterConvolution(self, event=None):
        """ Apply the convolution filter on data as suggested in Kalenda et al. (2014). """

        waveform_data = convolutionDifferenceFilter(self.current_waveform_raw)

        self.drawWaveform(channel_changed=2, waveform_data=waveform_data)


    def updatePlot(self, draw_waveform=True):
        """ Update the plot after changes. """

        self.make_picks_waveform_canvas.clear()

        # Mark the position of the current station on the map

        self.make_picks_station_choice.setCurrentIndex(self.current_station)

        for stn_mk in self.station_marker:
            if stn_mk != self.station_marker[self.current_station]:
                stn_mk.setPen((255, 255, 255))
                stn_mk.setBrush((255, 255, 255))
                stn_mk.setZValue(0)
            elif stn_mk in [self.station_marker[i] for i in self.ballistic_idx]:
                stn_mk.setPen((0, 0, 255))
                stn_mk.setBrush((0, 0, 255))
                stn_mk.setZValue(0)
            elif stn_mk in [self.station_marker[i] for i in self.fragmentation_idx]:
                stn_mk.setPen((0, 255, 0))
                stn_mk.setBrush((0, 255, 0))
                stn_mk.setZValue(0)
            else:
                stn_mk.setPen((255, 0, 0))
                stn_mk.setBrush((255, 0, 0))
                stn_mk.setZValue(1)

        for stn_mk in self.station_waveform:
            if stn_mk != self.station_waveform[self.current_station]:
                stn_mk.setPen((255, 255, 255))
                stn_mk.setZValue(0)
            else:
                stn_mk.setPen((255, 0, 0))
                stn_mk.setZValue(1)


        # Plot the waveform from the current station
        if draw_waveform:
            self.drawWaveform()

        # Set an arrow pointing to the current station on the waveform
        self.markStationWaveform()
        self.markCurrentStation()

        # Reset bandpass filter values to default
        # self.make_picks.set_val(self.bandpass_low_default)
        # self.bandpass_high_slider.set_val(self.bandpass_high_default)

        SolutionGUI.update(self)

    def exportCSV(self, event):
        """ Save picks to a CSV file. """

        dlg = QFileDialog.getSaveFileName(self, 'Save File')


        if '.csv' not in dlg[0]:
            file_name = dlg[0] + '.csv'
        else:
            file_name = dlg[0]

        # Open the output CSV
        with open(os.path.join(file_name), 'w') as f:

            # Write the header
            f.write('Pick group, Network, Code, Lat, Lon, Elev, Pick JD, Pick time, station_number \n')

            # Go through all picks
            for pick in self.pick_list:

                # Calculate Julian date of the pick time
                pick_jd = datetime2JD(self.setup.fireball_datetime + datetime.timedelta(seconds=pick.time))

                stn = pick.stn

                # Write the CSV entry
                f.write("{:d}, {:s}, {:s}, {:.6f}, {:.6f}, {:.2f}, {:.8f}, {:}, {:}\n".format(0, stn.network, \
                    stn.code, stn.position.lat, stn.position.lon, stn.position.elev, pick_jd, pick.time, pick.stn_no))

        errorMessage('Output to CSV!', 0, title='Exported!')

    def exportToAllTimes(self):

        sounding = parseWeather(self.setup)

        if self.setup.arrival_times_file != '':
            try:
                self.arrTimes = np.load(self.setup.arrival_times_file)
            except:
                errorMessage("WARNING: Unable to load allTimes_file {:} . Please check that file exists".format(setup.arrival_times_file), 1)
                self.arrTimes = self.calcAllTimes(self.stn_list, self.setup, sounding)
        else:  
            # Calculate all arrival times
            self.arrTimes = self.calcAllTimes(self.stn_list, self.setup, sounding)
        
        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = dlg[0]

        np.save(file_name, self.arrTimes)

        errorMessage("Saved Output File", 0)

    def addMakePicksWidgets(self):
        make_picks_master_tab = QWidget()
        make_picks_master = QVBoxLayout()
        make_picks_master_tab.setLayout(make_picks_master)

        self.make_picks_top_graphs = QHBoxLayout()
        self.make_picks_bottom_graphs = QVBoxLayout()
        make_picks_master.addLayout(self.make_picks_top_graphs)
        make_picks_master.addLayout(self.make_picks_bottom_graphs)

        self.make_picks_station_graph_view = pg.GraphicsLayoutWidget()
        self.make_picks_station_graph_canvas = self.make_picks_station_graph_view.addPlot()
        self.make_picks_top_graphs.addWidget(self.make_picks_station_graph_view)
        self.make_picks_station_graph_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        self.make_picks_map_graph_view = pg.GraphicsLayoutWidget()
        self.make_picks_map_graph_canvas = self.make_picks_map_graph_view.addPlot()
        self.make_picks_top_graphs.addWidget(self.make_picks_map_graph_view)
        self.make_picks_map_graph_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        self.make_picks_waveform_view = pg.GraphicsLayoutWidget()
        self.make_picks_waveform_canvas = self.make_picks_waveform_view.addPlot()
        self.make_picks_bottom_graphs.addWidget(self.make_picks_waveform_view)
        self.make_picks_waveform_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        make_picks_control_panel = QHBoxLayout()
        self.make_picks_bottom_graphs.addLayout(make_picks_control_panel)

        make_picks_station_group = QGroupBox("Station Navigation")
        make_picks_control_panel.addWidget(make_picks_station_group)

        station_group_layout = QGridLayout()
        make_picks_station_group.setLayout(station_group_layout)

        self.make_picks_station_choice = QComboBox()
        station_group_layout.addWidget(self.make_picks_station_choice, 0, 0, 1, 2)

        self.make_picks_channel_choice = QComboBox()
        station_group_layout.addWidget(self.make_picks_channel_choice, 1, 0, 1, 2)
        self.make_picks_channel_choice.currentIndexChanged.connect(partial(self.drawWaveform, 1))

        self.prev_stat = QPushButton('Prev')
        station_group_layout.addWidget(self.prev_stat, 2, 0, 1, 1)

        self.next_stat = QPushButton('Next')
        station_group_layout.addWidget(self.next_stat, 2, 1, 1, 1)

        launch = QPushButton('Load Station Data')
        station_group_layout.addWidget(launch, 3, 0, 1, 2)
        launch.clicked.connect(self.makePicks)

        make_picks_filter_group = QGroupBox("Waveform Filtering")
        make_picks_control_panel.addWidget(make_picks_filter_group)

        filter_group_layout = QGridLayout()
        make_picks_filter_group.setLayout(filter_group_layout)

        self.low_bandpass_label = QLabel('Low: 2 Hz')
        filter_group_layout.addWidget(self.low_bandpass_label, 0, 0, 1, 1)

        self.low_bandpass_slider = QSlider(Qt.Horizontal)
        filter_group_layout.addWidget(self.low_bandpass_slider, 0, 1, 1, 1)

        self.high_bandpass_label = QLabel("High: 8 Hz")
        filter_group_layout.addWidget(self.high_bandpass_label, 1, 0, 1, 1)

        self.high_bandpass_slider = QSlider(Qt.Horizontal)
        filter_group_layout.addWidget(self.high_bandpass_slider, 1, 1, 1, 1)

        self.low_bandpass_slider.setMinimum(0/self.bandpass_scale)
        self.low_bandpass_slider.setMaximum(5/self.bandpass_scale)
        self.low_bandpass_slider.setValue(2/self.bandpass_scale)
        self.low_bandpass_slider.setTickInterval(0.5)
        self.low_bandpass_slider.valueChanged.connect(partial(self.makeValueChange, self.low_bandpass_label, self.low_bandpass_slider))

        self.high_bandpass_slider.setMinimum(3/self.bandpass_scale)
        self.high_bandpass_slider.setMaximum(40/self.bandpass_scale)
        self.high_bandpass_slider.setValue(8/self.bandpass_scale)
        self.high_bandpass_slider.setTickInterval(0.5)
        self.high_bandpass_slider.valueChanged.connect(partial(self.makeValueChange, self.high_bandpass_label, self.high_bandpass_slider))


        self.filter_combo_box = QComboBox()
        filter_group_layout.addWidget(self.filter_combo_box, 2, 0, 1, 2)

        make_picks_picks_group = QGroupBox("Arrival Picks")
        make_picks_control_panel.addWidget(make_picks_picks_group)

        pick_group_layout = QGridLayout()
        make_picks_picks_group.setLayout(pick_group_layout)

        self.export_to_csv = QPushButton('Export to CSV')
        pick_group_layout.addWidget(self.export_to_csv)

        self.export_to_all_times = QPushButton('Export All Times')
        pick_group_layout.addWidget(self.export_to_all_times)

        self.tab_widget.addTab(make_picks_master_tab, 'Make Picks')         

    def changeRows(self, obj, change):
        n_rows = obj.rowCount()
        if change == 1:
            obj.insertRow(n_rows)
        else:
            obj.removeRow(n_rows-1)

    def fileSearch(self, filters, obj):
        # obj - where the text goes
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setNameFilters(filters)
        dlg.selectNameFilter(filters[0])
        #filenames = QStringList()

        dlg.exec_()

        filename = dlg.selectedFiles()

        obj.setText(filename[0])

    def folderSearch(self, obj):
        # dlg = QFileDialog()
        # dlg.getExistingDirectory(self, "Select Directory")

        # obj - where the text goes
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.Directory)
        #dlg.setNameFilters()
        # dlg.selectNameFilter(filters[0])
        #filenames = QStringList()

        dlg.exec_()

        filename = dlg.selectedFiles()

        obj.setText(filename[0])

    def perturbSetup(self):

        if self.setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(self.setup, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(self.setup, time=-1)

        else:
            sounding_u = []
            sounding_l = []

        if self.setup.perturb_method == 'ensemble':
            ensemble_file = self.setup.perturbation_spread_file
        else:
            ensemble_file = ''

        if self.setup.perturb_times == 0: self.setup.perturb_times = 1

        if not self.setup.perturb:
            self.setup.perturb_times = 1

        return np.array([sounding_l, sounding_u, ensemble_file])

    def perturbGenerate(self, ptb_n, dataset, perturb_data):
        
        sounding_l, sounding_u, ensemble_file = perturb_data[0], perturb_data[1], perturb_data[2]

        if ptb_n > 0:
            
            if self.setup.debug:
                print("STATUS: Perturbation {:}".format(ptb_n))

            # generate a perturbed sounding profile
            sounding_p = perturbation_method(self.setup, dataset, self.setup.perturb_method, \
                sounding_u=sounding_u, sounding_l=sounding_l, \
                spread_file=self.setup.perturbation_spread_file, lat=self.setup.lat_centre, lon=self.setup.lon_centre, \
                ensemble_file=ensemble_file, ensemble_no=ptb_n)

        else:
            sounding_p = dataset

        return sounding_p

    def supSearch(self):

        s_info, s_name, weights, ref_pos = convStationDat(self.setup.station_picks_file, self.setup, d_min=self.setup.weight_distance_min, d_max=self.setup.weight_distance_max)
        ref_pos = Position(ref_pos[0], ref_pos[1], ref_pos[2])
        
        perturb_data = self.perturbSetup()

        self.setup.ref_pos = ref_pos

        n_stations = len(s_name)

        xstn = s_info[0:n_stations, 0:3]

        dataset = parseWeather(self.setup)

        results = [None]*self.setup.perturb_times

        self.setup.single_point = None

        for ptb_n in range(self.setup.perturb_times):

            sounding_p = self.perturbGenerate(ptb_n, dataset, perturb_data)

            results[ptb_n] = psoSearch(s_info, weights, s_name, self.setup, sounding_p, manual=False)

            print("Error Function: {:5.2f} (Perturbation {:})".format(results[ptb_n].f_opt, ptb_n))
            
        self.scatterPlot(self.setup, results, n_stations, xstn, s_name, dataset, manual=False)

        self.residPlot(results, s_name, xstn, self.setup.working_directory, n_stations, manual=False)

        # set row count
        self.sup_results_table.setRowCount(n_stations + 1)

        # set column count
        self.sup_results_table.setColumnCount(5)

        self.sup_results_table.setRowCount(n_stations + 1)
        self.sup_results_table.setHorizontalHeaderLabels(['Station Name', "Latitude", "Longitude", "Elevation", "Residuals"])
        header = self.sup_results_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)

        self.sup_results_table.setItem(0, 0, QTableWidgetItem("Total (Time = ({:}s)".format(results[0].motc)))
        self.sup_results_table.setItem(0, 1, QTableWidgetItem(str(results[0].x_opt.lat)))
        self.sup_results_table.setItem(0, 2, QTableWidgetItem(str(results[0].x_opt.lon)))
        self.sup_results_table.setItem(0, 3, QTableWidgetItem(str(results[0].x_opt.elev)))
        self.sup_results_table.setItem(0, 4, QTableWidgetItem(str(results[0].f_opt)))

        for i in range(n_stations):
            self.sup_results_table.setItem(i+1, 0, QTableWidgetItem(s_name[i]))
            self.sup_results_table.setItem(i+1, 1, QTableWidgetItem(str(xstn[i][0])))
            self.sup_results_table.setItem(i+1, 2, QTableWidgetItem(str(xstn[i][1])))
            self.sup_results_table.setItem(i+1, 3, QTableWidgetItem(str(xstn[i][2])))
            self.sup_results_table.setItem(i+1, 4, QTableWidgetItem(str(results[0].r[i])))

    def supraSearch(self):

        try:
            s_info, s_name, weights, ref_pos = convStationDat(self.setup.station_picks_file, self.setup, d_min=self.setup.weight_distance_min, d_max=self.setup.weight_distance_max)
        except TypeError:
            errorMessage("Unable to use station picks file!", 2, detail="Error in file: '{:}'".format(self.setup.station_picks_file))
            return None

        if self.setup.manual_fragmentation_search == None:
            errorMessage("No manual fragmentation point defined!", 2, detail="Please specify a Manual Fragmentation Search in the Sources section of Variables")        
            return None

        ref_pos = Position(ref_pos[0], ref_pos[1], ref_pos[2])
     
        dataset = parseWeather(self.setup)
        
        self.perturbSetup()

        results = [None]*self.setup.perturb_times

        for ptb_n in range(self.setup.perturb_times):

            sounding_p = self.perturbGenerate(ptb_n, dataset)

            results[ptb_n] = psoSearch(s_info, weights, s_name, self.setup, sounding_p, manual=True)
            print("Error Function: {:5.2f} (Perturbation {:})".format(results[ptb_n].f_opt, ptb_n))
        
        n_stations = len(s_info)
        xstn = s_info[0:n_stations, 0:3]
            
        self.scatterPlot(self.setup, results, n_stations, xstn, s_name, dataset)

        self.residPlot(results, s_name, xstn, self.setup.working_directory, n_stations)

        # set row count
        self.tableWidget.setRowCount(n_stations + 1)
        self.tableWidget.setHorizontalHeaderLabels(['Station Name', "Latitude", "Longitude", "Elevation", "Residuals"])
        header = self.tableWidget.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)

        # set column count
        self.tableWidget.setColumnCount(5)

        self.tableWidget.setItem(0, 0, QTableWidgetItem("Total"))
        self.tableWidget.setItem(0, 1, QTableWidgetItem(str(results[0].x_opt.lat)))
        self.tableWidget.setItem(0, 2, QTableWidgetItem(str(results[0].x_opt.lon)))
        self.tableWidget.setItem(0, 3, QTableWidgetItem(str(results[0].x_opt.elev)))
        self.tableWidget.setItem(0, 4, QTableWidgetItem(str(results[0].f_opt)))

        for i in range(n_stations):
            self.tableWidget.setItem(i+1, 0, QTableWidgetItem(s_name[i]))
            self.tableWidget.setItem(i+1, 1, QTableWidgetItem(str(xstn[i][0])))
            self.tableWidget.setItem(i+1, 2, QTableWidgetItem(str(xstn[i][1])))
            self.tableWidget.setItem(i+1, 3, QTableWidgetItem(str(xstn[i][2])))
            self.tableWidget.setItem(i+1, 4, QTableWidgetItem(str(results[0].r[i])))


    def saveINI(self, write=True, autosave=False):

        """if write == True - used for saving data to setup obj, but not writing
        """
        if write and not autosave:
            dlg = QFileDialog.getSaveFileName(self, 'Save File')

        self.setup.fireball_name = self.fireball_name_edits.text()
        self.setup.get_data = tryBool(self.get_data_edits.currentText())
        self.setup.run_mode = self.run_mode_edits.currentText()
        self.setup.debug = tryBool(self.debug_edits.currentText())

        self.setup.working_directory = self.working_directory_edits.text()
        self.setup.arrival_times_file = self.arrival_times_edits.text()
        self.setup.sounding_file = self.sounding_file_edits.text()
        self.setup.perturbation_spread_file = self.perturbation_file_edits.text()
        self.setup.station_picks_file = self.station_picks_edits.text()
        self.setup.replot_points_file = self.points_name_edits.text()

        self.setup.lat_centre = tryFloat(self.lat_centre_edits.text())
        self.setup.lon_centre = tryFloat(self.lon_centre_edits.text())
        self.setup.deg_radius = tryFloat(self.deg_radius_edits.text())
        self.setup.fireball_datetime = self.fireball_datetime_edits.dateTime().toPyDateTime()
        self.setup.v_sound = tryFloat(self.v_sound_edits.text())

        self.setup.t0 = tryFloat(self.t0_edits.text())
        self.setup.v = tryFloat(self.v_edits.text())
        self.setup.azimuth = tryAngle(self.azim_edits.text())
        self.setup.zenith = tryAngle(self.zangle_edits.text())
        self.setup.lat_i = tryFloat(self.lat_i_edits.text())
        self.setup.lon_i = tryFloat(self.lon_i_edits.text())
        self.setup.elev_i = tryFloat(self.elev_i_edits.text())
        self.setup.lat_f = tryFloat(self.lat_f_edits.text())
        self.setup.lon_f = tryFloat(self.lon_f_edits.text())
        self.setup.elev_f = tryFloat(self.elev_f_edits.text())

        self.setup.pos_i = tryPosition(self.setup.lat_i, self.setup.lon_i, self.setup.elev_i)
        self.setup.pos_f = tryPosition(self.setup.lat_f, self.setup.lon_f, self.setup.elev_f)
        self.trajectory = tryTrajectory(self.setup.t0, self.setup.v, self.setup.azimuth, self.setup.zenith, self.setup.pos_i, self.setup.pos_f)

        self.setup.show_ballistic_waveform = tryBool(self.show_ballistic_waveform_edits.currentText())

        self.setup.fragmentation_point = trySupracenter(fromTable(self.fragmentation_point))
        self.setup.show_fragmentation_waveform = tryBool(self.show_fragmentation_waveform_edits.currentText())

        frag_lat  = tryFloat(self.lat_frag_edits.text())
        frag_lon  = tryFloat(self.lon_frag_edits.text())
        frag_elev = tryFloat(self.elev_frag_edits.text())
        frag_time = tryFloat(self.time_frag_edits.text())
        self.setup.manual_fragmentation_search = trySupracenter(tryPosition(frag_lat, frag_lon, frag_elev), frag_time)

        self.setup.v_fixed = tryFloat(self.v_fixed_edits.text())
        self.setup.enable_restricted_time = self.restricted_time_check.isChecked()
        self.setup.restricted_time = self.restricted_time_edits.dateTime().toPyDateTime()

        self.setup.azimuth_min = tryAngle(self.azimuth_min_edits.text())
        self.setup.azimuth_max = tryAngle(self.azimuth_max_edits.text())
        self.setup.zenith_min = tryAngle(self.zangle_min_edits.text())
        self.setup.zenith_max = tryAngle(self.zangle_max_edits.text())
        self.setup.lat_min = tryFloat(self.lat_min_edits.text())
        self.setup.lat_max = tryFloat(self.lat_max_edits.text())
        self.setup.lon_min = tryFloat(self.lon_min_edits.text())
        self.setup.lon_max = tryFloat(self.lon_max_edits.text())
        self.setup.elev_min = tryFloat(self.elev_min_edits.text())
        self.setup.elev_max = tryFloat(self.elev_max_edits.text())
        self.setup.t_min = tryFloat(self.t_min_edits.text())
        self.setup.t_max = tryFloat(self.t_max_edits.text())
        self.setup.v_min = tryFloat(self.v_min_edits.text())
        self.setup.v_max = tryFloat(self.v_max_edits.text())
        self.setup.weight_distance_min = tryFloat(self.weight_distance_min_edits.text())
        self.setup.weight_distance_max = tryFloat(self.weight_distance_max_edits.text())

        self.setup.pos_min = tryPosition(self.setup.lat_min, self.setup.lon_min, self.setup.elev_min)
        self.setup.pos_max = tryPosition(self.setup.lat_max, self.setup.lon_max, self.setup.elev_max)

        self.setup.supra_min = trySupracenter(self.setup.pos_min, self.setup.t_min)
        self.setup.supra_max = trySupracenter(self.setup.pos_max, self.setup.t_max)

        self.setup.traj_min = tryTrajectory(self.setup.t_min, self.setup.v_min, self.setup.azimuth_min, \
                                         self.setup.zenith_min, self.setup.pos_min, self.setup.pos_min)
        self.setup.traj_max = tryTrajectory(self.setup.t_max, self.setup.v_max, self.setup.azimuth_max, \
                                         self.setup.zenith_max, self.setup.pos_max, self.setup.pos_max)

        self.setup.enable_winds = tryBool(self.enable_winds_edits.currentText())
        self.setup.weather_type = self.weather_type_edits.currentText()

        self.setup.perturb_times = tryInt(self.perturb_times_edits.text())
        self.setup.observe_frag_no = tryInt(self.frag_no_edits.text())
        self.setup.perturb = tryBool(self.perturb_edits.currentText())
        self.setup.perturb_method = self.perturb_method_edits.currentText()

        self.setup.n_theta = tryInt(self.n_theta_edits.text())
        self.setup.n_phi = tryInt(self.n_phi_edits.text())
        self.setup.angle_precision = tryFloat(self.angle_precision_edits.text())
        self.setup.angle_error_tol = tryFloat(self.angle_error_tol_edits.text())

        self.setup.maxiter = tryInt(self.maxiter_edits.text())
        self.setup.swarmsize = tryInt(self.swarmsize_edits.text())
        self.setup.run_times = tryInt(self.run_times_edits.text())
        self.setup.minfunc = tryFloat(self.minfunc_edits.text())
        self.setup.minstep = tryFloat(self.minstep_edits.text())
        self.setup.phip = tryFloat(self.phip_edits.text())
        self.setup.phig = tryFloat(self.phig_edits.text())
        self.setup.omega = tryFloat(self.omega_edits.text())
        self.setup.pso_debug = tryBool(self.pso_debug_edits.currentText())

        self.setup.contour_res = tryInt(self.contour_res_edits.text())
        self.setup.high_f = self.high_f_edits.text()
        self.setup.high_b = self.high_b_edits.text()
        self.setup.rm_stat = self.rm_stat_edits.text()

        self.setup.stations = []
        temp_stats = fromTable(self.extra_point)
        for stn in temp_stats:
            self.setup.stations.append(Station(stn[0], stn[1], Position(stn[2], stn[3], stn[4]), stn[5], stn[6], stn[7]))

        if autosave:
            with open(self.setup_file, 'wb') as f:
                pickle.dump(self.setup, f)

        else:
            if write:
                if '.pkl' not in dlg[0]:
                    output = dlg[0] + '.pkl'
                else:
                    output = dlg[0]

                with open(output, 'wb') as f:
                    pickle.dump(self.setup, f)


        self.setup = saveDefaults(self.setup)

    def loadINI(self):
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setNameFilters(['Pickle File (*.pkl)'])
        #filenames = QStringList()

        dlg.exec_()

        filename = dlg.selectedFiles()

        if len(filename) == 0:
            errorMessage("No File Selected", 1)
            return None

        if '.pkl' not in filename[0]:
            errorMessage("No File Selected!", 1)
            return None
        
        self.setup_file = filename[0]

        with open(filename[0], 'rb') as f:
            self.setup = pickle.load(f)

        self.setup = saveDefaults(self.setup)

        self.fireball_name_edits.setText(self.setup.fireball_name)
        comboSet(self.get_data_edits, self.setup.get_data)
        comboSet(self.run_mode_edits, self.setup.run_mode)
        comboSet(self.debug_edits, self.setup.debug)

        self.working_directory_edits.setText(self.setup.working_directory)
        self.arrival_times_edits.setText(self.setup.arrival_times_file)
        self.sounding_file_edits.setText(self.setup.sounding_file)
        self.perturbation_file_edits.setText(self.setup.perturbation_spread_file)
        self.station_picks_edits.setText(self.setup.station_picks_file)
        self.points_name_edits.setText(self.setup.replot_points_file)

        self.lat_centre_edits.setText(str(self.setup.lat_centre))
        self.lon_centre_edits.setText(str(self.setup.lon_centre))
        self.deg_radius_edits.setText(str(self.setup.deg_radius))
        self.fireball_datetime_edits.setDateTime(self.setup.fireball_datetime)
        self.v_sound_edits.setText(str(self.setup.v_sound))

        self.t0_edits.setText(str(self.setup.t0))
        self.v_edits.setText(str(self.setup.v))
        self.azim_edits.setText(str(self.setup.azimuth.deg))
        self.zangle_edits.setText(str(self.setup.zenith.deg))
        self.lat_i_edits.setText(str(self.setup.lat_i))
        self.lon_i_edits.setText(str(self.setup.lon_i))
        self.elev_i_edits.setText(str(self.setup.elev_i))
        self.lat_f_edits.setText(str(self.setup.lat_f))
        self.lon_f_edits.setText(str(self.setup.lon_f))
        self.elev_f_edits.setText(str(self.setup.elev_f))
        comboSet(self.show_ballistic_waveform_edits, self.setup.show_ballistic_waveform)

        frag_list = []
        for element in self.setup.fragmentation_point:
            frag_list.append(element.toList())

        toTable(self.fragmentation_point, frag_list)
        comboSet(self.show_fragmentation_waveform_edits, self.setup.show_fragmentation_waveform)

        if self.setup.manual_fragmentation_search == None:
            self.lat_frag_edits.setText('')
            self.lon_frag_edits.setText('')
            self.elev_frag_edits.setText('')
            self.time_frag_edits.setText('')
        else: 
            self.lat_frag_edits.setText(str(self.setup.manual_fragmentation_search[0].position.lat))
            self.lon_frag_edits.setText(str(self.setup.manual_fragmentation_search[0].position.lon))
            self.elev_frag_edits.setText(str(self.setup.manual_fragmentation_search[0].position.elev))
            self.time_frag_edits.setText(str(self.setup.manual_fragmentation_search[0].time))

        self.v_fixed_edits.setText(str(self.setup.v_fixed))
        self.restricted_time_check.setChecked(self.setup.enable_restricted_time)
        self.restricted_time_edits.setDateTime(self.setup.restricted_time)

        self.azimuth_min_edits.setText(str(self.setup.azimuth_min.deg))
        self.azimuth_max_edits.setText(str(self.setup.azimuth_max.deg))
        self.zangle_min_edits.setText(str(self.setup.zenith_min.deg))
        self.zangle_max_edits.setText(str(self.setup.zenith_max.deg))
        self.lat_min_edits.setText(str(self.setup.lat_min))
        self.lat_max_edits.setText(str(self.setup.lat_max))
        self.lon_min_edits.setText(str(self.setup.lon_min))
        self.lon_max_edits.setText(str(self.setup.lon_max))
        self.elev_min_edits.setText(str(self.setup.elev_min))
        self.elev_max_edits.setText(str(self.setup.elev_max))
        self.t_min_edits.setText(str(self.setup.t_min))
        self.t_max_edits.setText(str(self.setup.t_max))
        self.v_min_edits.setText(str(self.setup.v_min))
        self.v_max_edits.setText(str(self.setup.v_max))
        self.weight_distance_min_edits.setText(str(self.setup.weight_distance_min))
        self.weight_distance_max_edits.setText(str(self.setup.weight_distance_max))

        comboSet(self.enable_winds_edits, self.setup.enable_winds)
        comboSet(self.weather_type_edits, self.setup.weather_type)

        self.perturb_times_edits.setText(str(self.setup.perturb_times))
        self.frag_no_edits.setText(str(self.setup.observe_frag_no))
        comboSet(self.perturb_edits, self.setup.perturb)
        comboSet(self.perturb_method_edits, self.setup.perturb_method)

        self.n_theta_edits.setText(str(self.setup.n_theta))
        self.n_phi_edits.setText(str(self.setup.n_phi))
        self.angle_precision_edits.setText(str(self.setup.angle_precision))
        self.angle_error_tol_edits.setText(str(self.setup.angle_error_tol))

        self.maxiter_edits.setText(str(self.setup.maxiter))
        self.swarmsize_edits.setText(str(self.setup.swarmsize))
        self.run_times_edits.setText(str(self.setup.run_times))
        self.minfunc_edits.setText(str(self.setup.minfunc))
        self.minstep_edits.setText(str(self.setup.minstep))
        self.phip_edits.setText(str(self.setup.phip))
        self.phig_edits.setText(str(self.setup.phig))
        self.omega_edits.setText(str(self.setup.omega))
        comboSet(self.pso_debug_edits, self.setup.pso_debug)

        self.contour_res_edits.setText(str(self.setup.contour_res))
        self.high_f_edits.setText(str(self.setup.high_f))
        self.high_b_edits.setText(str(self.setup.high_b))
        self.rm_stat_edits.setText(str(self.setup.rm_stat))

        toTableFromStn(self.extra_point, self.setup.stations)


    def scatterPlot(self, setup, results, n_stations, xstn, s_name, dataset, manual=True):

        """ Outputs a scatter plot of the search data and the optimal supracenter

        Arguments:
            single_point: [list] lat, lon, and height of a manual search point
            n_stations: [int] number of stations
            xstn: [list] lat, lon, and height of all the stations
            s_name: [list] name of all the stations
            r: [list] residuals to each station
            x_opt: [list] lat, lon, and height of the optimal supracenter
            reported_points: [list] name, lat, lon, and height of extra reference points to be plotted
            search: [list] min/max lat, lon, and height of the search boundary
            output_name: [string] folder name to save output files in
            ref_pos: [list] mean position of the stations to act as the origin for the local coordinate system
            sup: [ndarray] array of points searched with the search algorithm used for plotting
            errors: [ndarray] array of the errors in each of the points in sup for plotting

        Returns:
            min_search: [list] min lat, lon and height of the search area
            max_search: [list] max lat, lon and height of the search area
        """

        plt.style.use('dark_background')
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(20.9, 11.7)
        ax = fig.add_subplot(1, 1, 1, projection='3d')

        ### Labels
        ax.set_title("Supracenter Locations")
        ax.set_xlabel("Latitude (deg N)", linespacing=3.1)
        ax.set_ylabel("Longitude (deg E)", linespacing=3.1)
        ax.set_zlabel('Elevation (m)', linespacing=3.1)

        # plot station names and residuals
        for h in range(n_stations):

            # Convert station locations to geographic
            xstn[h, 0], xstn[h, 1], xstn[h, 2] = loc2Geo(self.setup.lat_centre, self.setup.lon_centre, 0, xstn[h, :])

            # Add station names
            ax.text(xstn[h, 0], xstn[h, 1], xstn[h, 2],  '%s' % (s_name[h]), size=10, zorder=1, color='w')

        for ptb in range(setup.perturb_times):
            for h in range(n_stations):

                xline = []
                yline = []
                zline = []
                try:
                    for line in results[ptb].trace[h]:
                        line[0], line[1], line[2] = loc2Geo(xstn[h, 0], xstn[h, 1], xstn[h, 2], [line[0], line[1], line[2]])

                        xline.append(line[0])
                        yline.append(line[1])
                        zline.append(line[2])
                except:            
                    if not perturb:
                        errorMessage('Cannot trace rays!', 2)
                        return None
                if manual:
                    if ptb == 0:
                        ax.plot3D(xline, yline, zline, 'white')
                    else:
                        ax.plot3D(xline, yline, zline, 'green')
                elif not setup.perturb:
                    if ptb == 0:
                        ax.plot3D(xline, yline, zline, 'white')
                else:
                    pass


        r = results[0].r
        x_opt = results[0].x_opt
        errors = results[0].errors
        sup = results[0].sup

        c = np.nan_to_num(r)

        lat_list = []
        lon_list = []
        elev_list = []

        # Add stations with color based off of residual
        ax.scatter(xstn[:, 0], xstn[:, 1], xstn[:, 2], c=abs(c), marker='^', cmap='viridis_r', depthshade=False)

        for ptb in range(setup.perturb_times):
            
            ax.scatter(x_opt.lat, x_opt.lon, x_opt.elev, c = 'r', marker='*')
            ax.scatter(results[ptb].x_opt.lat, results[ptb].x_opt.lon, results[ptb].x_opt.elev, c = 'g', marker='*', alpha=0.7)

            lat_list.append(results[ptb].x_opt.lat)   
            lon_list.append(results[ptb].x_opt.lon)  
            elev_list.append(results[ptb].x_opt.elev)   

        if self.setup.perturb:

            # Plot the surface
            a = np.nanmax(((x_opt.lat - lat_list)**2 + (x_opt.lon - lon_list)**2)**0.5)
            r = np.nanmax(abs(x_opt.elev - elev_list))

            x, y, z = sphereData(a, a, r, x_opt.lat, x_opt.lon, x_opt.elev)

            # uncertainty sphere
            ax.plot_wireframe(x, y, z, color='r', alpha=0.5)

        ax.text(x_opt.lat, x_opt.lon, x_opt.elev, '%s' % ('Supracenter'), zorder=1, color='w')

        if not manual:
            for i in range(len(sup)):
                sup[i, 0], sup[i, 1], sup[i, 2] = loc2Geo(self.setup.ref_pos.lat, self.setup.ref_pos.lon, self.setup.ref_pos.elev, sup[i, :])
            sc = ax.scatter(sup[:, 0], sup[:, 1], sup[:, 2], c=errors, cmap='inferno_r', depthshade=False)
            a = plt.colorbar(sc, ax=ax)
            a.set_label("Error in Supracenter (s)")

        if manual:
            try:
                self.plots.removeWidget(self.threelbar)
            except:
                pass
            self.plots.removeWidget(self.three_canvas)
            #self.three_canvas = FigureCanvas(Figure(figsize=(2, 2)))
            self.three_canvas = FigureCanvas(fig)
            self.three_canvas.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
            self.threelbar = NavigationToolbar(self.three_canvas, self)
            self.plots.addWidget(self.three_canvas)
            self.plots.addWidget(self.threelbar)
            self.three_canvas.draw()
        else:
            try:
                self.sup_plots.removeWidget(self.sup_threelbar)
            except:
                pass
            self.sup_plots.removeWidget(self.sup_three_canvas)
            #self.sup_three_canvas = FigureCanvas(Figure(figsize=(2, 2)))
            self.sup_three_canvas = FigureCanvas(fig)
            self.sup_three_canvas.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
            self.sup_threelbar = NavigationToolbar(self.sup_three_canvas, self)
            self.sup_plots.addWidget(self.sup_three_canvas)
            self.sup_plots.addWidget(self.sup_threelbar)
            self.sup_three_canvas.draw()

        ax.mouse_init()
        SolutionGUI.update(self)

    def residPlot(self, results_arr, s_name, xstn, output_name, n_stations, manual=True):
        """ outputs a 2D residual plot of the stations with the optimal supracenter

        Arguments:
            x_opt: [list] optimal supracenter position
            s_name: [list] list of station names
            xstn: [list] list of lat, lon, height positions of each station
            resid: [list] list of residuals to each station
            output_name: [string] folder to store the data in
            n_stations: [int] number of stations
        """
        
        x_opt = results_arr[0].x_opt
        resid = results_arr[0].r

        plt.style.use('dark_background')
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(20.9, 11.7)
        ax = fig.add_subplot(1, 1, 1)
        res = ax.scatter(xstn[:, 0], xstn[:, 1], c=abs(resid), marker='^', cmap='viridis_r', s=21)
        
        for h in range(n_stations):

            # Add station names
            ax.text(xstn[h, 0], xstn[h, 1],  '%s' % (s_name[h]), size=10, zorder=1, color='w')

            if not np.isnan(resid[h]):
                
                xline = []
                yline = []
                zline = []
                try:
                    for line in results_arr[h].trace:

                        line[0], line[1], line[2] = loc2Geo(xstn[h, 0], xstn[h, 1], xstn[h, 2], [line[0], line[1], line[2]])

                        xline.append(line[0])
                        yline.append(line[1])
                        zline.append(line[2])

                    ax.plot(xline, yline, c='w')
                except:            
                    if not perturb:
                        errorMessage('Cannot trace rays!', 2)
                        return None
        

        ax.text(x_opt.lat, x_opt.lon,  '%s' % ('Supracenter'), size=10, zorder=1, color='w')

        lat_list = []
        lon_list = []

        for ptb in range(self.setup.perturb_times):
            ax.scatter(x_opt.lat, x_opt.lon, c = 'r', marker='*', s=21)
            ax.scatter(results_arr[ptb].x_opt.lat, results_arr[ptb].x_opt.lon, c="g", marker='*', s=21, alpha=0.7)

            lat_list.append(results_arr[ptb].x_opt.lat)   
            lon_list.append(results_arr[ptb].x_opt.lon)  

        if self.setup.perturb:

            # Plot the surface
            a = np.nanmax(((x_opt.lat - lat_list)**2 + (x_opt.lon - lon_list)**2)**0.5)

            circle = plt.Circle((x_opt.lat, x_opt.lon), a, color='r', alpha=0.3)
            ax.add_artist(circle)

        ax.set_xlabel("Latitude (deg N)")
        ax.set_ylabel("Longitude (deg E)")
        ax.set_title("Station Residuals")

        c = plt.colorbar(res, ax=ax)
        c.set_label("Station Residuals (s)")

        if manual:
            self.plots.removeWidget(self.two_canvas)
            try:
                self.plots.removeWidget(self.twolbar)
            except:
                pass
            #self.two_canvas = FigureCanvas(Figure(figsize=(2, 2)))
            self.two_canvas = FigureCanvas(fig)
            self.two_canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
            self.twolbar = NavigationToolbar(self.two_canvas, self)
            self.plots.addWidget(self.two_canvas)    
            self.plots.addWidget(self.twolbar)
            self.two_canvas.draw()
        else:

            self.sup_plots.removeWidget(self.sup_two_canvas)
            try:
                self.sup_plots.removeWidget(self.sup_twolbar)
            except:
                pass
            #self.sup_two_canvas = FigureCanvas(Figure(figsize=(2, 2)))
            self.sup_two_canvas = FigureCanvas(fig)
            self.sup_two_canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
            self.sup_twolbar = NavigationToolbar(self.sup_two_canvas, self)
            self.sup_plots.addWidget(self.sup_two_canvas)    
            self.sup_plots.addWidget(self.sup_twolbar)
            self.sup_two_canvas.draw()
        SolutionGUI.update(self)

if __name__ == '__main__':
    print('#########################################')
    print('#     Western Meteor Python Library     #')
    print('# Seismic and Infrasonic Meteor Program #')
    print('#      /// Legendary Edition ///        #')
    print('#            Luke McFadden,             #')
    print('#              Denis Vida,              #') 
    print('#              Peter Brown              #')
    print('#                 2019                  #')
    print('#########################################')
    app = QApplication(sys.argv)

    gui = SolutionGUI()

    gui.showFullScreen()
    gui.show()

    app.exec_()
