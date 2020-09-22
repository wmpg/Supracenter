################################################
# Credits:
# Peter Brown - Supervisor
# Luke McFadden - General coding
# Denis Vida - Ballistic code, WMPL
# Wayne Edwards - Supracenter code
# Elizabeth Silber - Updated Supracenter code
# Gunter Stober - Advice on atmospheric profiles
# Stack Overflow - Frequent care and support
# Western Meteor Python Group
#################################################
# https://docs.obspy.org/tutorial/code_snippets/seismometer_correction_simulation.html#calculating-response-from-filter-stages-using-evalresp
# t_0 = impact time - fireball datetime

import os
import time
import datetime
import copy
import webbrowser

import zipfile
import pickle

from netCDF4 import Dataset
from PyQt5.QtWidgets import *
from functools import partial
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
import pyqtgraph.exporters

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Fireballs.SeismicTrajectory import timeOfArrival, trajSearch, getStationList, estimateSeismicTrajectoryAzimuth, plotStationsAndTrajectory, waveReleasePointWindsContour

from supra.Supracenter.slowscan2 import cyscan as slowscan
from supra.Supracenter.psoSearch import psoSearch
from supra.Supracenter.fetchCopernicus import copernicusAPI
from supra.Supracenter.cyscan2 import cyscan
from supra.Supracenter.cyscanVectors import cyscan as cyscanV
from supra.Supracenter.propegateBackwards import propegateBackwards

from supra.GUI.Tools.GUITools import *
from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.WidgetBuilder import *
from supra.GUI.Dialogs.AnnoteWindow import AnnoteWindow
from supra.GUI.Dialogs.Preferences import PreferenceWindow
from supra.GUI.Dialogs.Yields import Yield
from supra.GUI.Dialogs.FragStaff import FragmentationStaff
from supra.GUI.Dialogs.AllWaveformView import AllWaveformViewer
from supra.GUI.Dialogs.TrajInterp import TrajInterpWindow
from supra.GUI.Dialogs.StationList import StationList
from supra.GUI.Dialogs.ParticleMot import ParticleMotion
from supra.GUI.Dialogs.Polmap import Polmap
from supra.GUI.Dialogs.BandpassGUI import BandpassWindow
from supra.GUI.Dialogs.ReportDialog import ReportWindow
from supra.GUI.Tools.htmlLoader import htmlBuilder
from supra.GUI.Tools.Errors import errorCodes

from supra.GUI.Tabs.SupracenterSearch import supSearch

from supra.Stations.Filters import *
from supra.Stations.CalcAllTimes4 import calcAllTimes
from supra.Stations.StationObj import Polarization

from wmpl.Utils.TrajConversions import datetime2JD, jd2Date
from wmpl.Utils.Earth import greatCircleDistance

from supra.Utils.AngleConv import loc2Geo, chauvenet, angle2NDE
from supra.Utils.Formatting import *
from supra.Utils.Classes import Position, Constants, Pick, RectangleItem, Color, Plane, Annote
from supra.Utils.TryObj import *
from supra.Utils.pso import pso

from supra.Files.SaveObjs import Prefs, BAMFile
from supra.Files.SaveLoad import save, load

from supra.Atmosphere.Parse import parseWeather
from supra.Atmosphere.radiosonde import downloadRadio

from supra.Supracenter.l137 import estPressure

HEIGHT_SOLVER_DIV = 250
THEO = False

PEN = [(0     *255, 0.4470*255, 0.7410*255),        
       (0.8500*255, 0.3250*255, 0.0980*255),            
       (0.9290*255, 0.6940*255, 0.1250*255),          
       (0.4940*255, 0.1840*255, 0.5560*255),                
       (0.4660*255, 0.6740*255, 0.1880*255),                
       (0.3010*255, 0.7450*255, 0.9330*255),                
       (0.6350*255, 0.0780*255, 0.1840*255)]



def contourLoop(X, Y, ref_pos, dy, dx, T, i):

    A = Position(0, 0, 0)
    A.x = X[i]
    A.y = Y[i]
    A.z = 0
    A.pos_geo(ref_pos)

    # return data in a form readable by Rectangle Object
    return (A.lon, A.lat, dy, dx, T[i]) 


# Main Window
class SolutionGUI(QMainWindow):

    def __init__(self):
        super().__init__()
        
        ##############################
        # Load system-wide preferences
        ##############################
        self.prefs = Prefs()
        try:
            with open(os.path.join('supra', 'Misc', 'BAMprefs.bam'), 'rb') as f:
                self.prefs = pickle.load(f)
        except EOFError as e:
            # Prefs file missing? - Not sure how to handle this yet
            pass

        self.bam = BAMFile()
        self.color = Color()

        # Initialize all of the pyqt things in the GUI
        initMainGUI(self)
        initMainGUICosmetic(self)

        # Add widgets to the floating box
        self.addIniDockWidgets()


    def viewToolbar(self):

        # Toggles the toolbar
        self.ini_dock.toggleViewAction().trigger()


    def viewFullscreen(self):

        # Toggles fullscreen
        if self.windowState() & QtCore.Qt.WindowFullScreen:
            self.showNormal()
        else:
            self.showFullScreen()


    def quitApp(self):

        # Begins quit sequence
        reply = QMessageBox.question(self, 'Quit Program', 'Are you sure you want to quit?', QMessageBox.Yes, QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            qApp.quit()
        else:
            return None

        
    def openGit(self):

        webbrowser.open_new_tab("https://github.com/dvida/Supracenter")


    def openDocs(self):

        # docs are a locally stored html file
        webbrowser.open_new_tab(self.doc_file)
        
    def genReport(self):

        self.gr = ReportWindow(self.bam, self.prefs)
        self.gr.setGeometry(QRect(500, 400, 500, 400))
        self.gr.show()

    def stndownloadDialog(self):

        self.sd = StationList()
        self.sd.setGeometry(QRect(500, 400, 500, 400))
        self.sd.show()

    def preferencesDialog(self):

        self.p = PreferenceWindow()
        self.p.setGeometry(QRect(500, 400, 500, 400))
        self.p.show()

    def trajInterpDialog(self):

        self.t = TrajInterpWindow(self.bam.setup)
        self.t.setGeometry(QRect(500, 400, 500, 400))
        self.t.show()

    def csvLoad(self, table):
        """ Loads csv file into a table
        """

        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setNameFilters(['CSV File (*.csv)'])
        dlg.exec_()

        filename = dlg.selectedFiles()

        try:
            with open(filename[0]) as f:
                data_table = []
                next(f)
        
                for line in f:
                    a = line.split(',')
                    if len(a) != 9:
                        errorMessage('Wrong number of columns for a picks file!', 1, info='Make sure a picks file is imported!')
                        return None
                    data_table.append(a)


        except IsADirectoryError as e:
            errorMessage('Please select a valid file to load', 1, detail='{:}'.format(e))
            return None

        defTable(self.csv_table, 0, 9, headers=['Pick Group', 'Network', 'Code', 'Latitude', 'Longitude', 'Elevation', 'Pick JD', 'Pick Time', 'Station Number'])

        toTable(table, data_table)


    def csvSave(self, table):
        """  Saves a table to a csv
        """

        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = checkExt(dlg[0], '.csv')

        data_set = fromTable(table)
        # Open the output CSV
        with open(os.path.join(file_name), 'w') as f:

            # Write the header
            f.write('Pick group, Network, Code, Lat, Lon, Elev, Pick JD, Pick time, station_number \n')

            # Go through all picks
            for line in data_set:
                line[-1] = int(line[-1])
                # Write the CSV entry
                f.write("{:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}\n".format(*line))

        errorMessage('Output to CSV!', 0, title='Exported!', detail='Filename: {:}'.format(file_name))

    def supSearchSetup(self, manual):

        supSearch(self, self.bam, self.prefs, manual=manual, results_print=False)

    def psoTrajectory(self, station_list, sounding):

        ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)

        self.setup.pos_min.pos_loc(ref_pos)
        self.setup.pos_max.pos_loc(ref_pos)

        bounds = [
            (self.setup.pos_min.x, self.setup.pos_max.x), # X0
            (self.setup.pos_min.y, self.setup.pos_max.y), # Y0
            (self.setup.t_min, self.setup.t_max), # t0
            (self.setup.v_min, self.setup.v_max), # Velocity (m/s)
            (self.setup.azimuth_min.deg, self.setup.azimuth_max.deg),     # Azimuth
            (self.setup.zenith_min.deg, self.setup.zenith_max.deg)  # Zenith angle
            ]

        lower_bounds = [bound[0] for bound in bounds]
        upper_bounds = [bound[1] for bound in bounds]

        try:
            plane = fromTable(self.seis_plane)
            
            p = []

            for ii, point in enumerate(plane):

                p.append(Position(point[0], point[1], point[2]))
                p[ii].pos_loc(ref_pos)

            rest_plane = Plane(p[0].xyz, p[1].xyz, p[2].xyz)
        except:
            rest_plane = None
    
        if rest_plane == None:
            if self.prefs.debug:
                print('Free Search')
            x, fopt = pso(trajSearch, lower_bounds, upper_bounds, args=(station_list, sounding, ref_pos, self.setup, rest_plane), \
                maxiter=self.setup.maxiter, swarmsize=self.setup.swarmsize, \
                phip=self.setup.phip, phig=self.setup.phig, debug=False, omega=self.setup.omega, \
                particle_output=False)
        else:
            if self.prefs.debug:
                print('Plane Search')
            x, fopt = pso(trajSearch, lower_bounds, upper_bounds, ieqcons=[planeConst], args=(station_list, sounding, ref_pos, self.setup, rest_plane), \
                maxiter=self.setup.maxiter, swarmsize=self.setup.swarmsize, \
                phip=self.setup.phip, phig=self.setup.phig, debug=False, omega=self.setup.omega, \
                particle_output=False)

        print('Results:')
        print('X: {:.4f}'.format(x[0]))
        print('Y: {:.4f}'.format(x[1]))
        print('t: {:.4f}'.format(x[2]))
        print('v: {:.4f}'.format(x[3]))
        print('az: {:.4f}'.format(x[4]))
        print('ze: {:.4f}'.format(x[5]))
        print('err: {:.4f}'.format(fopt))

        geo = Position(0, 0, 0)
        geo.x = x[0]
        geo.y = x[1]
        geo.z = 0
        geo.pos_geo(ref_pos)

        print('Geometric Landing Point:')
        print(geo)

    def seisSearch(self):

        # Read station file

        try:
            station_list = getStationList(os.path.join(self.prefs.workdir, self.setup.fireball_name, self.setup.station_picks_file))
        except TypeError as e:
            errorMessage('Unexpected station list location!', 2, info="Can not find where 'station_picks_file' is!", detail='{:}'.format(e))
            return None

        class StationPick:
            def __init__(self):
                pass

        station_obj_list = []
        for i in range(len(station_list)):
            
            stnp = StationPick()

            stnp.group = int(station_list[i][0])
            stnp.network = station_list[i][1]
            stnp.code = station_list[i][2]
            stnp.position = Position(np.degrees(float(station_list[i][3])), \
                                     np.degrees(float(station_list[i][4])), \
                                     float(station_list[i][5]))
            stnp.position.pos_loc(Position(self.setup.lat_centre, self.setup.lon_centre, 0))
            stnp.time = float(station_list[i][7])


            station_obj_list.append([stnp.group, stnp.network, stnp.code, \
                                    stnp.position.x, stnp.position.y, stnp.position.z, \
                                    stnp.time])

        sounding = parseWeather(self.setup)

        self.psoTrajectory(station_obj_list, sounding)

    def rayTrace(self):

        A = Position(float(self.ray_lat_edits.text()), float(self.ray_lon_edits.text()), float(self.ray_height_edits.text()))
        B = Position(self.ray_pick_point[0], self.ray_pick_point[1], self.ray_pick_point[2])

        A.pos_loc(B)
        B.pos_loc(B)

        try:
            sounding = parseWeather(self.setup)
        except:
            errorMessage('Error reading weather profile in rayTrace', 2)
            return None

        if self.prefs.debug:
            print("Starting and End points of Ray Trace")
            print(A)
            print(B)

        if self.setup.perturb_times == 0:
            self.setup.perturb_times = 1

        trace_data =    [None]*self.setup.perturb_times
        trace_var =     [None]*self.setup.perturb_times
        t_arrival =     [None]*self.setup.perturb_times
        t_arrival_cy =  [None]*self.setup.perturb_times
        err =           [None]*self.setup.perturb_times

        #plt.style.use('dark_background')
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(5, 5)
        ax = fig.add_subplot(1, 1, 1, projection='3d')

        if self.setup.perturb_method == 'ensemble':
            ensemble_file = self.setup.perturbation_spread_file
        else:
            ensemble_file = ''

        x_var = []
        y_var = []
        z_var = []
        p_var = []
        t_var = []
        error_list = []
        for ptb_n in range(self.setup.perturb_times):
            trace_data = []
            trace_var = []
            if ptb_n > 0 and self.ray_enable_perts.isChecked():
                
                if self.prefs.debug:
                    print("STATUS: Perturbation {:}".format(ptb_n))

                # generate a perturbed sounding profile
                sounding_p = perturb(self.setup, sounding, self.setup.perturb_method, \
                    spread_file=self.setup.perturbation_spread_file, lat=self.setup.lat_centre, lon=self.setup.lon_centre, ensemble_file=ensemble_file, ensemble_no=ptb_n)
            else:

                # if not using perturbations on this current step, then return the original sounding profile
                sounding_p = sounding


            z_profile, _ = getWeather(np.array([A.x, A.y, A.z]), np.array([B.x, B.y, B.z]), \
                    self.setup.weather_type, A, copy.copy(sounding_p))

            z_profile = zInterp(B.z, A.z, z_profile, div=100)

            a, b, c, E, trace_data = slowscan(A.xyz, B.xyz, z_profile, wind=True, n_theta=self.setup.n_theta, n_phi=self.setup.n_theta, h_tol=self.setup.h_tol, v_tol=self.setup.v_tol)

            if trace_data == trace_data:
                if self.ray_enable_vars.isChecked():

                    last_k = 0
                    N = 15

 
                    m, n = np.shape(trace_var[0][0])

                    for i in range(m//N):
                        for j in range(n//N):
                            for line in trace_var:
                                k = line[3]

                                if k != last_k:
                                    #c = (0, 0, (t_var[0] - np.pi/2)/np.pi/2%1)

                                    ax.plot3D(x_var, y_var, z_var, c='r')
                                    x_var = []
                                    y_var = []
                                    z_var = []
                                    p_var = []
                                    t_var = []
                                x_var.append(line[0][i*N, j*N])
                                y_var.append(line[1][i*N, j*N])
                                z_var.append(line[2][i*N, j*N])
                                p_var.append(line[4][i*N, j*N])
                                t_var.append(line[5][i*N, j*N])

                                last_k = k
                    ax.plot3D(x_var, y_var, z_var, c='r')                
                if ptb_n == 0:
                    xline = []
                    yline = []
                    zline = []

                    try:
                        for line in trace_data:
                            #line[0], line[1], line[2] = loc2Geo(A.lat, A.lon, A.elev, [line[0], line[1], line[2]])

                            xline.append(line[0])
                            yline.append(line[1])
                            zline.append(line[2])

                        ax.plot3D(np.array(xline)/1000, np.array(yline)/1000, np.array(zline)/1000, 'black')
                        #ax.scatter(xline, yline, zline, 'blue', marker='o')
                        #ax.scatter(0, 0, 0, 'orange', marker='^')
                    except IndexError:
                        pass
                    except TypeError:
                        pass
      

                    # ax.set_xlim3d(B.x, A.x)
                    # ax.set_ylim3d(B.y, A.y)
                    # ax.set_zlim3d(B.z, A.z)

                    x_pts = [None]*len(xline)
                    y_pts = [None]*len(xline)
                    for i in range(len(xline)):
                        
                        x_pts[i], y_pts[i], _ = loc2Geo(B.lat, B.lon, B.elev, [xline[i], yline[i], zline[i]])

                    self.ray_canvas.plot(y_pts, x_pts, pen=(255, 255, 255), update=True)

                if self.ray_enable_perts.isChecked():
                    xline = []
                    yline = []
                    zline = []
                    if ptb_n > 0:

                        for line in trace_data:
                            #line[0], line[1], line[2] = loc2Geo(A.lat, A.lon, A.elev, [line[0], line[1], line[2]])

                            xline.append(line[0])
                            yline.append(line[1])
                            zline.append(line[2])
                        try:
                            ax.plot3D(np.array(xline)/1000, np.array(yline)/1000, np.array(zline)/1000)#'#15ff00')
                        except:
                            pass
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

                    if self.setup.weather_type == 'custom':
                        u = mags*np.sin(np.radians(dirs))
                        v = mags*np.cos(np.radians(dirs))
                    else:
                        u = mags*np.sin(dirs)
                        v = mags*np.cos(dirs)

                    c = t[:-1]
                    c = (c.ravel() - c.min()) / c.ptp()
                    c = np.concatenate((c, np.repeat(c, 2)))
                    c = plt.cm.seismic(c)

                    xline = []
                    yline = []
                    zline = []

                    max_mag = np.nanmax(mags)/1000


                    for ii, line in enumerate(trace_data):
                        #line[0], line[1], line[2] = loc2Geo(A.lat, A.lon, A.elev, [line[0], line[1], line[2]])

                        xline = line[0]
                        yline = line[1]
                        zline = line[2]

                        try:
                            ax.quiver(np.array(xline)/1000, np.array(yline)/1000, np.array(zline)/1000, u[ii]/max_mag, v[ii]/max_mag, 0, color=c)
                        except:
                            pass

        
        avg_error = np.mean(error_list)
        print("Mean error in computation from loss in speed: {:5.2f}".format(avg_error))

        # ax.set_xlim3d(B.x, A.x)
        # ax.set_ylim3d(B.y, A.y)
        # ax.set_zlim3d(B.z, A.z)

        ax.set_xlabel('x (km +East)')
        ax.set_ylabel('y (km +North)')
        ax.set_zlabel('z (km +Up)')

        ax.set_title('Ray Trace from {:} km to \n {:}'.format(A.z/1000, B))

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

    def trajSolver(self):

        A = Position(self.setup.lat_i, self.setup.lon_i, self.setup.elev_i)
        B = Position(self.setup.lat_f, self.setup.lon_f, self.setup.elev_f)

        if A.isNone() or B.isNone():
            errorMessage('Positions of the trajectory are None!', 2, detail='Please define both ends of the trajectory!')
            return None

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

        if self.prefs.debug:
            print("Solved Point:")
            print(pt)

        self.ray_lat_edits.setText(str(pt.lat))
        self.ray_lon_edits.setText(str(pt.lon))

        self.ray_pick_traj.setPoints(x=[pt.lon], y=[pt.lat], pen=(255, 0, 110))
        self.ray_canvas.addItem(self.ray_pick_traj, update=True)

        if self.ray_pick_point != [0, 0, 0]:
            self.rayTrace()

    def W_estGUI(self):
        """ Opens up yield estimater GUI
        """

        try:
            self.w = Yield(self.bam, self.prefs, self.current_station)
        except AttributeError as e:
            errorMessage('Not enough data for yield generator', 2, detail='{:}'.format(e))

        self.w.setGeometry(QRect(100, 100, 800, 200))
        
        self.w.show()

    def showContour(self, mode):

        print('Working on contour - This could take a while...')
        self.clearContour()

        ref_pos = Position(45.9789390626, 14.9929810526, 0)


        ### option to use a perturbation for the contour instead of nominal (change the 0 to the perturbation number)
        # sounding = self.perturbGenerate(0, sounding, self.perturbSetup())

        if mode == 'ballistic':
            if errorCodes(self.bam.setup, 'trajectory'):
                return None
            try:    
                points = self.bam.setup.trajectory.findPoints(gridspace=100, min_p=17000, max_p=50000)
            except AttributeError as e:
                errorMessage('Trajectory is not defined!', 2, detail='{:}'.format(e))
                return None
        elif mode == 'fragmentation':
            if errorCodes(self.bam.setup, 'fragmentation_point'):
                return None
            try:
                A = self.setup.fragmentation_point[0].position
                A.pos_loc(ref_pos)
                points = A.xyz
            except (TypeError, IndexError) as e:
                errorMessage('Fragmentation Point not defined correctly!', 1, info='Please define the fragmentation point in the setup toolbar', detail='{:}'.format(e))
                return None


        results = waveReleasePointWindsContour(self.bam, ref_pos, points, mode=mode)

        results = np.array(results)

        dx, dy = 0.01, 0.01
        
        X = results[:, 0]
        Y = results[:, 1]
        T = results[:, 3]
        data = []
        for i in range(len(X)):
            A = Position(0, 0, 0)
            A.x = X[i]
            A.y = Y[i]
            A.z = 0
            A.pos_geo(ref_pos)
            data.append((A.lon, A.lat, dy, dx, T[i]))
            # return data in a form readable by Rectangle Object

        self.contour_data_squares = RectangleItem(data)
        self.make_picks_map_graph_canvas.addItem(self.contour_data_squares)
        print('Contour Finished!')

    def clearContour(self):

        self.make_picks_map_graph_canvas.removeItem(self.contour_data_squares)

    def saveContour(self):
        filename = QFileDialog.getSaveFileName(self, 'Save File')

        np.save(filename[0], self.contour_data)

    def loadContour(self):
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setNameFilters(['Contour Numpy File (*.npy)'])

        dlg.exec_()

        filename = dlg.selectedFiles()

        try:
            self.contour_data = np.load(filename[0])
        except IsADirectoryError as e:
            errorMessage('Directory chosen when .npy file expected!', 2, detail='{:}'.format(e))
            return None
        print(self.contour_data)
        try:
            self.make_picks_map_graph_canvas.addItem(RectangleItem(self.contour_data))
        except IndexError as e:
            errorMessage('The selected contour file is of the wrong size!', 2, info='Are you sure that this is a contour file and not some other .npy file?', detail='{:}'.format(e))

    def checkForWorkDir(self):
        try:
            if not os.path.exists(self.prefs.workdir):
                os.makedirs(self.prefs.workdir)
            return True
        except (FileNotFoundError, TypeError) as e:
            errorMessage("No such file or directory: '{:}'".format(self.prefs.workdir), 2, info='Define a working directory in the toolbar on the side.', detail='{:}'.format(e))
            return False

    def loadRayGraph(self):

        if not self.checkForWorkDir():
            return None

        #Build seismic data path
        dir_path = os.path.join(self.prefs.workdir, self.setup.fireball_name)

        # Load the station and waveform files list
        data_file_path = os.path.join(dir_path, DATA_FILE)

        if os.path.isfile(data_file_path):
            
            stn_list = readStationAndWaveformsListFile(data_file_path, rm_stat=self.setup.rm_stat, debug=self.prefs.debug)

        else:
            errorMessage('Station and waveform data file not found! Download the waveform files first!', 2)
            return None

        stn_list = stn_list + self.setup.stations

        for stn in stn_list:
            if stn.code not in self.setup.rm_stat:
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

        if self.tog_picks.isChecked():

            mousePoint = self.ray_canvas.vb.mapToView(evt.pos())
            
            self.ray_pick.setPoints(x=[mousePoint.x()], y=[mousePoint.y()], pen=(255, 0, 110))
            self.ray_canvas.addItem(self.ray_pick, update=True)
            self.ray_pick_point = [mousePoint.y(), mousePoint.x(), 0]
            self.ray_pick_label.setText("Lat: {:10.4f} Lon: {:10.4f} Elev {:10.2f}".format(*self.ray_pick_point))


    def effectiveSoundSpeed(self, sounding):
        '''
        Returns the sound speed at every pressure level, also considering the winds and the k-vectors
        '''

        #############################
        # Get start and end positions
        #############################
        lat = [tryFloat(self.fatm_start_lat.text()), tryFloat(self.fatm_end_lat.text())]
        lon = [tryFloat(self.fatm_start_lon.text()), tryFloat(self.fatm_end_lon.text())]
        elev = [tryFloat(self.fatm_start_elev.text()), tryFloat(self.fatm_end_elev.text())]

        supra_pos = Position(lat[0], lon[0], elev[0])
        detec_pos = Position(lat[1], lon[1], elev[1])
        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

        supra_pos.pos_loc(ref_pos)
        detec_pos.pos_loc(ref_pos)

        ################################
        # Get eigen-path from ray-tracer
        ################################
        pts = cyscanV(supra_pos.xyz, detec_pos.xyz, sounding, \
                    wind=self.prefs.wind_en, n_theta=self.prefs.pso_theta, n_phi=self.prefs.pso_phi,
                    h_tol=self.prefs.pso_min_ang, v_tol=self.prefs.pso_min_dist)

        # The path taken in an isotropic atmosphere - straight line
        u = supra_pos.xyz - detec_pos.xyz
        nom_range = np.sqrt(u[0]**2 + u[1]**2 + u[2]**2)

        ray_range = 0
        for ii in range(len(pts[0]) - 1):
            k = np.array([pts[0][ii + 1] - pts[0][ii],\
                          pts[1][ii + 1] - pts[1][ii],\
                          pts[2][ii + 1] - pts[2][ii]])
            ray_range += np.sqrt(k[0]**2 + k[1]**2 + k[2]**2)


        for ii in range(len(pts[0]) - 1):
            k = np.array([pts[0][ii + 1] - pts[0][ii],\
                          pts[1][ii + 1] - pts[1][ii],\
                          pts[2][ii + 1] - pts[2][ii]])
            k /= np.sqrt(k[0]**2 + k[1]**2 + k[2]**2)

        c = sounding[:, 1]
        mags = sounding[:, 2]
        dirs = sounding[:, 3]
        u = mags*np.sin(dirs)
        v = mags*np.cos(dirs)
        w = np.array([u, v, 0])

        c_eff = c + np.dot(k, w)

        return c_eff, nom_range, ray_range
            

    def fatmPlot(self, sounding, perturbations):
        
        self.fatm_canvas.clear()
        self.fatm_canvas.setLabel('left', "Height", units='m', size='24pt')
        if self.fatm_variable_combo.currentText() == 'Sound Speed':
            X = sounding[:, 1]
            self.fatm_canvas.setLabel('bottom', "Sound Speed", units='m/s', size='24pt')
        elif self.fatm_variable_combo.currentText() == 'Wind Magnitude':
            X = sounding[:, 2]
            self.fatm_canvas.setLabel('bottom', "Wind Magnitude", units='m/s', size='24pt')
        elif self.fatm_variable_combo.currentText() == 'Wind Direction':
            X = np.degrees(sounding[:, 3])
            self.fatm_canvas.addItem(pg.InfiniteLine(pos=(310, 0), angle=90, pen=QColor(255, 0, 0)))
            self.fatm_canvas.setLabel('bottom', "Wind Direction", units='deg from N', size='24pt')
        elif self.fatm_variable_combo.currentText() == 'Effective Sound Speed':
            X, nom_range, nom_ray_range = self.effectiveSoundSpeed(sounding)
            self.fatm_canvas.setLabel('bottom', "Sound Speed", units='m/s', size='24pt')    
        # elif self.fatm_variable_combo.currentText() == 'U-Component of Wind':
        #     dirs = angle2NDE(np.degrees(sounding[:, 3]))
        #     mags = sounding[:, 2]
        #     X = mags*np.cos(np.radians(dirs))
        # elif self.fatm_variable_combo.currentText() == 'V-Component of Wind':
        #     dirs = angle2NDE(np.degrees(sounding[:, 3]))
        #     mags = sounding[:, 2]
        #     X = mags*np.sin(np.radians(dirs))
        else:
            return None
        Y = sounding[:, 0]

        
        self.fatm_canvas.plot(x=X, y=Y, pen='w')
        SolutionGUI.update(self)
        perts_range = []
        if len(perturbations) != 0:
            for ii, ptb in enumerate(perturbations):
                
                if self.fatm_variable_combo.currentText() == 'Sound Speed':
                    X = ptb[:, 1]
                elif self.fatm_variable_combo.currentText() == 'Wind Magnitude':
                    X = ptb[:, 2]
                elif self.fatm_variable_combo.currentText() == 'Wind Direction':
                    X = np.degrees(ptb[:, 3])
                elif self.fatm_variable_combo.currentText() == 'Effective Sound Speed':
                    X, nom_range, pert_ray_range = self.effectiveSoundSpeed(ptb)
                    perts_range.append(pert_ray_range)
                # elif self.fatm_variable_combo.currentText() == 'U-Component of Wind':
                #     dirs = angle2NDE(np.degrees(ptb[:, 3]))
                #     mags = ptb[:, 2]
                #     X = mags*np.cos(np.radians(dirs))
                # elif self.fatm_variable_combo.currentText() == 'V-Component of Wind':
                #     dirs = angle2NDE(np.degrees(ptb[:, 3]))
                #     mags = ptb[:, 2]
                #     X = mags*np.sin(np.radians(dirs))
                else:
                    return None
                Y = ptb[:, 0]
                self.fatm_canvas.plot(x=X, y=Y, pen='g')
                SolutionGUI.update(self)
        perts_range = np.array(perts_range)
        try:
            print('Isotropic Range:           {:.2f}         km'.format(nom_range/1000))
            print('Nominal Atmospheric Range: {:.2f}         km'.format(nom_ray_range/1000))
            print('Perturbation Range:        {:.2f} - {:.2f} km'.format(np.nanmin(perts_range/1000), np.nanmax(perts_range/1000)))
            print('Average Speed:             {:.2f}         m/s'.format(np.nanmean(X)))
        except:
            pass
        font=QtGui.QFont()
        font.setPixelSize(20)
        self.fatm_canvas.getAxis("bottom").tickFont = font
        self.fatm_canvas.getAxis("left").tickFont = font
        self.fatm_canvas.getAxis('bottom').setPen(self.color.WHITE) 
        self.fatm_canvas.getAxis('left').setPen(self.color.WHITE)
        # self.fatm_canvas.getLabel("bottom").tickFont = font
        # self.fatm_canvas.getLabel("left").tickFont = font
        # self.fatm_canvas.getLabel('bottom').setPen(self.color.WHITE) 
        # self.fatm_canvas.getLabel('left').setPen(self.color.WHITE)

        SolutionGUI.update(self)

    def fatmSaveAtm(self):
        
        if self.fatm_source_type.currentText() == "Copernicus Climate Change Service (ECMWF)":
            weather_type = 'ecmwf'
        elif self.fatm_source_type.currentText() == "Copernicus Climate Change Service (ECMWF) - Spread File":
            weather_type = 'spread'
        elif self.fatm_source_type.currentText() == "Radiosonde":
            weather_type = 'radio'

        self.bam.atmos.loadSounding(self.fatm_name_edits.text(), weather_type, \
            lat=self.bam.setup.lat_centre, lon=self.bam.setup.lon_centre, \
            rng=self.bam.setup.deg_radius, time=self.fatm_datetime_edits.dateTime())

        save(self)

    def fatmLoadAtm(self):

        lat = [tryFloat(self.fatm_start_lat.text()), tryFloat(self.fatm_end_lat.text())]
        lon = [tryFloat(self.fatm_start_lon.text()), tryFloat(self.fatm_end_lon.text())]
        elev = [tryFloat(self.fatm_start_elev.text()), tryFloat(self.fatm_end_elev.text())]

        self.fatmPlot(*self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev))


    def fatmFetch(self, download):

        if self.fatm_source_type.currentIndex() == 0:

            year = str(self.fatm_datetime_edits.dateTime().date().year())
            month = str(self.fatm_datetime_edits.dateTime().date().month())
            day = str(self.fatm_datetime_edits.dateTime().date().day())

            time_of = str("{:02d}".format(self.fatm_datetime_edits.dateTime().time().hour())) + \
            ':' + str("{:02d}".format(self.fatm_datetime_edits.dateTime().time().minute()))

            loc = self.fatm_name_edits.text()

            variables = []
            self.fatm_variable_combo.clear() 
            if self.fatm_temp.isChecked():
                variables.append('temperature')
            if self.fatm_u_wind.isChecked():
                variables.append('u_component_of_wind')
            if self.fatm_v_wind.isChecked():
                variables.append('v_component_of_wind')

            if download:
            
                if self.fatm_perts.isChecked():
                    print("Downloading Perturbation Ensemble")
                else:
                    print("Downloading Reanalysis")

                try:
                    copernicusAPI(variables, year, month, day, time_of, loc, ensemble=self.fatm_perts.isChecked())
                except Exception as e:
                    errorMessage("Error downloading weather data from CDS", 1, detail='{:}'.format(e))

                return None
            else:
                self.fatm_variable_combo.clear() 
                if self.fatm_temp.isChecked():
                    self.fatm_variable_combo.addItem('Temperature')
                if self.fatm_u_wind.isChecked():
                    self.fatm_variable_combo.addItem('U-Component of Wind')
                if self.fatm_v_wind.isChecked():
                    self.fatm_variable_combo.addItem('V-Component of Wind')
                if self.fatm_u_wind.isChecked() and self.fatm_v_wind.isChecked():
                    self.fatm_variable_combo.addItem('Wind Magnitude')
                    self.fatm_variable_combo.addItem('Wind Direction')
                self.fatmPlot()

        elif self.fatm_source_type.currentIndex() == 1:

            loc = self.fatm_name_edits.text()

            if download:
                
                if self.prefs.debug:
                    print("Looking for stations near: Lat: {:.4f} Lon: {:.4f}".format(self.setup.lat_centre, self.setup.lon_centre))
                
                downloadRadio(self.setup.lat_centre, self.setup.lon_centre, loc, debug=self.prefs.debug)

                with zipfile.ZipFile(loc, 'r') as zip_ref:
                    names = zip_ref.namelist()
                    zip_ref.extractall(os.path.join(self.prefs.workdir, self.setup.fireball_name))

                if self.prefs.debug:
                    print('File extracted {:}'.format(names[0]))
            else:
                self.fatm_variable_combo.addItem('Temperature')
                self.fatm_variable_combo.addItem('U-Component of Wind')
                self.fatm_variable_combo.addItem('V-Component of Wind')
                self.fatm_variable_combo.addItem('Wind Magnitude')
                self.fatm_variable_combo.addItem('Wind Direction')
                self.fatmPlot()

    def fatmPrint(self):
        
        filename = QFileDialog.getSaveFileName(self, 'Save File', '', 'Text File (*.txt)')

        self.setup.lat_centre = tryFloat(self.fatm_end_lat.text())
        self.setup.lon_centre = tryFloat(self.fatm_end_lon.text())
        self.setup.sounding_file = self.fatm_name_edits.text()
        atm_time = self.fatm_datetime_edits.dateTime().time().hour()

        variables = []

        if self.setup.weather_type == 'ecmwf':
            if self.fatm_temp.isChecked():
                variables.append('t')
            if self.fatm_u_wind.isChecked():
                variables.append('u')
            if self.fatm_v_wind.isChecked():
                variables.append('v')
        elif self.setup.weather_type == 'merra':
            if self.fatm_temp.isChecked():
                variables.append('T')
            if self.fatm_u_wind.isChecked():
                variables.append('U')
            if self.fatm_v_wind.isChecked():
                variables.append('V')

        if errorCodes(self.setup, 'sounding_file'):
            return None

        try:
            float(self.setup.lat_centre)
            float(self.setup.lon_centre)
        except TypeError as e:
            errorMessage('Lat centre and/or lon centre are not floats or are not defined!', 2, detail='{:}'.format(e))
            return None

        sounding = parseGeneralECMWF(self.setup.weather_type, self.setup.sounding_file, self.setup.lat_centre, \
                        self.setup.lon_centre, atm_time, variables)

        filename = checkExt(filename[0], '.txt')

        header = 'Pressures (hPa)'

        for element in variables:

            if element == variables[-1]:
                header = header + ', ' + element + '\n'
            else:
                header = header + ', ' + element


        with open(str(filename), 'w') as f:
            
            f.write(header)

            for line in sounding:
                info = ''
                for element in line:
                    if element == line[-1]:
                        info = info + str(element) + '\n'
                    else:
                        info = info + str(element) + ','
                f.write(info)

        if self.fatm_perts.isChecked():
            dataset = parseWeather(self.setup)
            for ptb_n in range(self.setup.perturb_times):
                sounding_p = self.perturbGenerate(ptb_n, dataset, self.perturbSetup())
                zProfile, _ = getWeather(np.array([self.setup.lat_centre, self.setup.lon_centre, 73721.58]), \
                                         np.array([self.setup.lat_centre, self.setup.lon_centre, 100]), \
                                         self.setup.weather_type, \
                                        [self.setup.lat_centre, self.setup.lon_centre, 100], \
                                        copy.copy(sounding_p), convert=False)
                zProfile = zInterp(100, 73721.58, zProfile, div=100)
                zProfile_ext = []
                for line in zProfile:
                    zProfile_ext.append((line[0]/1000, line[1], line[2], line[3], estPressure(line[0])/100))
                    
                    with open(str(filename) + str(ptb_n), 'w+') as f:
            
                        f.write(header)

                        for line in zProfile_ext:
                            info = ''
                            for element in line:
                                if element == line[-1]:
                                    info = info + str(element) + '\n'
                                else:
                                    info = info + str(element) + ','
                            f.write(info)


        errorMessage('Printed out sounding data', 0, title="Print Done")

    def refLoadTrajectory(self):
        self.latedit.setText(str(self.setup.trajectory.pos_f.lat))
        self.lonedit.setText(str(self.setup.trajectory.pos_f.lon))
        self.azedit.setText(str(self.setup.trajectory.azimuth.deg))
        self.zeedit.setText(str(self.setup.trajectory.zenith.deg))
        self.veedit.setText(str(self.setup.trajectory.v))
        self.vfedit.setText(str(self.setup.trajectory.v_f))
        self.tiedit.setText(str(self.setup.trajectory.t))

    def refLoad(self):
        """ Loads csv file into a table
        """

        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setNameFilters(['CSV File (*.csv)'])
        dlg.exec_()

        filename = dlg.selectedFiles()

        try:
            with open(filename[0]) as f:
                data_table = []
                next(f)
        
                for line in f:
                    a = line.split(',')
                    if len(a) != 5:
                        errorMessage('Wrong number of columns for file!', 1, info='Make sure the right file is imported!')
                        return None
                    data_table.append(a)


        except IsADirectoryError as e:
            errorMessage('Please select a valid file to load', 1, detail='{:}'.format(e))
            return None

        defTable(self.ref_table, 0, 5, headers=['Height', 'Latitude', 'Longitude', 'Time', 'δ Height'])

        toTable(self.ref_table, data_table)


    def refSave(self):
        """  Saves a table to a csv
        """

        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = checkExt(dlg[0], '.csv')

        data_set = fromTable(self.ref_table)
        # Open the output CSV
        with open(os.path.join(file_name), 'w') as f:

            # Write the header
            f.write('Height Latitude Longitude Time \n')

            # Go through all picks
            for line in data_set:
                line[-1] = int(line[-1])
                # Write the CSV entry
                f.write("{:.4f}, {:.4f}, {:.4f}, {:.4f}, {:.4f} \n".format(*line))

        errorMessage('Output to CSV!', 0, title='Exported!', detail='Filename: {:}'.format(file_name))


    def refHighStats(self):
        for ii, stn in enumerate(self.stn_list):

            if stn.code in self.setup.high_f and stn.code in self.setup.high_b:
                self.ref_stn_view[ii].setBackground(self.color.both)
            elif stn.code in self.setup.high_b:
                self.ref_stn_view[ii].setBackground(self.color.ballistic)
            elif stn.code in self.setup.high_f:
                self.ref_stn_view[ii].setBackground((51, 153, 51))

    def refClearStations(self):
        # Set a blank widget to remove all stations
        widget = QWidget()
        self.ref_waveforms.setWidget(widget)

    def refBallStations(self):
        self.refClearStations()
        self.refLoadStations(ballonly=True)

    def refFragStations(self):
        self.refClearStations()
        self.refLoadStations(fragonly=True)

    def refBothStations(self):
        self.refClearStations()
        self.refLoadStations(ballonly=True, fragonly=True)

    def refLoadStations(self, ballonly=False, fragonly=False):
        if not self.checkForWorkDir():
            return None

            #Build seismic data path
        dir_path = os.path.join(self.prefs.workdir, self.setup.fireball_name)

        # Load the station and waveform files list
        data_file_path = os.path.join(dir_path, DATA_FILE)
        self.stn_list = []
        if os.path.isfile(data_file_path):
            
            self.stn_list = readStationAndWaveformsListFile(data_file_path, rm_stat=self.setup.rm_stat)
            self.stn_list = self.stn_list + self.setup.stations
        else:
            errorMessage('Station and waveform data file not found! Download the waveform files first!', 2)
            sys.exit()

        # Need to check for duplicates here
        stn_list = self.stn_list

        widget = QWidget()
        widget.setStyleSheet('background-color: black;')
        layout = QVBoxLayout(widget)
        layout.setAlignment(Qt.AlignTop)

        # blankwidget = QWidget()
        # blankwidget.setStyleSheet('background-color: rgb(0, 100, 200);')

        self.ref_stn_view = []
        self.ref_stn_canvas = []
        self.waveform_data = [None]*len(stn_list)

        for index in range(len(stn_list)):
            
            stn = stn_list[index]

            if ballonly:
                if stn.code not in self.setup.high_b:
                    self.ref_stn_view.append(None)
                    self.ref_stn_canvas.append(None)
                    continue

            if fragonly:
                if stn.code not in self.setup.high_f:
                    self.ref_stn_view.append(None)
                    self.ref_stn_canvas.append(None)
                    continue

            stn = stn_list[index]
            station_layout = QVBoxLayout()
            layout.addLayout(station_layout)
            label_layout = QVBoxLayout()
            waveform_layout = QGridLayout()
            station_layout.addLayout(label_layout)
            station_layout.addLayout(waveform_layout)
            # station_layout.addWidget(blankwidget, 10, 0, 1, 20)

            label_layout.addWidget(QLabel('--- Station: {:2}-{:5} ---'.format(stn_list[index].network, stn_list[index].code)))

            self.ref_stn_view.append(pg.GraphicsLayoutWidget())
            self.ref_stn_canvas.append(self.ref_stn_view[index].addPlot())
            self.ref_stn_canvas[index].getAxis('bottom').setPen((255, 255, 255)) 
            self.ref_stn_canvas[index].getAxis('left').setPen((255, 255, 255)) 
            waveform_layout.addWidget(self.ref_stn_view[index], index, 0)
            _, _ = self.discountDrawWaveform(self.setup, index, self.ref_stn_canvas[index])


        self.ref_waveforms.setWidget(widget)
        self.ref_waveforms.setWidgetResizable(True)


    def makeRefTraj(self):

        self.ref_traj_lat = float(self.latedit.text())
        self.ref_traj_lon = float(self.lonedit.text())
        self.ref_traj_az = float(self.azedit.text())
        self.ref_traj_ze = float(self.zeedit.text())
        self.ref_traj_t = float(self.tiedit.text())
        self.ref_traj_v = float(self.veedit.text())
        self.ref_traj_vf = float(self.vfedit.text())

        self.ref_traj = Trajectory(self.ref_traj_t, self.ref_traj_v, \
                                    zenith=Angle(self.ref_traj_ze), \
                                    azimuth=Angle(self.ref_traj_az), \
                                    pos_f=Position(self.ref_traj_lat, self.ref_traj_lon, 0), v_f=self.ref_traj_vf)

    def refSyncHeights(self):

        self.makeRefTraj()

        frags = fromTable(self.ref_table)

        new_frags = []

        for f in frags:
            P = self.ref_traj.findGeo(f[0])
            T = self.ref_traj.findTime(f[0])

            new_frags.append([P.elev, P.lat, P.lon, T, 0])

        toTable(self.ref_table, new_frags)


    def refPlotTraj(self):
        
        self.makeRefTraj()

        ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)

        sounding = parseWeather(self.setup)

        self.ref_traj.pos_f.pos_loc(ref_pos)

        points = self.ref_traj.findPoints(gridspace=100, min_p=17000, max_p=50000)

        for ii, stn in enumerate(self.stn_list):

            stn.position.pos_loc(ref_pos)

            b_time = timeOfArrival(stn.position.xyz, self.ref_traj, self.setup, points, sounding=sounding, \
                            travel=False, fast=False, ref_loc=ref_pos, theo=False, div=37)


            try:
                self.ref_stn_canvas[ii].addItem(pg.InfiniteLine(pos=(b_time, 0), angle=90, pen=QColor(0, 0, 255)))
                # self.setColortoRow(self.ref_table, jj, QColor("blue"))

            except:
                # no f_time
                pass

        SolutionGUI.update(self)

    def refPlotHeights(self):

        alpha = 100

        colors = [QColor(0, 255, 26), QColor(3, 252, 176), QColor(252, 3, 3), QColor(176, 252, 3), QColor(255, 133, 3),
                QColor(149, 0, 255), QColor(76, 128, 4), QColor(82, 27, 27), QColor(101, 128, 125), QColor(5, 176, 249)]

        colors_alpha = [QColor(0, 255, 26, alpha), QColor(3, 252, 176, alpha), QColor(252, 3, 3, alpha), QColor(176, 252, 3, alpha), QColor(255, 133, 3, alpha),
                QColor(149, 0, 255, alpha), QColor(76, 128, 4, alpha), QColor(82, 27, 27, alpha), QColor(101, 128, 125, alpha), QColor(5, 176, 249, alpha)]

        self.makeRefTraj()
        frags = fromTable(self.ref_table)

        ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)

        sounding = parseWeather(self.setup)

        for jj, f in enumerate(frags):
            S = Position(f[1], f[2], f[0])
            S.pos_loc(ref_pos)
            T = f[3]
            U = f[4]
            

            if U == 0:
                for ii, stn in enumerate(self.stn_list):
                    stn.position.pos_loc(ref_pos)

                    # Cut down atmospheric profile to the correct heights, and interp
                    zProfile, _ = getWeather(np.array([S.lat, S.lon, S.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), self.setup.weather_type, \
                            [ref_pos.lat, ref_pos.lon, ref_pos.elev], copy.copy(sounding), convert=False)
                    # Travel time of the fragmentation wave
                    f_time, frag_azimuth, frag_takeoff, frag_err = cyscan(S.xyz, stn.position.xyz, zProfile, wind=True, \
                        n_theta=self.setup.n_theta, n_phi=self.setup.n_phi, h_tol=self.setup.h_tol, v_tol=self.setup.v_tol)

                    try:
                        self.ref_stn_canvas[ii].addItem(pg.InfiniteLine(pos=(f_time, 0), angle=90, pen=colors[jj]))
                        # self.setColortoRow(self.ref_table, jj, QColor("blue"))
                        for j in range(5):
                            self.ref_table.item(jj, j).setBackground(colors[jj])
                    except:
                        # no f_time
                        pass
            else:
                P_upper = self.ref_traj.findGeo(f[0] + U)
                T_upper = self.ref_traj.findTime(f[0] + U)

                P_lower = self.ref_traj.findGeo(f[0] - U)
                T_lower = self.ref_traj.findTime(f[0] - U)


                for ii, stn in enumerate(self.stn_list):

                    stn.position.pos_loc(ref_pos)

                    # Cut down atmospheric profile to the correct heights, and interp
                    zProfile_u, _ = getWeather(np.array([P_upper.lat, P_upper.lon, P_upper.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), self.setup.weather_type, \
                            [ref_pos.lat, ref_pos.lon, ref_pos.elev], copy.copy(sounding), convert=False)
                    
                    # Cut down atmospheric profile to the correct heights, and interp
                    zProfile_l, _ = getWeather(np.array([P_lower.lat, P_lower.lon, P_lower.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), self.setup.weather_type, \
                            [ref_pos.lat, ref_pos.lon, ref_pos.elev], copy.copy(sounding), convert=False) 
                    
                    # Travel time of the fragmentation wave
                    f_time_u, frag_azimuth, frag_takeoff, frag_err = cyscan(P_upper.xyz, stn.position.xyz, zProfile_u, wind=True, \
                        n_theta=self.setup.n_theta, n_phi=self.setup.n_phi, h_tol=self.setup.h_tol, v_tol=self.setup.v_tol)

                    f_time_l, frag_azimuth, frag_takeoff, frag_err = cyscan(P_lower.xyz, stn.position.xyz, zProfile_l, wind=True, \
                        n_theta=self.setup.n_theta, n_phi=self.setup.n_phi, h_tol=self.setup.h_tol, v_tol=self.setup.v_tol)

                    try:

                        self.ref_stn_canvas[ii].addItem(pg.LinearRegionItem(values=(f_time_l, f_time_u), pen=colors[jj], brush=colors_alpha[jj], movable=False))

                        # self.ref_stn_canvas[ii].addItem(pg.InfiniteLine(pos=(f_time, 0), angle=90, pen=colors[jj]))
                        # self.setColortoRow(self.ref_table, jj, QColor("blue"))
                        for j in range(5):
                            self.ref_table.item(jj, j).setBackground(colors[jj])
                    except:
                        # no f_time
                        pass
        
        SolutionGUI.update(self)

    def discountDrawWaveform(self, setup, station_no, canvas):
        # Extract current station
        stn = self.stn_list[station_no]

        # Get the miniSEED file path
        mseed_file_path = os.path.join(self.prefs.workdir, setup.fireball_name, stn.file_name)

        # Try reading the mseed file, if it doesn't work, skip to the next frame
        try:
            mseed = obspy.read(mseed_file_path)

        except (TypeError, FileNotFoundError):

            print('mseed file could not be read:', mseed_file_path)
            return None, None


        # Unpact miniSEED data
        delta = mseed[0].stats.delta
        start_datetime = mseed[0].stats.starttime.datetime
        end_datetime = mseed[0].stats.endtime.datetime

        stn.offset = (start_datetime - setup.fireball_datetime).total_seconds()

        waveform_data = mseed[0].data

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

        # Get bandpass filter values
        bandpass_low = float(2)
        bandpass_high = float(8)


        # Init the butterworth bandpass filter
        butter_b, butter_a = butterworthBandpassFilter(bandpass_low, bandpass_high, \
            1.0/self.current_waveform_delta, order=6)

        # Filter the data
        waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(self.current_waveform_raw))

        # Plot the waveform
        self.waveform_data[station_no] = pg.PlotDataItem(x=time_data, y=waveform_data, pen='w')
        canvas.addItem(self.waveform_data[station_no])
        
        ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)

        try:
            t_arrival = stn.position.pos_distance(ref_pos)/self.prefs.avg_sp_sound
            canvas.setXRange(t_arrival-100, t_arrival+100, padding=1)

        except TypeError:
            t_arrival = None

        #canvas.setLabel('bottom', "Time after {:}".format(setup.fireball_datetime), units='s')
        # canvas.setLabel('left', "Signal Response")

        return np.min(waveform_data), np.max(waveform_data)

    def loadFrags(self):
        
        frag_list = []
        frag_name = fileSearch(['CSV (*.csv)'], None)
        try:
            with open(frag_name, 'r') as f:

                for line in f:
                    a = line.strip('\n').split(',')
                    a = [float(i) for i in a]

                    frag_list.append(a)

                f.close()
        except IsADirectoryError as e:
            errorMessage('Invalid CSV File!', 2, info='Directory given when .csv file was expected!', detail='{:}'.format(e))
            return None
        toTable(self.fragmentation_point, frag_list)


    def addIniDockWidgets(self):

        all_vars = QGroupBox()
        self.ini_dock.setWidget(all_vars)

        dock_layout = QVBoxLayout()
        all_vars.setLayout(dock_layout)

        ini_tabs = QTabWidget()
        dock_layout.addWidget(ini_tabs)

        self.load_button = QPushButton('Load')
        dock_layout.addWidget(self.load_button)
        #self.load_button.clicked.connect(partial(loadGUI, self))
        self.load_button.clicked.connect(partial(load, self))

        self.save_button = QPushButton('Save')
        dock_layout.addWidget(self.save_button)
        # self.save_button.clicked.connect(partial(saveGUI, self, True))
        self.save_button.clicked.connect(partial(save, self))

        tab1 = QWidget()
        tab1_content = QGridLayout()
        tab1.setLayout(tab1_content)
        ini_tabs.addTab(tab1, "General")

        self.fireball_name_label, self.fireball_name_edits = createLabelEditObj('Fireball Name:', tab1_content, 1, tool_tip='fireball_name')
        self.get_data_label, self.get_data_edits = createComboBoxObj('Get Data: ', tab1_content, 2, items=['False', 'True'], tool_tip='get_data')
        self.run_mode_label, self.run_mode_edits = createComboBoxObj('Run Mode: ', tab1_content, 3, items=['Search', 'Replot', 'Manual'], tool_tip='run_mode')
        self.arrival_times_label, self.arrival_times_edits, self.arrival_times_buton = createFileSearchObj('Arrival Times:', tab1_content, 6, width=1, h_shift=0, tool_tip='arrival_times_file')
        self.arrival_times_buton.clicked.connect(partial(fileSearch, ['Numpy Array (*.npy)'], self.arrival_times_edits))
        # self.sounding_file_label, self.sounding_file_edits, self.sounding_file_buton = createFileSearchObj('Sounding File:', tab1_content, 7, width=1, h_shift=0, tool_tip='sounding_file')
        # self.sounding_file_buton.clicked.connect(partial(fileSearch, ['NetCDF (*.nc)', 'HDF (*.HDF)', 'CSV (*.csv)', 'TXT (*.txt)'], self.sounding_file_edits))
        # self.perturbation_file_label, self.perturbation_file_edits, self.perturbation_file_buton = createFileSearchObj('Perturbation', tab1_content, 8, width=1, h_shift=0, tool_tip='perturbation_spread_file')
        # self.perturbation_file_buton.clicked.connect(partial(fileSearch, ['NetCDF (*.nc)'], self.perturbation_file_edits))
        self.station_picks_label, self.station_picks_edits, self.station_picks_buton = createFileSearchObj('Station Picks File: ', tab1_content, 9, width=1, h_shift=0, tool_tip='station_picks_file')
        self.station_picks_buton.clicked.connect(partial(fileSearch, ['CSV (*.csv)', 'Text File (*.txt)'], self.station_picks_edits))
        self.points_name_label, self.points_name_edits, self.points_name_buton = createFileSearchObj('Replot Points File: ', tab1_content, 10, width=1, h_shift=0, tool_tip='points_name')
        self.points_name_buton.clicked.connect(partial(fileSearch, ['CSV (*.csv)'], self.points_name_edits))
        self.lat_centre_label, self.lat_centre_edits = createLabelEditObj('Latitude Center:', tab1_content, 11, tool_tip='lat_centre', validate='float')
        self.lon_centre_label, self.lon_centre_edits = createLabelEditObj('Longitude Center:', tab1_content, 12, tool_tip='lon_centre', validate='float')
        self.deg_radius_label, self.deg_radius_edits = createLabelEditObj('Degrees in Search Radius:', tab1_content, 13, tool_tip='deg_radius', validate='float')
        self.fireball_datetime_label, self.fireball_datetime_edits = createLabelDateEditObj("Fireball Datetime", tab1_content, 14, tool_tip='fireball_datetime')
        self.t0_label, self.t0_edits = createLabelEditObj('t0:', tab1_content, 16, width=3, tool_tip='t0', validate='float')
        self.v_label, self.v_edits = createLabelEditObj('v:', tab1_content, 17, width=3, tool_tip='v', validate='float')
        self.azim_label, self.azim_edits = createLabelEditObj('azim:', tab1_content, 18, width=3, tool_tip='azim', validate='float')
        self.zangle_label, self.zangle_edits = createLabelEditObj('zangle:', tab1_content, 19, width=3, tool_tip='zangle', validate='float')
        self.lat_i_label, self.lat_i_edits = createLabelEditObj('lat_i:', tab1_content, 20, tool_tip='lat_i', validate='float')
        self.lon_i_label, self.lon_i_edits = createLabelEditObj('lon_i:', tab1_content, 21, tool_tip='lon_i', validate='float')
        self.elev_i_label, self.elev_i_edits = createLabelEditObj('elev_i:', tab1_content, 22, tool_tip='elev_i', validate='float')
        self.lat_f_label, self.lat_f_edits = createLabelEditObj('lat_f:', tab1_content, 23, tool_tip='lat_f', validate='float')
        self.lon_f_label, self.lon_f_edits = createLabelEditObj('lon_f:', tab1_content, 24, tool_tip='lon_f', validate='float')
        self.elev_f_label, self.elev_f_edits = createLabelEditObj('elev_f:', tab1_content, 25, tool_tip='elev_f', validate='float')
        self.vf_label, self.vf_edits = createLabelEditObj('v_f:', tab1_content, 27, width=3, validate='float')

        tab2 = QWidget()
        tab2_content = QGridLayout()
        tab2.setLayout(tab2_content)
        ini_tabs.addTab(tab2, "Sources")

        self.fragmentation_point_label = QLabel("Fragmentation Point(s)")
        self.fragmentation_point = QTableWidget(0, 4)
        tab2_content.addWidget(self.fragmentation_point, 1, 1, 1, 4)
        tab2_content.addWidget(self.fragmentation_point_label, 0, 1, 1, 2)
        self.fragmentation_csv = QPushButton("Load Points")
        tab2_content.addWidget(self.fragmentation_csv, 0, 3, 1, 1)
        self.fragmentation_csv.clicked.connect(self.loadFrags)
        self.fragmentation_point.setHorizontalHeaderLabels(['Latitude', 'Longitude', 'Elevation', 'Time'])
        header = self.fragmentation_point.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        self.fragmentation_point_label.setToolTip(toolTime('fragmentation_point'))
        self.fragmentation_point_add = QPushButton("+")
        tab2_content.addWidget(self.fragmentation_point_add, 2, 3, 1, 2)
        self.fragmentation_point_add.clicked.connect(partial(changeRows, self.fragmentation_point, 1))
        self.fragmentation_point_add.setToolTip("Add row")
        self.fragmentation_point_min = QPushButton("-")
        tab2_content.addWidget(self.fragmentation_point_min, 2, 1, 1, 2)
        self.fragmentation_point_min.clicked.connect(partial(changeRows, self.fragmentation_point, -1))
        self.fragmentation_point_min.setToolTip("Remove row")
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
        self.v_fixed_label, self.v_fixed_edits = createLabelEditObj('v_fixed:', tab2_content, 7, width=3, tool_tip='v_fixed', validate='float')
        self.restricted_time_check = QCheckBox("Enable Restricted Time: ")
        tab2_content.addWidget(self.restricted_time_check, 9, 4, 1, 1)
        self.restricted_time_label, self.restricted_time_edits = createLabelDateEditObj("Restricted Time: ", tab2_content, 9, width=2, tool_tip='restricted_time')
        self.azimuth_min_label, self.azimuth_min_edits = createLabelEditObj('azimuth_min:', tab2_content, 10, tool_tip='azimuth_min', validate='float')
        self.azimuth_max_label, self.azimuth_max_edits = createLabelEditObj('azimuth_max:', tab2_content, 10, h_shift=2, tool_tip='azimuth_max', validate='float')
        self.zangle_min_label, self.zangle_min_edits = createLabelEditObj('zangle_min:', tab2_content, 11, tool_tip='zenith_min', validate='float')
        self.zangle_max_label, self.zangle_max_edits = createLabelEditObj('zangle_max:', tab2_content, 11, h_shift=2, tool_tip='zenith_max', validate='float')
        self.lat_min_label, self.lat_min_edits = createLabelEditObj('lat_min', tab2_content, 12, tool_tip='x_min', validate='float')
        self.lat_max_label, self.lat_max_edits = createLabelEditObj('lat_max', tab2_content, 12, h_shift=2, tool_tip='x_max', validate='float')
        self.lon_min_label, self.lon_min_edits = createLabelEditObj('lon_min', tab2_content, 13, tool_tip='y_min', validate='float')
        self.lon_max_label, self.lon_max_edits = createLabelEditObj('lon_max', tab2_content, 13, h_shift=2, tool_tip='y_max', validate='float')
        self.elev_min_label, self.elev_min_edits = createLabelEditObj('elev_min:', tab2_content, 14, tool_tip='z_min', validate='float')
        self.elev_max_label, self.elev_max_edits = createLabelEditObj('elev_max:', tab2_content, 14, h_shift=2, tool_tip='z_max', validate='float')
        self.t_min_label, self.t_min_edits = createLabelEditObj('t_min:', tab2_content, 15, tool_tip='t_min', validate='float')
        self.t_max_label, self.t_max_edits = createLabelEditObj('t_max:', tab2_content, 15, h_shift=2, tool_tip='t_max', validate='float')
        self.v_min_label, self.v_min_edits = createLabelEditObj('v_min:', tab2_content, 16, tool_tip='v_min', validate='float')
        self.v_max_label, self.v_max_edits = createLabelEditObj('v_max:', tab2_content, 16, h_shift=2, tool_tip='v_max', validate='float')
        self.weight_distance_min_label, self.weight_distance_min_edits = createLabelEditObj('weight_distance_min', tab2_content, 17, tool_tip='weight_distance_min', validate='float')
        self.weight_distance_max_label, self.weight_distance_max_edits = createLabelEditObj('weight_distance_max', tab2_content, 17, h_shift=2, tool_tip='weight_distance_max', validate='float')

        tab3 = QWidget()
        tab3_content = QGridLayout()
        tab3.setLayout(tab3_content)
        ini_tabs.addTab(tab3, "Atmosphere")

        self.frag_no_label, self.frag_no_edits = createLabelEditObj('Fragmentation Number', tab3_content, 4, tool_tip='fragno', validate='int')

        tab4 = QWidget()
        tab4_content = QGridLayout()
        tab4.setLayout(tab4_content)
        ini_tabs.addTab(tab4, "Tweaks")

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
        self.extra_point_add.clicked.connect(partial(changeRows, self.extra_point, 1))
        self.extra_point_min = QPushButton("-")
        tab5_content.addWidget(self.extra_point_min, 2, 1, 1, 1)
        self.extra_point_min.clicked.connect(partial(changeRows, self.extra_point, -1))

        tab6 = QWidget()
        tab6_content = QGridLayout()
        tab6.setLayout(tab6_content)
        ini_tabs.addTab(tab6, "Notes")


    def atmPlotProfile(self, lat, lon, var_typ='t', perturb='none'):
        
        consts = Constants()

        if self.setup.weather_type == 'none':
            errorMessage('Weather type is set to "none", no weather can be displayed', 1)
            return None

        dataset = parseWeather(self.setup)

        if self.setup.weather_type == 'ecmwf':
            sounding = findECMWFSound(lat, lon, dataset)
        elif self.setup.weather_type == 'binary':
            sounding = findAus(lat, lon, dataset)
        elif self.setup.weather_type == 'custom':
            sounding = dataset

        self.var_typ = var_typ
        self.atm_canvas.setLabel('left', "Height", units='m')
        if self.var_typ == 't':
            #(consts.GAMMA*consts.R/consts.M_0*temperature[:])**0.5
            X = sounding[:, 1]
            self.atm_canvas.setLabel('bottom', "Speed of Sound", units='m/s')
        elif self.var_typ == 'm':
            X = sounding[:, 2]
            self.atm_canvas.setLabel('bottom', "Wind Magnitude", units='m/s')
        elif self.var_typ == 'd':
            X = sounding[:, 3]
            self.atm_canvas.setLabel('bottom', "Wind Direction", units='deg E from N')
        else:
            errorMessage('Error reading var_typ in atmPlotProfile', 2)
            return None
        Y = sounding[:, 0]

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
            for ptb_n in range(1, self.setup.perturb_times):

                if self.prefs.debug:
                    print("STATUS: Perturbation {:}".format(ptb_n))

                # generate a perturbed sounding profile
                sounding_p = perturbation_method(self.setup, dataset, self.setup.perturb_method, \
                    sounding_u=sounding_u, sounding_l=sounding_l, \
                    spread_file=self.setup.perturbation_spread_file, lat=self.setup.lat_centre, lon=self.setup.lon_centre, ensemble_file=ensemble_file, ensemble_no=ptb_n)
                sounding_p = findECMWFSound(lat, lon, sounding_p)
                
                if self.var_typ == 't':
                    X = sounding_p[:, 1]
                elif self.var_typ == 'm':
                    X = sounding_p[:, 2]
                elif self.var_typ == 'd':
                    X = sounding_p[:, 3]
                else:
                    print('error, atmPlotProfile')

                Y = sounding_p[:, 0]
                
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
    
    
    def makeValueChange(self, obj, slider):
        
        if obj == self.low_bandpass_label:
            obj.setText('Low: {:8.2f} Hz'.format(slider.value()*self.bandpass_scale))
        elif obj == self.high_bandpass_label:
            obj.setText('High: {:8.2f} Hz'.format(slider.value()*self.bandpass_scale))
        else:
            errorMessage('Bad atm slider pass in makeValueChange', 2)

        self.updatePlot()

    def makeStationObj(self, lst):

        new_lst = []

        for line in lst:
            pos = Position(line[2], line[3], line[4])
            stn = Station(line[0], line[1], pos, line[5], line[6], line[7])
            new_lst.append(stn)

        return new_lst
    
    def makePicks(self):

        if not self.checkForWorkDir():
            return None

        #     #Build seismic data path
        # self.dir_path = os.path.join(self.prefs.workdir, self.bam.setup.fireball_name)

        # # Load the station and waveform files list
        # data_file_path = os.path.join(self.dir_path, DATA_FILE)
        # if os.path.isfile(data_file_path):
            
        #     stn_list = readStationAndWaveformsListFile(data_file_path, rm_stat=self.bam.setup.rm_stat)

        # else:
        #     errorMessage('Station and waveform data file not found! Download the waveform files first!', 2)
        #     return None

        # if self.setup.stations is not None:
        #     stn_list = stn_list + self.makeStationObj(self.setup.stations)

        # Init the constants
        self.bam.setup.search_area = [self.bam.setup.lat_centre - self.bam.setup.deg_radius, 
                                      self.bam.setup.lat_centre + self.bam.setup.deg_radius,
                                      self.bam.setup.lon_centre - self.bam.setup.deg_radius,
                                      self.bam.setup.lon_centre + self.bam.setup.deg_radius]


        # sounding = parseWeather(self.setup)

        # if len(stn_list) == 0:
        #     errorMessage('No Stations to load', 2)
        #     return None 

        try:
            #turn coordinates into position objects
            self.bam.setup.traj_i = Position(self.bam.setup.lat_i, self.bam.setup.lon_i, self.bam.setup.elev_i)
            self.bam.setup.traj_f = Position(self.bam.setup.lat_f, self.bam.setup.lon_f, self.bam.setup.elev_f)
        except:
            self.bam.setup.traj_i = Position(0, 0, 0)
            self.bam.setup.traj_f = Position(0, 0, 0)
            errorMessage("Warning: Unable to build trajectory points", 1)

        self.waveformPicker()


    def waveformPicker(self, waveform_window=600):
        """

        Arguments:
            data_list: [list]

        Keyword arguments:
            waveform_window: [int] Number of seconds for the waveform window.
            difference_filter_all: [bool] If True, the Kalenda et al. (2014) difference filter will be applied
                on the data plotted in the overview plot of all waveforms.
        """

        self.v_sound = self.prefs.avg_sp_sound
        self.t0 = self.bam.setup.t0


        # Filter out all stations for which the mseed file does not exist
        filtered_stn_list = []

        names = []
        lats = []
        lons = []
  

        self.lat_centre = self.bam.setup.lat_centre
        self.lon_centre = self.bam.setup.lon_centre

        self.waveform_window = waveform_window

        self.current_station = 0
        self.current_waveform_raw = None
        self.current_waveform_delta = None
        self.current_waveform_processed = None

        # List of picks
        self.pick_list = []

        self.pick_group = 0

        # Define a list of colors for groups
        self.pick_group_colors = ['w', 'g', 'm', 'c', 'y']

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

        ### Sort stations by distance from source ###

        # Calculate distances of station from source
        self.source_dists = []


        for stn in self.bam.stn_list:

            stat_name, stat_lat, stat_lon = stn.metadata.code, stn.metadata.position.lat, stn.metadata.position.lon

            names.append(stat_name)
            lats.append(stat_lat)
            lons.append(stat_lon)

            # Calculate the distance in kilometers
            dist = greatCircleDistance(np.radians(self.bam.setup.lat_centre), np.radians(self.bam.setup.lon_centre), \
                np.radians(stat_lat), np.radians(stat_lon))

            self.source_dists.append(dist)

        # Get sorted arguments
        dist_sorted_args = np.argsort(self.source_dists)

        # Sort the stations by distance
        self.bam.stn_list = [self.bam.stn_list[i] for i in dist_sorted_args]
        self.source_dists = [self.source_dists[i] for i in dist_sorted_args]
        #############################################
        
        # Init the plot framework
        self.initPlot()

        if not hasattr(self, 'make_picks_gmap_view'):
            self.make_picks_gmap_view = QWebView()
            self.make_picks_top_graphs.addWidget(self.make_picks_gmap_view)
            self.make_picks_gmap_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        # Extract coordinates of the reference station
        gmap_filename = htmlBuilder(self.bam.setup, self.prefs, self.bam.stn_list)
        self.make_picks_gmap_view.load(QUrl().fromLocalFile(gmap_filename))

        self.make_picks_map_graph_canvas.setLabel('bottom', "Longitude", units='deg E')
        self.make_picks_map_graph_canvas.setLabel('left', "Latitude", units='deg N')

        for ii, stn in enumerate(self.bam.stn_list):


            self.station_marker[ii].setPoints(x=[stn.metadata.position.lon], y=[stn.metadata.position.lat], pen=(255, 255, 255), brush=(255, 255, 255), symbol='t')
            self.make_picks_map_graph_canvas.addItem(self.station_marker[ii], update=True)
            txt = pg.TextItem("{:}".format(stn.metadata.code))
            txt.setPos(stn.metadata.position.lon, stn.metadata.position.lat)
            self.make_picks_map_graph_canvas.addItem(txt)

        if self.prefs.frag_en:
            
            # Fragmentation plot
            for i, line in enumerate(self.bam.setup.fragmentation_point):
                self.make_picks_map_graph_canvas.scatterPlot(x=[float(line.position.lon)], y=[float(line.position.lat)],\
                    pen=(0 + i*255/len(self.bam.setup.fragmentation_point), 255 - i*255/len(self.bam.setup.fragmentation_point), 0), symbol='+')

        # Plot source location
        self.make_picks_map_graph_canvas.scatterPlot(x=[self.bam.setup.lon_centre], y=[self.bam.setup.lat_centre], symbol='+', pen=(255, 255, 0))

        # Manual trajectory search
        if self.prefs.ballistic_en:

            try:

                if self.bam.setup.trajectory.pos_i.isNone():
                    raise TypeError
                # Plot the trajectory with the bottom point known
                self.make_picks_map_graph_canvas.plot([self.bam.setup.trajectory.pos_i.lon, self.bam.setup.trajectory.pos_f.lon],\
                                                      [self.bam.setup.trajectory.pos_i.lat, self.bam.setup.trajectory.pos_f.lat],\
                                                        pen=(0, 0, 255))


                # Plot intersection with the ground
                self.make_picks_map_graph_canvas.scatterPlot(x=[self.bam.setup.trajectory.pos_f.lon], \
                                                             y=[self.bam.setup.trajectory.pos_f.lat], \
                                                                symbol='+', pen=(0, 0, 255))
            except TypeError as e:
                errorMessage('Trajectory is not defined!', 1, info='If not defining a trajectory, then turn off show ballistic waveform', detail='{:}'.format(e))
                self.prefs.ballistic_en = False


        self.bam.stn_list = calcAllTimes(self.bam, self.prefs)
        save(self)
        SolutionGUI.update(self)

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


    def initPlot(self):
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

        self.station_marker = [None]*len(self.bam.stn_list)
        self.station_waveform = [None]*len(self.bam.stn_list)
        for ii, stn in enumerate(self.bam.stn_list):
            self.make_picks_station_choice.addItem("{:}-{:}".format(stn.metadata.network, stn.metadata.code))
            self.station_marker[ii] = pg.ScatterPlotItem()
            self.station_waveform[ii] = pg.PlotCurveItem()

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
        for i in range(len(self.bam.stn_list)):
            lats.append(self.bam.stn_list[i].metadata.position.lat)
            lons.append(self.bam.stn_list[i].metadata.position.lon)
        
        self.ballistic_idx = []
        self.fragmentation_idx = []


        # Go though all stations and waveforms
        bad_stats = []

        for idx, stn in enumerate(self.bam.stn_list):

            sys.stdout.write('\rPlotting: {:} {:}              '.format(stn.metadata.network, stn.metadata.code))
            sys.stdout.flush()
            time.sleep(0.001)


            mseed = stn.stream                
  

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

            stn.offset = (start_datetime - self.bam.setup.fireball_datetime - datetime.timedelta(minutes=5)).total_seconds()

            # Skip stations with no data
            if len(waveform_data) == 0:
                continue

            waveform_data = convolutionDifferenceFilter(waveform_data)

            # Calculate the distance from the source point to this station (kilometers)
            station_dist = greatCircleDistance(np.radians(self.bam.setup.lat_centre), np.radians(self.bam.setup.lon_centre), stn.metadata.position.lat_r, stn.metadata.position.lon_r)

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

                try:
                    # Time of arrival
                    toa = station_dist/(prefs.avg_sp_sound/1000) + self.bam.setup.t0
                except:
                    toa = station_dist/(310/1000)

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

            if stn.metadata.code in self.bam.setup.high_f:
                self.fragmentation_idx.append(idx)
            if stn.metadata.code in self.bam.setup.high_b:
                self.ballistic_idx.append(idx)


        toa_line_time = np.linspace(0, max_time, 10)

        # Plot the constant sound speed line (assumption is that the release happened at t = 0)
        try:
            self.make_picks_station_graph_canvas.plot((toa_line_time)*prefs.avg_sp_sound, (toa_line_time + setup.t0), pen=(255, 0, 0))
        except:
            self.make_picks_station_graph_canvas.plot((toa_line_time)*310, (toa_line_time), pen=(255, 0, 0))
        print('')
        
        self.make_picks_station_graph_canvas.setLabel('bottom', "Distance", units='m')
        self.make_picks_station_graph_canvas.setLabel('left', "Time", units='s')

        SolutionGUI.update(self)

    def keyPressEvent(self, event):

        if event.key() == QtCore.Qt.Key_Alt:
            self.alt_pressed = True

        if event.key() == QtCore.Qt.Key_D:
            self.incrementStation()


        if event.key() == QtCore.Qt.Key_A:
            self.decrementStation()


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

        if event.key() == QtCore.Qt.Key_I:
            try:
                self.invertGraph()
            except:
                pass

        if event.key() == QtCore.Qt.Key_Up:
            self.group_no += 1
            self.group_no = self.group_no%2

            if self.group_no == 0:
                self.tog_picks.changeImg(1)
            else:
                self.tog_picks.changeImg(5)

            print("Current Pick Group: {:}".format(self.group_no))

        if event.key() == QtCore.Qt.Key_Down:
            self.group_no -= 1
            self.group_no = self.group_no%2

            if self.group_no == 0:
                self.tog_picks.changeImg(1)
            else:
                self.tog_picks.changeImg(5)

            print("Current Pick Group: {:}".format(self.group_no))

    def keyReleaseEvent(self, event):

        if event.key() == QtCore.Qt.Key_Alt:
            self.alt_pressed = False

    def incrementStation(self, event=None):
        """ Increments the current station index. """

        self.make_picks_waveform_canvas.clear()

        self.current_station += 1

        if self.current_station >= len(self.bam.stn_list):
            self.current_station = 0

        # while self.checkExists() == False:
        #     self.current_station += 1
        #     if self.current_station >= len(self.bam.stn_list):
        #         self.current_station = 0


        self.updatePlot()


    def decrementStation(self, event=None):
        """ Decrements the current station index. """

        self.make_picks_waveform_canvas.clear()

        self.current_station -= 1

        if self.current_station < 0:
            self.current_station = len(self.bam.stn_list) - 1

        # while self.checkExists() == False:
        #     self.current_station -= 1
        #     if self.current_station < 0:
        #         self.current_station = len(self.bam.stn_list) - 1

        self.updatePlot()


    def checkExists(self):
        """
        Checks if the current waveform is readable
        """

        # Extract current station
        stn = self.bam.stn_list[self.current_station]

        # Get the miniSEED file path
        mseed_file_path = os.path.join(self.dir_path, stn.file_name)

        try:
            
            if os.path.isfile(mseed_file_path):
                pass
            else:
                if self.prefs.debug:
                    print('File {:s} does not exist!'.format(mseed_file_path))
                return False

        except TypeError as e:
            if self.prefs.debug:
                print('Opening file {:s} failed with error: {:s}'.format(mseed_file_path, str(e)))
            return False

        try:
            obspy.read(mseed_file_path)

        except TypeError:
            if self.prefs.debug:
                print('mseed file could not be read:', mseed_file_path)
            return False

        return True


    def mouseClicked(self, evt):

        ############################
        # Conditional Clicks HERE
        ############################

        stn = self.bam.stn_list[self.current_station]
        channel = self.make_picks_channel_choice.currentText()
        if self.tog_picks.isChecked():
            mousePoint = self.make_picks_waveform_canvas.vb.mapToView(evt.pos())

            self.make_picks_waveform_canvas.scatterPlot(x=[mousePoint.x()], y=[0], pen=self.colors[self.group_no], brush=self.colors[self.group_no], update=True)

            pick = Pick(mousePoint.x(), stn, self.current_station, stn, self.group_no)
            self.pick_list.append(pick)
            print("New pick object made: {:} {:} {:}".format(mousePoint.x(), stn.metadata.code, self.current_station))

            if self.show_height.isChecked():
                stat_picks = []
                for pick in self.pick_list:
                    if pick.stn == stn:
                        stat_picks.append(pick)


                self.w = FragmentationStaff(self.bam.setup, [stn, self.current_station, stat_picks])
                self.w.setGeometry(QRect(100, 100, 900, 900))
                self.w.show()

            elif self.solve_height.isChecked():
                ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)
                P = self.setup.trajectory.trajInterp(div=HEIGHT_SOLVER_DIV)
                stn = self.stn_list[self.current_station]
                stn.position.pos_loc(ref_pos)
                dataset = parseWeather(self.setup)
                C = []
                max_steps = len(P)*self.setup.perturb_times + 1
                count = 0
                loadingBar("Trying Heights", 0, max_steps)
                A = self.setup.trajectory.pos_i
                B = self.setup.trajectory.pos_f

                A.pos_loc(B)
                B.pos_loc(B)

                # Get prediction of time of the meteor, so the timing of each fragmentation can be known
                length_of_meteor = np.sqrt((A.x - B.x)**2 + (A.y - B.y)**2 + (A.z - B.z)**2)
                time_of_meteor = length_of_meteor/self.setup.trajectory.v
                for ii, point in enumerate(P):
                    point.pos_loc(ref_pos)
                    for ptb_n in range(self.setup.perturb_times):
                        
                        self.sounding = self.perturbGenerate(ptb_n, dataset, self.perturbSetup())
                        zProfile, _ = getWeather(np.array([point.lat, point.lon, point.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), self.setup.weather_type, \
                                [ref_pos.lat, ref_pos.lon, ref_pos.elev], self.sounding, convert=False)
                        
                        #zProfile = zInterp(stn.position.z, point.z, zProfile, div=37)

                        f_time, _, _, _ = cyscan(np.array([point.x, point.y, point.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True, \
                            n_theta=self.setup.n_theta, n_phi=self.setup.n_phi, h_tol=self.setup.h_tol, v_tol=self.setup.v_tol)

                        correction = time_of_meteor - A.z/self.setup.trajectory.pos_i.elev*(time_of_meteor)
                        C.append(f_time + correction)
                        count += 1
                        loadingBar("Trying Heights", count, max_steps)
                C = np.array(C)

                idx = np.nanargmin(np.abs(C - pick.time))
                print("Error in time: {:.2f} s".format(np.abs(C - pick.time)[idx]))
                height_idx = idx//self.setup.perturb_times
                pert_idx = idx%self.setup.perturb_times

                self.position.append(P[height_idx])

                self.x = AllWaveformViewer(self.setup, self.stn_list, self.position, pert_idx)
                self.x.setGeometry(QRect(100, 100, 900, 900))
                self.x.show()

        elif self.tog_rm_picks.isChecked():

            self.make_picks_waveform_canvas.clear()
            for ii, pick in enumerate(self.pick_list):
                if pick.stn_no == self.current_station:
                    self.pick_list.pop(ii)
                    if self.prefs.debug:
                        print('Pick removed!')

                self.make_picks_waveform_canvas.scatterPlot(x=[pick.time], y=[0], pen=self.colors[self.group_no], brush=self.colors[self.group_no], update=True)
            self.drawWaveform(station_no=self.current_station)

        elif self.gnd_mot_picks.isChecked():

            # Open ground motion Dialog
            current_chn_start = channel[0:2]
            channel_opts = [self.make_picks_channel_choice.itemText(i) for i in range(self.make_picks_channel_choice.count())]
            
            # Check if zne/z12 is available
            count = 0

            for chn in channel_opts:
                if current_chn_start in chn:
                    count += 1

            if count == 3:
                
                self.gr = ParticleMotion(self.make_picks_map_graph_canvas, self.bam, stn, channel, t_arrival=self.source_dists[self.current_station]/(310/1000), group_no=self.group_no)
                self.gr.setGeometry(QRect(100, 100, 1200, 700))
                self.gr.show()
            elif count < 3:
                errorMessage("Not enough channel data for particle motion!", 2, \
                        detail="Three orthogonal stations are needed to do particle motion!")
            else:
                errorMessage("If you are seeing this, then somehow more than 3 channels have been selected",\
                         2, detail="")

        elif self.bandpass_picks.isChecked():


            # Open bandpass GUI
            self.bp = BandpassWindow(self.bam, stn, channel, t_arrival=self.source_dists[self.current_station]/(310/1000))
            self.bp.setGeometry(QRect(100, 100, 1200, 700))
            self.bp.show()

        elif self.polmap_picks.isChecked():

            ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)
            points = []

            # Calculate all points here
            for stn in self.bam.stn_list:

                if not hasattr(stn, "polarization"):
                    stn.polarization = Polarization()

                if len(stn.polarization.azimuth) > 0: 

                    D = propegateBackwards(ref_pos, stn, self.bam)

                    for line in D:
                        if not np.isnan(line[0]): 
                            P = Position(0, 0, 0)
                            P.x = line[0]
                            P.y = line[1]
                            P.z = line[2]
                            P.pos_geo(ref_pos)
                            S = Supracenter(P, line[3])
                            points.append([S, stn.color])

            # Pass grid to polmap
            self.pm = Polmap(self.bam, points)
            self.pm.setGeometry(QRect(100, 100, 1200, 700))
            self.pm.show()



        elif self.alt_pressed:

            # Create annotation
            mousePoint = self.make_picks_waveform_canvas.vb.mapToView(evt.pos())

            # pick = Pick(mousePoint.x(), self.stn_list[self.current_station], self.current_station, self.stn_list[self.current_station], self.group_no)


            self.a = AnnoteWindow(mousePoint.x(), self.stn_list[self.current_station])
            self.a.setGeometry(QRect(400, 500, 400, 500))
            self.a.show()
            
            self.drawWaveform()
            self.alt_pressed = False
            

    def addAnnotes(self):

        pass
        # for an in self.bam.stn_list[self.current_station].annotations:
        #     line = pg.InfiniteLine(pos=(an.time, 0), angle=90, pen=an.color, label=an.title)
        #     self.make_picks_waveform_canvas.addItem(line, update=True)
        #     # line.make_picks_waveform_canvas.mpl_connect('button_press_event', self.onClick)
        

    def onClick(self):
        print('It Worked!')

    def drawWaveform(self, channel_changed=0, waveform_data=None, station_no=0):
        """ Draws the current waveform from the current station in the waveform window. Custom waveform 
            can be given an drawn, which is used when bandpass filtering is performed. 

        """

        station_no = self.current_station

        # Clear waveform axis
        self.make_picks_waveform_canvas.clear()

        # Extract current station
        stn = self.bam.stn_list[station_no]

        # Get the miniSEED file path
        
        # Try reading the mseed file, if it doesn't work, skip to the next frame
        try:
            mseed = stn.stream

        except TypeError:
            if self.prefs.debug:
                print('mseed file could not be read:', mseed_file_path)
            return None

        if channel_changed == 0:
            # Populate channel list
            self.make_picks_channel_choice.blockSignals(True)
            self.make_picks_channel_choice.clear()
            for i in range(len(mseed)):
                self.make_picks_channel_choice.addItem(mseed[i].stats.channel)
            self.make_picks_channel_choice.blockSignals(False)
        
        current_channel = self.make_picks_channel_choice.currentIndex()
        chn_selected = self.make_picks_channel_choice.currentText()

        resp = stn.response

        st = mseed

        # A second stream containing channels with the response
        st2 = mseed.select(inventory=resp.select(channel=chn_selected))

        # Use st2 if able to, else use st
        st2 = st2.select(channel=chn_selected)
        st = st.select(channel=chn_selected)
        # Unpact miniSEED data

        if len(st2) > 0 and resp is not None and self.rm_resp.isChecked():
            st = st2


            #TODO - bug, this shouldn't run every time the waveform is shown
            #Obspy says that this is because the response is removed on the actual data, use .copy() 
            

            st = st[0].remove_response(inventory=resp, output="VEL")
            st.remove_sensitivity(resp) 
            rm_resp = True
        else:
            st = st[0]
            rm_resp = False

        st.detrend()

        delta = st.stats.delta
        start_datetime = st.stats.starttime.datetime
        end_datetime = st.stats.endtime.datetime

        stn.offset = (start_datetime - self.bam.setup.fireball_datetime).total_seconds()

        # Check if the waveform data is already given or not
        if waveform_data is None or channel_changed != 2:
            waveform_data = mseed[current_channel].data

            # Store raw data for bookkeeping on first open
            self.current_waveform_raw = waveform_data

        self.current_waveform_delta = delta
        self.current_waveform_time = np.arange(0, mseed[current_channel].stats.npts / mseed[current_channel].stats.sampling_rate, \
             delta)
        # self.current_waveform_time = np.arange(0, (end_datetime - start_datetime).total_seconds(), \
        #     delta)

        # Construct time array, 0 is at start_datetime
        time_data = np.copy(self.current_waveform_time)

        # Cut the waveform data length to match the time data
        waveform_data = waveform_data[:len(time_data)]
        time_data = time_data[:len(waveform_data)] + stn.offset

        # Store currently plotted waveform
        self.current_waveform_processed = waveform_data

        # Calculate the time of arrival assuming constant propagation with the given speed of sound
        try:
            t_arrival = self.source_dists[self.current_station]/(self.v_sound/1000) + self.t0
        except:
            t_arrival = self.source_dists[self.current_station]/(310/1000)
        

        # Plot the waveform
        # self.current_station_waveform = pg.PlotDataItem(x=time_data, y=waveform_data, pen='w')
        # 2.370980392E10

        self.current_station_waveform = pg.PlotDataItem(x=time_data, y=waveform_data, pen='w')
        self.make_picks_waveform_canvas.addItem(self.current_station_waveform)
        # self.make_picks_waveform_canvas.plot(x=[t_arrival, t_arrival], y=[np.min(waveform_data), np.max(waveform_data)], pen=pg.mkPen(color=(255, 0, 0), width=2))
        self.make_picks_waveform_canvas.setXRange(t_arrival-100, t_arrival+100, padding=1)

        self.make_picks_waveform_canvas.setLabel('bottom', "Time after {:} s".format(self.bam.setup.fireball_datetime))

        if rm_resp:
            self.make_picks_waveform_canvas.setLabel('left', "Ground Motion", units='m/s')
        else:
            self.make_picks_waveform_canvas.setLabel('left', "Signal Response")

        # self.make_picks_waveform_canvas.setLabel('left', pg.LabelItem("Overpressure", size='20pt'), units='Pa')
        self.make_picks_waveform_canvas.plot(x=[-10000, 10000], y=[0, 0], pen=pg.mkPen(color=(100, 100, 100)))


        font=QtGui.QFont()
        font.setPixelSize(20)
        self.make_picks_waveform_canvas.getAxis("bottom").tickFont = font
        self.make_picks_waveform_canvas.getAxis("left").tickFont = font
        self.make_picks_waveform_canvas.getAxis('bottom').setPen(self.color.WHITE) 
        self.make_picks_waveform_canvas.getAxis('left').setPen(self.color.WHITE)

        for pick in self.pick_list:
            if pick.stn_no == self.current_station:
                self.make_picks_waveform_canvas.scatterPlot(x=[pick.time], y=[0], pen='r', update=True)

        SolutionGUI.update(self)
        # Initialize variables
        b_time = 0

        # Extract coordinates of the reference station
        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

        # Calculate ground distances
        try:
            stn.stn_ground_distance(self.bam.setup.trajectory.pos_f)

            print('####################')
            print("Current Station: {:}-{:}".format(stn.network, stn.code))
            print("Ground Distance: {:7.3f} km".format(stn.ground_distance/1000))
        except:

            pass

        # If manual ballistic search is on
        if self.prefs.ballistic_en and self.show_ball.isChecked():

            b_time = stn.times.ballistic[0][0][0]

            if b_time == b_time:
                self.make_picks_waveform_canvas.plot(x=[b_time, b_time], y=[np.min(waveform_data), np.max(waveform_data)], pen=pg.mkPen(color=(0, 0, 255), width=2) , label='Ballistic')
                print("Ballistic Arrival: {:.3f} s".format(b_time))
            else:
                print("No Ballistic Arrival")

            data, remove = chauvenet(stn.times.ballistic[0][1][0])
            try:
                print('Perturbation Arrival Range: {:.3f} - {:.3f}s'.format(np.nanmin(data), np.nanmax(data)))
                print('Removed points {:}'.format(remove))
            except ValueError:
                print('No Perturbation Arrivals')

            for i in range(len(data)):
                if self.show_perts.isChecked():
                    try:
                        self.make_picks_waveform_canvas.plot(x=[data[i]]*2, \
                         y=[np.min(waveform_data), np.max(waveform_data)], pen=pg.mkPen(color=(0, 0, 255), style=QtCore.Qt.DotLine) )
                    except:
                        pass
            # Fragmentation Prediction

            # If manual fragmentation search is on
        if self.prefs.frag_en and self.show_frags.isChecked():

            for i, frag in enumerate(self.bam.setup.fragmentation_point):

                f_time = stn.times.fragmentation[i][0][0]

                v_time = (frag.position.elev - stn.metadata.position.elev)/310
                h_time = frag.position.ground_distance(stn.metadata.position)/1000
                p_time = h_time + v_time
                print('++++++++++++++++')
                print('Fragmentation {:} ({:6.2f} km)'.format(i+1, frag.position.elev/1000))
                frag.position.pos_loc(stn.metadata.position)
                stn.metadata.position.pos_loc(stn.metadata.position)
                xyz_range = np.sqrt((frag.position.x - stn.metadata.position.x)**2 + \
                                    (frag.position.y - stn.metadata.position.y)**2 + \
                                    (frag.position.z - stn.metadata.position.z)**2)
                print('Range {:7.3f} km'.format(xyz_range/1000))
                if not np.isnan(f_time):
                    # Plot Fragmentation Prediction
                    self.make_picks_waveform_canvas.plot(x=[f_time]*2, y=[np.min(waveform_data), np.max(waveform_data)], pen=pg.mkPen(color=self.pick_group_colors[(i+1)%4], width=2), label='Fragmentation')
                    
                    if self.show_prec.isChecked():
                        # Plot Precursor Arrivals
                        self.make_picks_waveform_canvas.plot(x=[p_time]*2, y=[np.min(waveform_data), np.max(waveform_data)], pen=(210, 235, 52), label='Fragmentation')


                    stn.stn_distance(frag.position)
                    #print("Range: {:7.3f} km".format(stn.distance/1000))                   
                    print('Arrival: {:.3f} s'.format(f_time))

                else:
                    pass
                    print('No Fragmentation {:} ({:6.2f} km) Arrival'.format(i+1, frag.position.elev/1000))

                data, remove = self.obtainPerts(stn.times.fragmentation, i)
                try:
                    print('Perturbation Arrival Range: {:.3f} - {:.3f}s'.format(np.nanmin(data), np.nanmax(data)))
                    print('Removed points {:}'.format(remove))
                except ValueError:
                    print('No Perturbation Arrivals')
                for j in range(len(data)):
                    if self.show_perts.isChecked():
                        
                        try:
                            if not np.isnan(data[j]):
                                self.make_picks_waveform_canvas.plot(x=[data[j]]*2, y=[np.min(waveform_data),\
                                    np.max(waveform_data)], alpha=0.3,\
                                    pen=pg.mkPen(color=self.pick_group_colors[(i+1)%4], style=QtCore.Qt.DotLine), zorder=3)
                        except IndexError:
                            errorMessage("Error in Arrival Times Index", 2, detail="Check that the arrival times file being used aligns with stations and perturbation times being used. A common problem here is that more perturbation times were selected than are available in the given Arrival Times Fireball. Try setting perturbation_times = 0 as a first test. If that doesn't work, try not using the Arrival Times file selected in the toolbar.")
                            return None

        self.addAnnotes()
    def obtainPerts(self, data, frag):
        data_new = []

        for i in range(len(data[frag][1])):
            data_new.append(data[frag][1][i][0])
        data, remove = chauvenet(data_new)

        return data, remove

    def showSpectrogram(self, event=None):
        """ Show the spectrogram of the waveform in the current window. """


        # ### Show the spectrogram ###
        wave_arr = self.current_waveform_raw
        fig = plt.figure()
        ax_spec = fig.add_subplot(111)

        ax_spec.specgram(wave_arr, Fs=20, cmap=plt.cm.inferno)
   
        ax_spec.set_xlabel('Time (s)')
        ax_spec.set_ylabel('Frequency (Hz)')

        fig.show()

    def filterBandpass(self, event=None):
        """ Run bandpass filtering using values set on sliders. """

        # Get bandpass filter values
        bandpass_low = float(self.low_bandpass_edits.text())
        bandpass_high = float(self.high_bandpass_edits.text())


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
        self.drawWaveform(channel_changed=2, waveform_data=waveform_data, station_no=self.current_station)


    def filterConvolution(self, event=None):
        """ Apply the convolution filter on data as suggested in Kalenda et al. (2014). """

        waveform_data = convolutionDifferenceFilter(self.current_waveform_raw)

        self.drawWaveform(channel_changed=2, waveform_data=waveform_data, station_no=self.current_station)


    def updatePlot(self, draw_waveform=True):
        """ Update the plot after changes. """

        if errorCodes(self, 'current_station', debug=self.prefs.debug):
            return None

        self.make_picks_waveform_canvas.clear()

        # Mark the position of the current station on the map
        self.make_picks_station_choice.setCurrentIndex(self.current_station)

        for stn_mk in (i for i in self.station_marker if i is not None):
            if stn_mk == self.station_marker[self.current_station]:
                stn_mk.setPen((255, 0, 0))
                stn_mk.setBrush((255, 0, 0))
                stn_mk.setZValue(1)
            elif stn_mk in [self.station_marker[i] for i in self.ballistic_idx]:
                stn_mk.setPen((0, 0, 255))
                stn_mk.setBrush((0, 0, 255))
                stn_mk.setZValue(0)
            elif stn_mk in [self.station_marker[i] for i in self.fragmentation_idx]:
                stn_mk.setPen((0, 255, 0))
                stn_mk.setBrush((0, 255, 0))
                stn_mk.setZValue(0)
            else:
                stn_mk.setPen((255, 255, 255))
                stn_mk.setBrush((255, 255, 255))
                stn_mk.setZValue(0)


        for stn_mk in (i for i in self.station_waveform if i is not None):
            if stn_mk != self.station_waveform[self.current_station]:
                stn_mk.setPen((255, 255, 255))
                stn_mk.setZValue(0)
            else:
                stn_mk.setPen((255, 0, 0))
                stn_mk.setZValue(1)

        # Plot the waveform from the current station
        if draw_waveform:
            self.drawWaveform(station_no=self.current_station)

        # self.showTitle()

        SolutionGUI.update(self)

    def exportCSV(self, event):
        """ Save picks to a CSV file. """

        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = checkExt(dlg[0], '.csv')

        # Open the output CSV
        try:
            with open(os.path.join(file_name), 'w') as f:

                # Write the header
                f.write('Pick group, Network, Code, Lat, Lon, Elev, Pick JD, Pick time, station_number \n')

                # Go through all picks
                for pick in self.pick_list:

                    # Calculate Julian date of the pick time
                    pick_jd = datetime2JD(self.bam.setup.fireball_datetime + datetime.timedelta(seconds=pick.time))

                    stn = pick.stn

                    # Write the CSV entry
                    f.write("{:d}, {:s}, {:s}, {:.6f}, {:.6f}, {:.2f}, {:.8f}, {:}, {:}\n".format(0, stn.metadata.network, \
                        stn.metadata.code, stn.metadata.position.lat, stn.metadata.position.lon, stn.metadata.position.elev, pick_jd, pick.time, pick.stn_no))
        except FileNotFoundError as e:
            errorMessage('Could not find file!', 2, detail='{:}'.format(e))
            return None

        errorMessage('Output to CSV!', 0, title='Exported!')

    def exportToAllTimes(self):
        
        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = dlg[0]

        np.save(file_name, self.arrTimes)

        errorMessage("Saved Output File", 0)

    def exportImage(self):

        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = dlg[0]

        exporter = pg.exporters.SVGExporter(self.make_picks_waveform_view.scene())

        file_name = file_name + '.svg'
        exporter.export(file_name)


    def select(self, inverted):

        if errorCodes(self, 'current_station_waveform'):
            return None

        if inverted:
            self.make_picks_waveform_view.setBackground(self.color.BLACK)
            self.current_station_waveform.setPen(self.color.WHITE)
            self.make_picks_waveform_canvas.getAxis('bottom').setPen(self.color.WHITE) 
            self.make_picks_waveform_canvas.getAxis('left').setPen(self.color.WHITE)
        else:
            self.make_picks_waveform_view.setBackground(self.color.WHITE)
            self.current_station_waveform.setPen(self.color.BLACK)
            self.make_picks_waveform_canvas.getAxis('bottom').setPen(self.color.BLACK) 
            self.make_picks_waveform_canvas.getAxis('left').setPen(self.color.BLACK)

    def invertGraph(self):

        self.select(self.inverted)

        self.inverted = not self.inverted
        
    def showTitle(self):

        if errorCodes(self, 'current_station_waveform'):
            return None

        stn = self.stn_list[self.current_station]
        if self.showtitled:
            self.make_picks_waveform_canvas.setTitle('')
        else:
            if self.inverted:
                self.make_picks_waveform_canvas.setTitle('{:}-{:} [{:}] {:}'.format(stn.network, stn.code, stn.channel, stn.position), color=(0, 0, 0))
            else:
                self.make_picks_waveform_canvas.setTitle('{:}-{:} [{:}] {:}'.format(stn.network, stn.code, stn.channel, stn.position), color=(255, 255, 255))

        self.showtitled = not self.showtitled

    def perturbSetup(self):
        """ Pulls the correct file names for the perturbing function to read, depending on the perturbation type
        """

        if self.bam.setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(self.bam.setup, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(self.bam.setup, time=-1)

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

    def perturbGenerate(self, ptb_n, dataset, perturb_data, line=False):
        """ Generates a perturbed cubic atmospheric profile (similar to 'dataset') based off of the perturb data 
        """

        sounding_l, sounding_u, ensemble_file = perturb_data[0], perturb_data[1], perturb_data[2]

        # Perturbed soundings
        if ptb_n > 0:
            

            sounding_p = perturbation_method(self.setup, dataset, self.setup.perturb_method, \
                sounding_u=sounding_u, sounding_l=sounding_l, \
                spread_file=self.setup.perturbation_spread_file, lat=self.setup.lat_centre, lon=self.setup.lon_centre, \
                ensemble_file=ensemble_file, ensemble_no=ptb_n, line=line)

        # Nominal sounding
        else:
            sounding_p = dataset


        return sounding_p




    def supraSearch(self):

        # Error Parsing
        ###############

        if self.setup.lat_centre is None or self.setup.lon_centre is None:
            errorMessage("Lat center or Lon center are not defined!", 1, info="Please define both to use this function")
            return None

        # if self.setup.search_min

        # Parse the picked station data
        s_info, s_name, weights, ref_pos = self.getStationData()

        if s_info is None:
            return None

        # Check if manual search is defined
        if self.setup.manual_fragmentation_search is None:
            errorMessage("No manual fragmentation point defined!", 2, detail="Please specify a Manual Fragmentation Search in the Sources section of Variables")        
            return None

        # Collect Results
        #################

        ref_pos = Position(ref_pos[0], ref_pos[1], ref_pos[2])
     
        # Generate weather profile (cubic)
        dataset = parseWeather(self.setup)

        # Initialize results array
        results = [None]*(self.setup.perturb_times + 1)

        # Run through all perturbations
        for ptb_n in range(self.setup.perturb_times + 1):
            
            # for ptb_n == 0, nominal sounding, ptb_n > 1, perturbed sounding
            sounding_p = self.perturbGenerate(ptb_n, dataset, self.perturbSetup())

            # Return results for a specific atmosphere
            results[ptb_n] = psoSearch(s_info, weights, s_name, self.setup, sounding_p, ref_pos, manual=True)

            print("Error Function: {:5.2f} (Perturbation {:})".format(results[ptb_n].f_opt, ptb_n))
            print("Opt: Latitude: {:.4f} Longitude: {:.4f} Elevation: {:.2f} Mean Error: {:.4f}"\
                .format(results[ptb_n].x_opt.lat, results[ptb_n].x_opt.lon, results[ptb_n].x_opt.elev, results[ptb_n].motc))
        
        n_stations = len(s_info)
        xstn = s_info[:n_stations, :3]
            
        # Display Results
        #################

        self.scatterPlot(self.setup, results, n_stations, xstn, s_name, dataset)

        self.residPlot(results, s_name, xstn, self.prefs.workdir, n_stations)

        defTable(self.tableWidget, n_stations + 1, 5, headers=['Station Name', "Latitude", "Longitude", "Elevation", "Residuals"])
        
        setTableRow(self.tableWidget, 0, terms=["Total", results[0].x_opt.lat, results[0].x_opt.lon, results[0].x_opt.elev, results[0].f_opt])

        for i in range(n_stations):
            setTableRow(self.tableWidget, i + 1, terms=[s_name[i], xstn[i][0], xstn[i][1], xstn[i][2], results[0].r[i]])



if __name__ == '__main__':

    app = QApplication(sys.argv)

    splash_pix = QPixmap(os.path.join('supra', 'Fireballs','docs', '_images', 'wmpl.png'))
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()

    app.processEvents()

    gui = SolutionGUI()

    gui.showFullScreen()
    gui.showMaximized()
    gui.show()

    splash.finish(gui)

    sys.exit(app.exec_())