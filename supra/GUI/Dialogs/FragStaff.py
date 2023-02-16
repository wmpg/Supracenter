import csv
import os

import numpy as np
import matplotlib.pyplot as plt

import pyqtgraph as pg
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *


from supra.Utils.AngleConv import chauvenet
from supra.Utils.Classes import Position

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT
from supra.GUI.Tools.GUITools import *

from supra.Lightcurve.light_curve import processLightCurve, readLightCurve

from supra.Stations.ProcessStation import procTrace, procStream, findChn, findDominantPeriodPSD

from supra.Geminus.overpressure2 import overpressureihmod_Ro
from supra.Geminus.geminusSearch import periodSearch, presSearch

class FragmentationStaff(QWidget):

    """ A visualization of arrival times vs. heights. Makes it easier to see what height along
    the trajectory (fragmentation or ballistic) the arrival is from.
    """


    def __init__(self, setup, pack, bam):

        QWidget.__init__(self)

        self.buildGUI()

        # Take important values from main window class
        stn, self.current_station, self.pick_list, self.channel = pack

        # The first pick created is the one all times are based around (reference time)
        nom_pick = self.pick_list[0]

        main_layout = QGridLayout()
        
        # Pass setup value
        self.setup = setup
        self.bam = bam

        ###########
        # Build GUI
        ###########

        self.height_plot = MatplotlibPyQT()

        main_layout.addWidget(self.height_plot, 1, 1, 1, 100)

        export_button = QPushButton('Export')
        main_layout.addWidget(export_button, 3, 1, 1, 25)
        export_button.clicked.connect(self.export)

        X = []
        Y = []
        Y_M = []

        # Create a local coordinate system with the bottom of the trajectory as the reference
        A = self.setup.trajectory.pos_i
        B = self.setup.trajectory.pos_f

        A.pos_loc(B)
        B.pos_loc(B)

        # All points are placed here
        self.dots_x = []
        self.dots_y = []

        #########################
        # Light Curve Plot
        #########################
        if len(self.setup.light_curve_file) > 0 and hasattr(self.setup, "light_curve_file"):
            self.height_plot.ax1 = self.height_plot.figure.add_subplot(211)
            self.height_plot.ax2 = self.height_plot.figure.add_subplot(212, sharex=self.height_plot.ax1)
            # lc_plot = MatplotlibPyQT()
            # main_layout.addWidget(lc_plot, 1, 1, 1, 100)
            # self.light_curve_view = pg.GraphicsLayoutWidget()
            # self.light_curve_canvas = self.light_curve_view.addPlot()



            light_curve = readLightCurve(self.setup.light_curve_file)

            if light_curve is None:
                print("Light curve is None")


            light_curve_list = processLightCurve(light_curve)

            if light_curve_list is None:
                print("Light curve list is None")


            for L in light_curve_list:

                self.height_plot.ax1.plot(L.h, L.I, label=L.station)


                # light_curve_curve = pg.ScatterPlotItem(x=L.M, y=L.t)
                # self.light_curve_canvas.addItem(light_curve_curve)

            self.height_plot.ax1.legend()
            # plt.gca().invert_yaxis()


            # main_layout.addWidget(self.light_curve_view, 1, 101, 1, 10)

            # blank_spacer = QWidget()
            # main_layout.addWidget(blank_spacer, 2, 101, 2, 10)


            self.height_plot.ax1.set_xlim((-10, 100))
            self.height_plot.ax1.set_xlabel("Height [km]")
            self.height_plot.ax1.set_ylabel("Intensity")
        else:
            self.height_plot.ax2 = self.height_plot.figure.add_subplot(111)

        #########################
        # Station Plot
        #########################


        try:

            # Bandpass the waveform from 2 - 8 Hz (not optimal, but does a good job 
            # in showing arrivals clearly for most cases)
            st, resp, gap_times = procStream(stn, ref_time=self.setup.fireball_datetime)
            st = findChn(st, self.channel)
            waveform_data, time_data = procTrace(st, ref_datetime=self.setup.fireball_datetime,\
                    resp=resp, bandpass=[2, 8])

            # Scale the data so that the maximum is brought out to SCALE_LEN

            SCALE_LEN = 10 # km
            max_val = 0
            for ii in range(len(waveform_data)):
                wave = waveform_data[ii]
                for point in wave:
                    if abs(point) > max_val:
                        max_val = abs(point)

            scaling = SCALE_LEN/max_val

            # Plot all waveform data segments (for gaps in data)
            for ii in range(len(waveform_data)):
                self.height_plot.ax2.plot(waveform_data[ii]*scaling, time_data[ii] - nom_pick.time)

        except ValueError:
            print("Could not filter waveform!")



        #########################
        # Light Curve Plot
        #########################

        if len(self.setup.light_curve_file) > 0 or not hasattr(self.setup, "light_curve_file"):


            light_curve = readLightCurve(self.setup.light_curve_file)

            light_curve_list = processLightCurve(light_curve)

            if light_curve_list is not None:


                for L in light_curve_list:
                    self.height_plot.ax1.scatter(L.h, L.I, label=L.station)
                    # light_curve_curve = pg.ScatterPlotItem(x=L.M, y=L.t)
                    # self.light_curve_canvas.addItem(light_curve_curve)

                
                self.height_plot.ax1.grid(alpha=0.2)
                self.height_plot.ax1.legend()

            # plt.gca().invert_yaxis()


            # main_layout.addWidget(self.light_curve_view, 1, 101, 1, 10)

            # blank_spacer = QWidget()
            # main_layout.addWidget(blank_spacer, 2, 101, 2, 10)


        
        ########################
        # Generate Hyperbola
        ########################

        # D_0 = A

        # stn.metadata.position.pos_loc(B)

        # theta = self.setup.trajectory.zenith.rad
        # h_0 = A.z

        # h = np.arange(0, 100000)
        # v = self.setup.trajectory.v
        # k = stn.metadata.position - D_0
        # n = Position(0, 0, 0)
        # n.x, n.y, n.z = self.setup.trajectory.vector.x, self.setup.trajectory.vector.y, self.setup.trajectory.vector.z
        # n.pos_geo(B)
        # c = 350

        # T = (h - h_0)/(-v*np.cos(theta)) + (k - n*((h - h_0)/(-np.cos(theta)))).mag()/c - nom_pick.time
        
        # # estimate_plot = pg.PlotDataItem(x=h, y=T)
        # # self.height_canvas.addItem(estimate_plot, update=True)

        # self.height_plot.ax2.scatter(h/1000, T)


        #######################
        # Plot nominal points
        #######################
        NO_WIND_SPEED = 330
        # base_points = pg.ScatterPlotItem()
        angle_off = []
        no_wind_points = []
        u = np.array([self.setup.trajectory.vector.x,
                      self.setup.trajectory.vector.y,
                      self.setup.trajectory.vector.z])
        for i in range(len(self.setup.fragmentation_point)):

            f_time = stn.times.fragmentation[i][0][0]

            X = self.setup.fragmentation_point[i].position.elev
            Y = f_time - nom_pick.time
            
            self.dots_x.append(X)
            self.dots_y.append(Y)


            travel_dis = self.setup.fragmentation_point[i].position.pos_distance(stn.metadata.position)

            travel_time = travel_dis/NO_WIND_SPEED

            no_wind_points.append([self.setup.fragmentation_point[i].position.elev, \
                         travel_time - nom_pick.time + self.setup.fragmentation_point[i].time])

            az = stn.times.fragmentation[i][0][1]
            tf = stn.times.fragmentation[i][0][2]

            az = np.radians(az)
            tf = np.radians(180 - tf)
            v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

            angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))


            print("Points", X, Y, np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))

            # base_points.addPoints(x=[X], y=[Y], pen=(255, 0, 238), brush=(255, 0, 238), symbol='o')

        
        ptb_colors = [(0, 255, 26, 150), (3, 252, 176, 150), (252, 3, 3, 150), (176, 252, 3, 150), (255, 133, 3, 150),
                      (149, 0, 255, 150), (76, 128, 4, 150), (82, 27, 27, 150), (101, 128, 125, 150), (5, 176, 249, 150)]
        

        # base_points.setZValue(1)
        # self.height_canvas.addItem(base_points, update=True)

        #########################
        # Plot Precursor Points
        #########################

        # pre_points = pg.ScatterPlotItem()
        
        for i in range(len(self.setup.fragmentation_point)):
        
            PRE_SPEED = 3100

            X = self.setup.fragmentation_point[i].position.elev
            
            # Distance between frag point and the ground below it
            v_dist = X - stn.metadata.position.elev
            
            # Horizontal distance betweent the new ground point and the stn
            h_dist = self.setup.fragmentation_point[i].position.ground_distance(stn.metadata.position)
            
            # Speed of wave in air
            v_time = v_dist/NO_WIND_SPEED

            # Speed of wave in ground
            h_time = h_dist/PRE_SPEED
            # Total travel time
            Y = v_time + h_time - nom_pick.time
            
            if i == 0:
                self.height_plot.ax2.scatter(np.array(X)/1000, np.array(Y), c='y', label="Precursor Points (Speed = {:.2f} km/s)".format(PRE_SPEED/1000))
            else:
                self.height_plot.ax2.scatter(np.array(X)/1000, np.array(Y), c='y')
            # pre_points.addPoints(x=[X], y=[Y], pen=(210, 235, 52), brush=(210, 235, 52), symbol='o')

        # self.height_canvas.addItem(pre_points, update=True)

        #########################
        # Perturbation points
        #########################
        # prt_points = pg.ScatterPlotItem()
        for i in range(len(self.setup.fragmentation_point)):
            data, remove = self.obtainPerts(stn.times.fragmentation, i)
            azdata, remove = self.obtainPerts(stn.times.fragmentation, i, pt=1)
            tfdata, remove = self.obtainPerts(stn.times.fragmentation, i, pt=2)
            Y = []
            X = self.setup.fragmentation_point[i].position.elev
            for pt, az, tf in zip(data, azdata, tfdata):
                
                Y = (pt - nom_pick.time)

                self.dots_x.append(X)
                self.dots_y.append(Y)

                az = np.radians(az)
                tf = np.radians(180 - tf)
                v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

                angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))

        colour_angle = abs(90 - np.array(angle_off))
        sc = self.height_plot.ax2.scatter(np.array(self.dots_x)/1000, np.array(self.dots_y), c=colour_angle, cmap='viridis_r', label="Nominal Atmosphere Arrivals")
        cbar = self.height_plot.figure.colorbar(sc, orientation="horizontal", pad=0.2)
        cbar.ax.set_xlabel("Difference from 90 deg [deg]")
        no_wind_points = np.array(no_wind_points)
        self.height_plot.ax2.scatter(no_wind_points[:, 0]/1000, no_wind_points[:, 1], c='w', label="Isotropic Atmosphere Arrivals (Speed = {:.2f} km/s)".format(NO_WIND_SPEED/1000))

        for pick in self.pick_list:
            if pick.group == 0:
                self.height_plot.ax2.axhline(pick.time - nom_pick.time, c='g', label="User Picks")
                # self.height_canvas.addItem(pg.InfiniteLine(pos=(0, pick.time - nom_pick.time), angle=0, pen=QColor(0, 255, 0)))
            else:
                self.height_plot.ax2.axhline(pick.time - nom_pick.time, c='b')
                # self.height_canvas.addItem(pg.InfiniteLine(pos=(0, pick.time - nom_pick.time), angle=0, pen=QColor(0, 0, 255)))


        self.dots = np.array([self.dots_x, self.dots_y])

        #####################
        # Angle Calculation
        #####################

        u = np.array([self.setup.trajectory.vector.x,
                      self.setup.trajectory.vector.y,
                      self.setup.trajectory.vector.z])

        angle_off = []
        X = []
        for i in range(len(self.setup.fragmentation_point)):
            az = stn.times.fragmentation[i][0][1]
            tf = stn.times.fragmentation[i][0][2]

            az = np.radians(az)
            tf = np.radians(180 - tf)
            v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

            angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))
            X.append(self.setup.fragmentation_point[i].position.elev)
        angle_off = np.array(angle_off)

        ###############################
        # Find optimal ballistic angle
        ###############################
        try:
            best_indx = np.nanargmin(abs(angle_off - 90))
            print("Optimal Ballistic Height {:.2f} km with angle of {:.2f} deg".format(X[best_indx]/1000, angle_off[best_indx]))
            
            self.height_plot.ax2.axvline(X[best_indx]/1000, c='b', label="Optimal Ballistic Solution")
            # self.angle_canvas.addItem(pg.InfiniteLine(pos=(X[best_indx], 0), angle=90, pen=QColor(0, 0, 255)))
            # self.height_canvas.addItem(pg.InfiniteLine(pos=(X[best_indx], 0), angle=90, pen=QColor(0, 0, 255)))
            # self.angle_canvas.scatterPlot(x=X, y=angle_off, pen=(255, 255, 255), symbol='o', brush=(255, 255, 255))

            best_arr = []
            angle_arr = []
            
        except ValueError:
            best_indx = None

        angle_off = 0
        height = None
        for i in range(len(self.setup.fragmentation_point)):
            for j in range(len(stn.times.fragmentation[i][1])):
                az = stn.times.fragmentation[i][1][j][1]
                tf = stn.times.fragmentation[i][1][j][2]
                az = np.radians(az)
                tf = np.radians(180 - tf)
                v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

                angle_off_new = np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v)))))

                # self.angle_canvas.scatterPlot(x=[self.setup.fragmentation_point[i].position.elev], y=[angle_off_new], symbol='o')

                if abs(angle_off_new - 90) < abs(angle_off - 90) and not np.isnan(angle_off_new):
                    angle_off = angle_off_new

                    height = self.setup.fragmentation_point[i].position.elev

        if height is not None:
            # self.angle_canvas.addItem(pg.InfiniteLine(pos=(height, 0), angle=90, pen=QColor(0, 0, 255)))
            self.height_plot.ax2.axvline(height/1000, c='b', label="Optimal Ballistic Solution")
            # self.height_canvas.addItem(pg.InfiniteLine(pos=(height, 0), angle=90, pen=QColor(0, 0, 255)))
        
        # self.angle_canvas.addItem(pg.InfiniteLine(pos=(0, 90), angle=0, pen=QColor(255, 0, 0)))

        ########################
        # Fit Hyperbola
        ########################

        try:
            x_vals = []
            y_vals = []

            for pp, point in enumerate(self.dots_y):
                if not np.isnan(point):
                    y_vals.append(point)
                    x_vals.append(self.dots_x[pp])

            x_vals, y_vals = np.array(x_vals), np.array(y_vals)

            def hypfunc(x, a, b, h, k):
                return b*np.sqrt(1 + ((x - h)/a)**2) + k
            def invhypfunc(x, a, b, h, k):
                return a*np.sqrt(((x - k)/b)**2 - 1) + h

            # Hyperbola in the form y = kx (since we invert the x values)
            from scipy.optimize import curve_fit

            popt, pcov = curve_fit(hypfunc, x_vals, y_vals)
         
            h = np.linspace(0, 100000, 1000)

            self.height_plot.ax2.plot(h/1000, hypfunc(h, *popt), c='m')
        except TypeError:
            print("Could not generate hyperbola fit!")
        except RuntimeError:
            print("Optimal hyperbola not found!")
        except ValueError:
            print("No Arrivals!")

        ########################
        # Overpressure vs Time
        ########################
        try:
            st, resp, gap_times = procStream(stn, ref_time=self.setup.fireball_datetime)
            st = findChn(st, self.channel)
            waveform_data, time_data = procTrace(st, ref_datetime=self.setup.fireball_datetime,\
                    resp=resp, bandpass=[2, 8])



        except ValueError:
            print("Could not filter waveform!")

        # waveform_peaks = []
        # time_peaks = []
        # N = 10

        # for wave, time_pack in zip(waveform_data, time_data):
        #     for ii in range(0, len(wave), N):
        #         peak = np.nanmax(np.abs(wave[ii:ii+N]))

        #         time_peak = np.mean(time_pack[ii:ii+N] - nom_pick.time)

        #         waveform_peaks.append(peak)
        #         time_peaks.append(time_peak)

        # ##################
        # # Relate Time to Height
        # ##################

        # height_peaks = invhypfunc(time_peaks, *popt)

        # ##################
        # # Get Bounds for Heights
        # ##################
        # min_height = 17000
        # max_height = 40000
        # h_indicies = np.where(np.logical_and(height_peaks>=min_height, height_peaks<=max_height))

        # new_heights = []
        # new_press = []
        # for hh in h_indicies[0]:

        #     new_heights.append(height_peaks[hh])
        #     new_press.append(waveform_peaks[hh])

        # plt.subplot(3, 1, 1)
    
        # plt.plot(new_heights, new_press)
        # plt.xlabel("Height [m]")
        # plt.ylabel("Overpressure [Pa]")


        # spacer = 5*60
        # shifter = 10
        # periods = []
        # period_times = []

        # for wave, time_pack in zip(waveform_data, time_data):
        #     for ii in range(len(wave)//shifter):

        #         temp_waveform = wave[shifter*ii:shifter*ii+spacer]
        #         temp_time = time_pack[shifter*ii:shifter*ii+spacer]

        #         period = findDominantPeriodPSD(temp_waveform, st.stats.sampling_rate)

        #         periods.append(period)
        #         period_times.append(time_pack[shifter*ii])

        # plt.subplot(3, 1, 2)
        # height_periods = invhypfunc(period_times, *popt)
        # h_indicies = np.where(np.logical_and(height_periods>=min_height, height_periods<=max_height))

        # new_heights = []
        # new_period = []
        # for hh in h_indicies[0]:

        #     new_heights.append(height_periods[hh])
        #     new_period.append(periods[hh])


        # plt.plot(new_heights, new_period)
        # plt.xlabel("Height [m]")
        # plt.ylabel("Dominant Period [s]")        


        # plt.subplot(3, 1, 3)
        # ws_list = []
        # lin_list = []
        # h_list = []
        # N_times=1
        # ws_list_p = []
        # lin_list_p = []
        # h_list_p = []
        # for ii in range(len(new_heights)//N_times):
        #     ws, lin = self.periodSearch(new_heights[ii*N_times], stn, new_period[ii*N_times])
        #     print("Phase 1: {:}/{:}".format(ii+1, len(new_heights)//N_times))
        #     ws_list_p.append(ws)
        #     lin_list_p.append(lin)
        #     h_list_p.append(new_heights[ii*N_times])
        # plt.scatter(h_list_p, ws_list_p, label='Period Weak-Shock', alpha=0.4)
        # plt.scatter(h_list_p,lin_list_p, label='Period Linear', alpha=0.4)
        # # for ii in range(len(new_heights)//N_times):
        # #     ws, lin = self.presSearch(new_heights[ii*N_times], stn, new_press[ii*N_times])
        # #     print("Phase 2: {:}/{:}".format(ii+1, len(new_heights)//N_times))
        # #     ws_list.append(ws)
        # #     lin_list.append(lin)
        # #     h_list.append(new_heights[ii*N_times])
        # # plt.scatter(h_list, ws_list, label='Pressure Weak-Shock', alpha=0.4)
        # # plt.scatter(h_list,lin_list, label='Pressure Linear', alpha=0.4)

        # plt.xlabel("Height [m]")
        # plt.ylabel("Relaxation Radius [m]")
        # plt.legend()
        # plt.show()   

        # fit_hyper = pg.PlotDataItem(x=h, y=hypfunc(h, *popt))
        # self.height_canvas.addItem(fit_hyper, update=True, pen='m')
        
        #####################
        # Build plot window
        #####################

        # 25 deg tolerance window
        # phigh = pg.PlotCurveItem([np.nanmin(X), np.nanmax(X)], [65, 65], pen = 'g')           
        # plow = pg.PlotCurveItem([np.nanmin(X), np.nanmax(X)], [115, 115], pen = 'g')                  
        # pfill = pg.FillBetweenItem(phigh, plow, brush = (0, 0, 255, 100))
        # self.angle_canvas.addItem(phigh)
        # self.angle_canvas.addItem(plow)
        # self.angle_canvas.addItem(pfill)

        # Build axes
        # self.height_canvas.setTitle('Fragmentation Height Prediction of Given Pick', color=(0, 0, 0))
        # self.angle_canvas.setTitle('Angles of Initial Acoustic Wave Path', color=(0, 0, 0))
        # self.height_canvas.setLabel('left', 'Difference in Time from {:.2f}'.format(nom_pick.time), units='s')
        # self.angle_canvas.setLabel('left', 'Angle Away from Trajectory of Initial Acoustic Wave Path', units='deg')
        # self.height_canvas.setLabel('bottom', 'Height of Solution', units='m')
        # self.angle_canvas.setLabel('bottom', 'Height of Solution', units='m')
        tol = 5 #seconds

        self.height_plot.ax2.set_xlim((-10, 100))
        self.height_plot.ax2.set_ylim((np.min([0, np.nanmin(self.dots_y)]) - tol, tol + np.max([0, np.nanmax(self.dots_y)])))
        # print(len(X), len(Y))
        # self.height_plot.ax2.scatter(np.array(X)/1000, np.array(Y), c='m')

        self.height_plot.ax2.set_xlabel("Height [km]")
        self.height_plot.ax2.set_ylabel("Time after {:.2f} [s]".format(nom_pick.time))
        self.height_plot.ax2.legend()
        self.height_plot.ax2.grid(alpha=0.2)

        # self.height_plot.ax1.show()
        # self.height_plot.ax2.show()
        self.height_plot.figure.tight_layout()
        self.height_plot.figure.subplots_adjust(hspace=0)
        self.height_plot.show()
        # self.height_canvas.setLimits(xMin=B.elev, xMax=A.elev, yMin=-40, yMax=100, minXRange=1000, maxXRange=33000, minYRange=2, maxYRange=140)

        # Fonts
        # font= QFont()
        # font.setPixelSize(20)
        # self.height_canvas.getAxis("bottom").tickFont = font
        # self.height_canvas.getAxis("left").tickFont = font
        # self.height_canvas.getAxis('bottom').setPen((255, 255, 255)) 
        # self.height_canvas.getAxis('left').setPen((255, 255, 255))
        # self.angle_canvas.getAxis("bottom").tickFont = font
        # self.angle_canvas.getAxis("left").tickFont = font
        # self.angle_canvas.getAxis('bottom').setPen((255, 255, 255)) 
        # self.angle_canvas.getAxis('left').setPen((255, 255, 255))
        self.setLayout(main_layout)

    def buildGUI(self):
        self.setWindowTitle('Fragmentation Staff')

        app_icon = QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)
        
        theme(self)

    def export(self):

        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = dlg[0]

        exporter = pg.exporters.SVGExporter(self.plot_view.scene())

        file_name = file_name + '.svg'
        exporter.export(file_name)

    def obtainPerts(self, data, frag, pt=0):
        data_new = []

        for i in range(len(data[frag][1])):
            data_new.append(data[frag][1][i][pt])
        data, remove = chauvenet(data_new)

        return data, remove

    def linkSeis(self):

        if self.linked_seis:
            self.seis_canvas.setYLink(None)   
        
        else:
            self.seis_canvas.setYLink(self.height_canvas)   
            
    def presSearch(self, h, stn, op):


        traj = self.setup.trajectory


        source = traj.findGeo(h)

        self.sounding_pres = None
        source_list = [source.lat, source.lon, source.elev/1000]


        stat = stn
        stat_pos = stat.metadata.position
        stat = [stat_pos.lat, stat_pos.lon, stat_pos.elev/1000]

        v = traj.v/1000

        theta = 90 - traj.zenith.deg

        dphi = np.degrees(np.arctan2(stat_pos.lon - source.lon, stat_pos.lat - source.lat)) - traj.azimuth.deg

        # Switch 3 doesn't do anything in this version of overpressure.py
        sw = [1, 0, 1]
        
        lat =   [source.lat, stat_pos.lat]
        lon =   [source.lon, stat_pos.lon]
        elev =  [source.elev, stat_pos.elev]

        sounding, _ = self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev)

        pres = 10*101.325*np.exp(-0.00012*sounding[:, 0])
        
        sounding_pres = np.zeros((sounding.shape[0], sounding.shape[1]+1))
        sounding_pres[:, :-1] = sounding
        sounding_pres[:, -1] = pres
        sounding_pres[:, 1] -= 273.15

        sounding_pres = np.flip(sounding_pres, axis=0)

        gem_inputs = [source_list, stat, v, theta, dphi, sounding_pres, sw, True, True]
        Ro_ws, Ro_lin = presSearch(op, gem_inputs, paths=False)

        return Ro_ws, Ro_lin

            
    def periodSearch(self, h, stn, op):


        traj = self.setup.trajectory


        source = traj.findGeo(h)

        self.sounding_pres = None
        source_list = [source.lat, source.lon, source.elev/1000]


        stat = stn
        stat_pos = stat.metadata.position
        stat = [stat_pos.lat, stat_pos.lon, stat_pos.elev/1000]

        v = traj.v/1000

        theta = 90 - traj.zenith.deg

        dphi = np.degrees(np.arctan2(stat_pos.lon - source.lon, stat_pos.lat - source.lat)) - traj.azimuth.deg

        # Switch 3 doesn't do anything in this version of overpressure.py
        sw = [1, 0, 1]
        
        lat =   [source.lat, stat_pos.lat]
        lon =   [source.lon, stat_pos.lon]
        elev =  [source.elev, stat_pos.elev]

        sounding, _ = self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev)

        pres = 10*101.325*np.exp(-0.00012*sounding[:, 0])
        
        sounding_pres = np.zeros((sounding.shape[0], sounding.shape[1]+1))
        sounding_pres[:, :-1] = sounding
        sounding_pres[:, -1] = pres
        sounding_pres[:, 1] -= 273.15

        sounding_pres = np.flip(sounding_pres, axis=0)

        gem_inputs = [source_list, stat, v, theta, dphi, sounding_pres, sw, True, True]
        Ro_ws, Ro_lin = periodSearch(op, gem_inputs, paths=False)

        return Ro_ws, Ro_lin

if __name__ == '__main__':

    pass