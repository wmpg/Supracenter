
import os
import numpy as np
import pickle
import obspy
import scipy.signal
import scipy.odr 

from obspy.signal.polarization import polarization_analysis
from obspy.core.utcdatetime import UTCDateTime
from obspy.signal.filter import remez_fir
from obspy import read_inventory
from obspy.signal.rotate import rotate2zne

from functools import partial

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *

from supra.Files.SaveObjs import Prefs

from supra.Utils.Classes import Position, Color

from supra.Stations.Filters import *

from supra.Stations.StationObj import Polarization

from supra.Supracenter.finalangle import finalanglecheck
from supra.Supracenter.anglescanrev import anglescanrev

import matplotlib.pyplot as plt
from scipy.fft import fft

def jurkevicSum(z, n, e):

    N = np.min([len(z), len(n), len(e)])
    X = np.array([z[:N], n[:N], e[:N]])

    S = np.matmul(X, np.transpose(X))

    return S


def jurkevicPol(z, n, e):

    S = jurkevicSum(z, n, e)

    w, v = scipy.linalg.eig(S)
    
    u = np.array(v)[:,0]
    
    az = np.arctan2(u[1]*np.sign(u[0]), u[2]*np.sign(u[0]))

    return np.degrees(az)%360

def bandpassFunc(waveform_data, low, high, delta):
        # Init the butterworth bandpass filter
    butter_b, butter_a = butterworthBandpassFilter(low, high, \
        1.0/delta, order=2)


    # Filter the data
    waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(waveform_data))

    return waveform_data

def jurkevicWindows(z, n, e, stats, window_size=5, window_overlap=1.0, bandpass=[2, 8]):

    z = bandpassFunc(z, bandpass[0], bandpass[1], stats.delta)
    n = bandpassFunc(n, bandpass[0], bandpass[1], stats.delta)
    e = bandpassFunc(e, bandpass[0], bandpass[1], stats.delta)


    npts = len(z)

    pps = stats.sampling_rate

    num_pts_in_window = int(window_size*pps)

    movement = int(window_overlap*num_pts_in_window)

    start = 0

    az_list = []
    t = []

    # 0.2 from class, 1.0 from paper
    tukey_window = scipy.signal.tukey(num_pts_in_window, 0.2)

    while start + num_pts_in_window < npts:


        z_part = z[start:start+num_pts_in_window]*tukey_window
        n_part = n[start:start+num_pts_in_window]*tukey_window
        e_part = e[start:start+num_pts_in_window]*tukey_window

        az_list.append(jurkevicPol(z_part,\
                                    n_part,\
                                    e_part)%180)
        t.append(start/pps)
        start += movement


    return az_list, t

# class AzimuthWindow(QWidget):

#     def __init__(self):
#         QWidget.__init__(self)

#         self.buildGUI()

#     def buildGUI(self):
#         self.setWindowTitle('Particle Motion')

#         app_icon = QtGui.QIcon()
#         app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
#         self.setWindowIcon(app_icon)
        
#         p = self.palette()
#         p.setColor(self.backgroundRole(), Qt.black)
#         self.setPalette(p)

#         theme(self)

#         layout = QGridLayout()
#         self.setLayout(layout)

#         _, self.azimuthedits = createLabelEditObj("Azimuth", layout, 0, validate='float')
#         _, self.azimutherroredits = createLabelEditObj("Azimuth Error", layout, 1, validate='float')
#         self.close_button = createButton("Mark", layout, 2, 1, self.onClose, args=[])

#     def onClose(self):

#         self.azimuth = float(self.azimuthedits.text())
#         self.azimuth_err = float(self.azimutherroredits.text())
#         # self.close()

class ParticleMotion(QWidget):

    def __init__(self, stn_map_canvas, bam, stn, channel, t_arrival=0, group_no=0):

        QWidget.__init__(self)
        

        self.bam = bam
        self.stn = stn
        self.channel = channel
        self.t_arrival = t_arrival
        self.stn_map_canvas = stn_map_canvas
        self.group_no = group_no

        self.buildGUI()
        

    def buildGUI(self):

        self.setWindowTitle('Particle Motion')

        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)
        
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QHBoxLayout()
        self.setLayout(layout)
     
        left_pane = QGridLayout()
        right_pane = QGridLayout()
        rightest_pane = QHBoxLayout()

        layout.addLayout(left_pane)
        layout.addLayout(right_pane)
        layout.addLayout(rightest_pane)

        self.added = False
        ################
        # LEFT PANE
        ################

        stn_label = QLabel("Station: {:}-{:} {:} | Channel: {:}*".format(self.stn.metadata.network, \
                    self.stn.metadata.code, self.stn.metadata.name, self.channel[0:2]))

        self.chn_label = [QLabel(), QLabel(), QLabel()]

        self.zne_canvas = [None]*3

        self.zne_plot_view = pg.GraphicsLayoutWidget()
        self.zne_canvas[0] = self.zne_plot_view.addPlot()
        self.zne_plot_view.nextRow()
        self.zne_canvas[1] = self.zne_plot_view.addPlot()
        self.zne_plot_view.nextRow()        
        self.zne_canvas[2]= self.zne_plot_view.addPlot()

        self.zne_plot_view.setBackground((0,0,0))
        self.zne_canvas[0].getAxis('bottom').setPen((255, 255, 255)) 
        self.zne_canvas[0].getAxis('left').setPen((255, 255, 255))
        self.zne_canvas[1].getAxis('bottom').setPen((255, 255, 255)) 
        self.zne_canvas[1].getAxis('left').setPen((255, 255, 255)) 
        self.zne_canvas[2].getAxis('bottom').setPen((255, 255, 255)) 
        self.zne_canvas[2].getAxis('left').setPen((255, 255, 255))         

        self.zne_canvas[1].setXLink(self.zne_canvas[0])
        self.zne_canvas[2].setXLink(self.zne_canvas[0])

        self.zne_canvas[1].setYLink(self.zne_canvas[0])
        self.zne_canvas[2].setYLink(self.zne_canvas[0])

        left_pane.addWidget(stn_label, 0, 1)
        left_pane.addWidget(self.zne_plot_view, 1, 1, 3, 1)
        
        for i in range(len(self.chn_label)):
            left_pane.addWidget(self.chn_label[i], i+1, 0)

        _, self.low_edits  = createLabelEditObj("Low", left_pane, 5, width=1, h_shift=0, tool_tip='', validate='float', default_txt='2')
        _, self.high_edits = createLabelEditObj("High", left_pane, 6, width=1, h_shift=0, tool_tip='', validate='float', default_txt='8')

        self.low_edits.returnPressed.connect(self.selectorPlot)
        self.high_edits.returnPressed.connect(self.selectorPlot)


        self.load_button = createButton('Load', left_pane, 7, 0, self.loadButton)

        #################
        # RIGHT PANE
        #################

        self.particle_motion_view = pg.GraphicsLayoutWidget()
        self.particle_motion_canvas = self.particle_motion_view.addPlot()

        self.waveform_plot_view = pg.GraphicsLayoutWidget()
        self.waveform_canvas = self.waveform_plot_view.addPlot()

        self.waveform_plot_view.setBackground((0,0,0))
        self.waveform_canvas.getAxis('bottom').setPen((255, 255, 255)) 
        self.waveform_canvas.getAxis('left').setPen((255, 255, 255))

        self.waveform_fft_view = pg.GraphicsLayoutWidget()
        self.waveform_fft_canvas = self.waveform_fft_view.addPlot()

        self.waveform_fft_view.setBackground((0,0,0))
        self.waveform_fft_canvas.getAxis('bottom').setPen((255, 255, 255)) 
        self.waveform_fft_canvas.getAxis('left').setPen((255, 255, 255))

        self.waveform_hist_view = pg.GraphicsLayoutWidget()


        self.plot_type = QComboBox()
        self.plot_type.addItem('Azimuth')
        self.plot_type.addItem('Incidence')
        self.plot_type.addItem('Azimuth Window')
        self.plot_type.addItem('Jurkevic Azimuth')
        self.plot_type.addItem('Azimuth Colourmap')
        self.plot_type.currentIndexChanged.connect(self.selectorPlot)

        right_pane.addWidget(self.waveform_fft_view, 0, 0)
        self.waveform_fft_canvas.setLogMode(True, True)

        right_pane.addWidget(self.waveform_plot_view, 1, 0)
        self.waveform_canvas.setXLink(self.particle_motion_canvas)
        
        rightest_pane.addWidget(self.waveform_hist_view)

        right_pane.addWidget(self.particle_motion_view, 2, 0)
        right_pane.addWidget(self.plot_type, 3, 0)

        self.mark_button = createButton('Mark', right_pane, 4, 0, self.markWaveform)

        h_label = QLabel("Height")
        right_pane.addWidget(h_label, 5, 0)

        self.ray_height = QLineEdit("50000")
        right_pane.addWidget(self.ray_height, 5, 1)

        win_len_label = QLabel("Window Length")
        right_pane.addWidget(win_len_label, 6, 0)

        self.win_len_e = QLineEdit("2")
        right_pane.addWidget(self.win_len_e, 6, 1)

        win_frac_label = QLabel("Window Fraction")
        right_pane.addWidget(win_frac_label, 7, 0)

        self.win_frac_e = QLineEdit("0.05")
        right_pane.addWidget(self.win_frac_e, 7, 1)

        self.win_len_e.returnPressed.connect(self.selectorPlot)
        self.win_frac_e.returnPressed.connect(self.selectorPlot)

    def rotate2ZNE(self):
        """
        Rotates non ZNE stations to ZNE
        """

        orientation_list = []

        for i in range(len(self.condensed_stream)):
            orientation_list.append(self.condensed_stream[i].stats.channel[2])

        if "E" not in orientation_list and "N" not in orientation_list:

            stat_orientation = []
            stat_waveform = []
            stat_lens = []

            for i in range(len(self.condensed_stream)):

                stat_orientation.append(self.stn.response.get_orientation(self.condensed_stream[i].get_id()))
                stat_waveform.append(np.copy(self.condensed_stream[i].data))
                stat_lens.append(len(stat_waveform[i]))

            # wavelengths must all be the same length
            min_len = np.min(stat_lens)

            waveform_array = rotate2zne(stat_waveform[0][0:min_len], stat_orientation[0]['azimuth'], stat_orientation[0]['dip'],\
                                        stat_waveform[1][0:min_len], stat_orientation[1]['azimuth'], stat_orientation[1]['dip'],\
                                        stat_waveform[2][0:min_len], stat_orientation[2]['azimuth'], stat_orientation[2]['dip'])

            for i in range(len(self.condensed_stream)):

                self.condensed_stream[i].data = waveform_array[i]

    def loadButton(self):
        self.loadWaveform()

    def loadWaveform(self):
        
        font=QtGui.QFont()
        font.setPixelSize(24)

        self.selector = [None]*3
        self.hilberted = [None]*3

        self.condensed_stream = self.stn.stream.select(channel="{:}*".format(self.channel[0:2]))

        self.rotate2ZNE()

        for i in range(len(self.condensed_stream)):

            st = self.condensed_stream[i]
            stn = self.stn

            delta = st.stats.delta
            start_datetime = st.stats.starttime.datetime
            end_datetime = st.stats.endtime.datetime

            stn.offset = (start_datetime - self.bam.setup.fireball_datetime).total_seconds()

            self.current_waveform_delta = delta
            self.current_waveform_time = np.arange(0, st.stats.npts / st.stats.sampling_rate, \
                 delta)

            time_data = np.copy(self.current_waveform_time)
            
            st.detrend()

            resp = stn.response
            st2 = st.copy()
            st2 = st2.remove_response(inventory=resp, output="DISP")
            # st2.remove_sensitivity(resp) 

            waveform_data = st2.data

            waveform_data = waveform_data[:len(time_data)]
            time_data = time_data[:len(waveform_data)] + stn.offset

            self.current_waveform_processed = waveform_data

            #####
            # Bandpass
            #####
            low  = float(self.low_edits.text())
            high = float(self.high_edits.text())

            if low > 0 and high > 0:
                # Init the butterworth bandpass filter
                butter_b, butter_a = butterworthBandpassFilter(low, high, \
                    1.0/self.current_waveform_delta, order=2)

                # Filter the data
                waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(waveform_data))

            else:
                # Init the butterworth bandpass filter
                butter_b, butter_a = butterworthBandpassFilter(2, 8, \
                    1.0/self.current_waveform_delta, order=2)

                # Filter the data
                waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(waveform_data))                

            self.current_station_waveform = pg.PlotDataItem(x=time_data, y=waveform_data, pen='w')
            self.zne_canvas[i].addItem(self.current_station_waveform)
            self.hilberted[i] = np.abs(scipy.signal.hilbert(waveform_data))
            # self.zne_canvas[i].addItem(pg.PlotDataItem(x=time_data, y= self.hilberted[i], pen='r'))
            # self.zne_canvas[i].addItem(pg.PlotDataItem(x=time_data, y=-self.hilberted[i], pen='r'))
            self.zne_canvas[i].setXRange(self.t_arrival-100, self.t_arrival+100, padding=1)
            self.zne_canvas[i].setLabel('bottom', "Time after {:} s".format(self.bam.setup.fireball_datetime))
            self.zne_canvas[i].setLabel('left', "Displacement", units='m')

            self.zne_canvas[i].plot(x=[-10000, 10000], y=[0, 0], pen=pg.mkPen(color=(100, 100, 100)))

            self.chn_label[i].setText("{:}".format(st.stats.channel[2]))
            
            self.chn_label[i].setFont(font) 

            self.selector[i] = pg.LinearRegionItem()
            self.selector[i].setZValue(self.t_arrival)
            self.selector[i].sigRegionChanged.connect(partial(self.selectorChanged, i))
            self.selector[i].sigRegionChangeFinished.connect(partial(self.selectorChanged, i))
            self.selector[i].sigRegionChangeFinished.connect(self.selectorPlot)

            self.zne_canvas[i].addItem(self.selector[i])

    def makePropLine(self, D, alpha=255):

        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)
        lats = []
        lons = []
        for line in D:
            temp = Position(0, 0, 0)
            temp.x = line[0]
            temp.y = line[1]
            temp.z = line[2]
            temp.pos_geo(ref_pos)
            if not np.isnan(temp.lat) and not np.isnan(temp.lon):
                lats.append(temp.lat)
                lons.append(temp.lon)

        lats.sort()
        lons.sort()

        start_pt = pg.PlotCurveItem()
        start_pt.setData(x=lons, y=lats)
        start_pt.setPen((255, 85, 0, alpha))

        return start_pt

    def propegateBack(self, offset):

        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

        S = self.stn.metadata.position
        
        S.pos_loc(ref_pos)

        sounding, perturbations = self.bam.atmos.getSounding(lat=[S.lat, S.lat], lon=[S.lon, S.lon], heights=[S.elev, float(self.ray_height.text())])

        D = []
        # Use 25 angles between 90 and 180 deg

        ### To do this more right, calculate D with bad winds, and then use D to find new winds and then recalc D

        for zenith in np.linspace(1, 89, 25):
            # D = anglescanrev(S.xyz, self.azimuth + offset, zenith, sounding, wind=True)
            # D = anglescanrev(S.xyz, (self.azimuth + offset + 180)%360, zenith, sounding, wind=True)

            D.append(anglescanrev(S.xyz, self.azimuth + offset, zenith, sounding, wind=True))
            D.append(anglescanrev(S.xyz, (self.azimuth + offset + 180)%360, zenith, sounding, wind=True))
        # pt, err = finalanglecheck(self.bam, self.bam.setup.trajectory, self.stn.metadata.position, self.azimuth)

        start_pt = self.makePropLine(np.array(D))


        return start_pt

    def markWaveform(self):
        
        if self.group_no == 0:
            az_line = pg.InfiniteLine(pos=(self.stn.metadata.position.lon, \
                        self.stn.metadata.position.lat), angle=90-self.azimuth, pen=QColor(0, 255, 0))
        else:
            az_line = pg.InfiniteLine(pos=(self.stn.metadata.position.lon, \
                        self.stn.metadata.position.lat), angle=90-self.azimuth, pen=QColor(0, 0, 255))

        self.stn_map_canvas.addItem(az_line)
        
        # Nominal curve
        nom_curve = self.propegateBack(0)

        # Azimuthal Error Bounds
        ph = self.propegateBack(self.az_error)
        pl = self.propegateBack(-self.az_error)
        pfill = pg.FillBetweenItem(ph, pl, brush=(225, 85, 0, 100))
        self.stn_map_canvas.addItem(ph)
        self.stn_map_canvas.addItem(pl)
        self.stn_map_canvas.addItem(pfill)


        if not hasattr(self.stn, "polarization"):
            self.stn.polarization = Polarization()

        if self.plot_type.currentText() != 'Azimuth Colourmap':
            
            self.stn.polarization.azimuth.append(self.azimuth)
            self.stn.polarization.azimuth_error.append(self.az_error)
        
        else:

            az =     float(input("Please enter in azimuth value: "))
            az_err = float(input("Please enter in azimuth error: "))

            self.stn.polarization.azimuth.append(az)
            self.stn.polarization.azimuth_error.append(az_err)


        # Time is start of the region
        self.stn.polarization.time.append(self.selector[0].getRegion()[0])

        a = Color()
        self.stn.color = a.generate() 

        print(self.stn.polarization.azimuth, self.stn.polarization.azimuth_error)

        self.close()

    def drawWaveform(self):
        pass

    def selectorChanged(self, i):

        for j in range(3):
            self.selector[j].blockSignals(True)
            if j != i:
                self.selector[j].setRegion(self.selector[i].getRegion())
            self.selector[j].blockSignals(False)

    def rescalePlot(self):

        x_rng = self.particle_motion_canvas.getAxis('bottom').range
        y_rng = self.particle_motion_canvas.getAxis('left').range

        base = np.min([x_rng[0], y_rng[0]])
        top  = np.max([x_rng[1], y_rng[1]])

        self.particle_motion_canvas.setXRange(base, top, padding=0)
        self.particle_motion_canvas.setYRange(base, top, padding=0)

    def selectorPlot(self):
        
        roi = self.selector[0].getRegion()

        if 1 < np.abs(roi[1] - roi[0]) < 30 or \
             (1 < np.abs(roi[1] - roi[0]) and \
                (self.plot_type.currentText() == 'Azimuth Window' or \
                    self.plot_type.currentText() == 'Jurkevic Azimuth') or \
                    self.plot_type.currentText() == 'Azimuth Colourmap'):
            
            ### TEMP
            stime = UTCDateTime(self.bam.setup.fireball_datetime) + roi[0]
            etime = UTCDateTime(self.bam.setup.fireball_datetime) + roi[1]      

            win_len  = float(self.win_len_e.text())
            win_frac = float(self.win_frac_e.text())

            try:

                pol_res = polarization_analysis(self.condensed_stream, win_len, win_frac, \
                                                    float(self.low_edits.text()), float(self.high_edits.text()), stime, etime, adaptive=True)
            except ValueError:
                pol_res = {}
                pol_res['azimuth'] = np.nan
                pol_res['timestamp'] = np.nan
            except IndexError:
                pol_res = {}
                pol_res['azimuth'] = np.nan
                pol_res['timestamp'] = np.nan
                print("Too many indicies for array - Not sure on this error yet")

            self.particle_motion_canvas.clear()
            self.waveform_canvas.clear()
            self.waveform_fft_canvas.clear()

            st_data = [None]*3
            ti_data = [None]*3
            raw_data = [None]*3
            noise_data = [None]*3

            for i in range(len(self.condensed_stream)):

                st = self.condensed_stream[i].copy()
                stn = self.stn

                delta = st.stats.delta
                start_datetime = st.stats.starttime.datetime
                end_datetime = st.stats.endtime.datetime

                stn.offset = (start_datetime - self.bam.setup.fireball_datetime).total_seconds()

                self.current_waveform_delta = delta
                self.current_waveform_time = np.arange(0, st.stats.npts / st.stats.sampling_rate, \
                     delta)

                time_data = np.copy(self.current_waveform_time)
        
                st.detrend()

                resp = stn.response
                st2 = st.copy()
                st2 = st2.remove_response(inventory=resp, output="DISP")
                # st2.remove_sensitivity(resp) 

                waveform_data = st2.data

                waveform_data = waveform_data[:len(time_data)]
                time_data = time_data[:len(waveform_data)] + stn.offset

                self.current_waveform_processed = waveform_data

                number_of_pts_per_s = st.stats.sampling_rate

                len_of_region = roi[1] - roi[0]

                num_of_pts_in_roi = len_of_region*number_of_pts_per_s

                num_of_pts_in_offset = np.abs(number_of_pts_per_s*stn.offset)

                num_of_pts_to_roi = roi[0]*number_of_pts_per_s

                pt_0 = int(num_of_pts_in_offset + num_of_pts_to_roi)
                pt_1 = int(pt_0 + num_of_pts_in_roi)

                raw_data[i] = waveform_data[pt_0:pt_1]
                noise_data[i] = waveform_data[0:(pt_1-pt_0)]

                #####
                # Bandpass
                #####
                low  = float(self.low_edits.text())
                high = float(self.high_edits.text())

                if low > 0 and high > 0:
                    # Init the butterworth bandpass filter
                    butter_b, butter_a = butterworthBandpassFilter(low, high, \
                        1.0/self.current_waveform_delta, order=2)

                    # Filter the data
                    waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(waveform_data))


                # waveform_data = st.data

                # butter_b, butter_a = butterworthBandpassFilter(float(self.low_edits.text()), float(self.high_edits.text()), \
                # 1.0/self.current_waveform_delta, order=6)

                # # Filter the data
                # waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(waveform_data))

                st_data[i] = waveform_data[pt_0:pt_1]
                ti_data[i] = time_data[pt_0:pt_1]

            e_r, n_r, z_r = raw_data[0], raw_data[1], raw_data[2]
            e_n, n_n, z_n = noise_data[0], noise_data[1], noise_data[2]
            e, n, z = st_data[0], st_data[1], st_data[2]
            et, nt, zt = ti_data[0], ti_data[1], ti_data[2]        

            e = np.array(e)
            n = np.array(n)

            def fit_func(beta, x):
                return beta[0]*x


            data = scipy.odr.Data(e, n)
            model = scipy.odr.Model(fit_func)
            odr = scipy.odr.ODR(data, model, beta0=[1.0])
            out = odr.run()
            az_slope = out.beta[0]
            az_error = out.sd_beta[0]
            azimuth = np.arctan2(1.0, az_slope)
            az_error = np.degrees(1.0 / ((1.0 ** 2 + az_slope ** 2) * azimuth) * az_error)

            z = np.array(z)
            r = np.sqrt(n**2 + e**2)*np.cos(azimuth - np.arctan2(n, e))

            data = scipy.odr.Data(r, z)
            model = scipy.odr.Model(fit_func)
            odr = scipy.odr.ODR(data, model, beta0=[1.0])
            out = odr.run()
            in_slope = out.beta[0]
            in_error = out.sd_beta[0]

            x_h, y_h = np.abs(scipy.signal.hilbert(st_data[0])), np.abs(scipy.signal.hilbert(st_data[1]))

            data = scipy.odr.Data(x_h, y_h)
            model = scipy.odr.Model(fit_func)
            odr = scipy.odr.ODR(data, model, beta0=[1.0])
            out = odr.run()
            h_az_slope = out.beta[0]
            h_az_error = out.sd_beta[0]

            h_azimuth = np.arctan2(1.0, h_az_slope)
            h_az_error = np.degrees(1.0/((1.0**2 + h_az_slope**2)*h_azimuth)*h_az_error)

            incidence = np.arctan2(1.0, in_slope)

            in_error = np.degrees(1.0/((1.0**2 + in_slope**2)*incidence)*in_error)

            azimuth = np.degrees(azimuth)
            incidence = np.degrees(incidence)

            # Fit to lines using orthoganal distance regression to a linear function
            if self.group_no == 0:
                pen = QColor(0, 255, 0)
                brush = QColor(0, 255, 0, 125)
            else:
                pen = QColor(0, 0, 255)
                brush = QColor(0, 0, 255, 125)

            if self.plot_type.currentText() == 'Azimuth':
                p_mot_plot = pg.PlotDataItem()
                p_mot_plot.setData(x=e, y=n)
                self.particle_motion_canvas.addItem(pg.InfiniteLine(pos=(0, 0), angle=90-azimuth, pen=pen))
                self.particle_motion_canvas.setLabel('bottom', "Channel: {:}".format(self.condensed_stream[0].stats.channel))
                self.particle_motion_canvas.setLabel('left', "Channel: {:}".format(self.condensed_stream[1].stats.channel))

                self.particle_motion_canvas.addItem(p_mot_plot)
                self.particle_motion_canvas.addItem(pg.TextItem(text='Azimuth = {:.2f}° ± {:.2f}°'.format(azimuth, az_error), color=(255, 0, 0), \
                                        anchor=(0, 0)))
                self.particle_motion_canvas.setXRange(np.min(e), np.max(e), padding=0)
                self.particle_motion_canvas.setYRange(np.min(n), np.max(n), padding=0)


            elif self.plot_type.currentText() == 'Incidence':
                p_mot_plot = pg.PlotDataItem()
                p_mot_plot.setData(x=r, y=z)
                self.particle_motion_canvas.addItem(pg.InfiniteLine(pos=(0, 0), angle=90-incidence, pen=pen))
                self.particle_motion_canvas.setLabel('bottom', "Horizontal in Direction of Azimuth")
                self.particle_motion_canvas.setLabel('left', "Channel: {:}".format(self.condensed_stream[2].stats.channel))
                self.particle_motion_canvas.addItem(p_mot_plot)
                self.particle_motion_canvas.addItem(pg.TextItem(text='Incidence = {:.2f}° ± {:.2f}°'.format(incidence, in_error), color=(255, 0, 0), \
                                        anchor=(0, 0)))
                
                self.particle_motion_canvas.setXRange(np.min(r), np.max(r), padding=0)
                self.particle_motion_canvas.setYRange(np.min(z), np.max(z), padding=0)

            elif self.plot_type.currentText() == 'Jurkevic Azimuth':
                # fix this
                bndps = [[0.001, 10], [0.1, 20], [1, 30]]
                
                for bb, b in enumerate(bndps):
                    az_list, t_list = jurkevicWindows(z_r, n_r, e_r, st.stats, window_size=win_len, window_overlap=win_frac, bandpass=b)
                
                    p_mot_plot = pg.ScatterPlotItem()
                    p_mot_plot.setData(x=t_list, y=az_list, pen=pg.intColor(bb), \
                                brush=pg.intColor(bb), name="Bandpass: {:} - {:} Hz".format(b[0], b[1]))
                    self.particle_motion_canvas.addItem(p_mot_plot)

                    # print("{:} -> Bandpass: {:} - {:} Hz".format(pg.intColor(bb), b[0], b[1]))
                    

                self.particle_motion_canvas.setLabel('bottom', "Time")
                self.particle_motion_canvas.setLabel('left', "Azimuth")
                
                self.particle_motion_canvas.setXRange(0, np.max(t_list), padding=0)
                self.particle_motion_canvas.setYRange(0, 180, padding=0)
                self.particle_motion_canvas.setLimits(xMin=0, xMax=np.max(t_list), yMin=0, yMax=180)


            elif self.plot_type.currentText() == 'Azimuth Colourmap':

                bandpass_low = float(self.low_edits.text())
                bandpass_high = float(self.high_edits.text())

                bins = 15
                bin_overlap = 0.5


                # this equation is wrong
                size_of_bin = ((bandpass_high - bandpass_low) - (bins - 1)*bin_overlap)/bins

                lows  = np.linspace(bandpass_low, bandpass_high-size_of_bin, bins)
                highs = np.linspace(bandpass_low+size_of_bin, bandpass_high, bins)



                bndps = []

                for i in range(len(lows)):
                    bndps.append([lows[i], highs[i]])

                total_list = []
                for bb, b in enumerate(bndps):
                    az_list, t_list = jurkevicWindows(z_r, n_r, e_r, st.stats, window_size=win_len, window_overlap=1.0, bandpass=b)
                    total_list.append(az_list)

                img = pg.ImageItem()
                self.particle_motion_canvas.addItem(img)

                az_img = np.array(total_list)
                img.setImage(np.transpose(az_img))
                img.scale(t_list[-1]/np.size(az_img, axis=1),\
                            lows[-1]/np.size(az_img, axis=0))
                img.translate(0, bandpass_low)
                self.hist = pg.HistogramLUTItem()

                # Link the histogram to the image
                self.hist.setImageItem(img)
                if self.added:
                    self.waveform_hist_view.clear()
                    

                self.waveform_hist_view.addItem(self.hist)    
                self.added = True
                # Fit the min and max levels of the histogram to the data available
                self.hist.setLevels(np.min(az_img), np.max(az_img))
                # This gradient is roughly comparable to the gradient used by Matplotlib
                # You can adjust it and then save it using hist.gradient.saveState()
                self.hist.gradient.restoreState(
                        {'mode': 'rgb',
                         'ticks': [(0.5, (0, 182, 188, 255)),
                                   (1.0, (246, 111, 0, 255)),
                                   (0.0, (75, 0, 113, 255))]})
                self.particle_motion_canvas.setLimits(xMin=0, xMax=t_list[-1]+bandpass_low, yMin=bandpass_low, yMax=bandpass_high)
                self.hist.setHistogramRange(0, 180)
                # self.particle_motion_canvas.setLookupTable(self.hist)

            else:
                az_window = pol_res['azimuth']

                # times are in unix time, casting UTCDatetime to float makes it also unix
                t_window = pol_res['timestamp'] - float(stime)
                az_win_error = pol_res['azimuth_error']

                p_mot_plot = pg.ScatterPlotItem()
                p_mot_plot.setData(x=t_window, y=az_window)
                p_mot_plot_err = pg.ErrorBarItem()
                p_mot_plot_err.setData(x=t_window, y=az_window, height=az_win_error)
                azimuth = np.mean(az_window)
                az_error = np.std(az_window)
                self.particle_motion_canvas.setLabel('bottom', "Time")
                self.particle_motion_canvas.setLabel('left', "Azimuth")
                self.particle_motion_canvas.addItem(pg.InfiniteLine(pos=(0, azimuth), angle=0, pen=pen))
                self.particle_motion_canvas.addItem(pg.TextItem(text='Azimuth = {:.2f}° ± {:.2f}°'.format(azimuth, az_error), color=(255, 0, 0), \
                        anchor=(0, 0)))
                self.particle_motion_canvas.addItem(p_mot_plot_err)
                self.particle_motion_canvas.addItem(p_mot_plot)
                self.particle_motion_canvas.setXRange(0, np.max(t_window), padding=0)
                self.particle_motion_canvas.setYRange(0, 180, padding=0)
                # self.particle_motion_canvas.getViewBox().setLimits(xMin=0, xMax=15,   
                #              yMin=range_[1][0], yMax=range_[1][1]) 
                


            self.azimuth = azimuth
            self.az_error = az_error

            roi_waveform = pg.PlotDataItem()
            pts = np.linspace(pt_0/st.stats.sampling_rate + stn.offset - roi[0], \
                                pt_1/st.stats.sampling_rate + stn.offset - roi[0], num=len(z))
            roi_waveform.setData(x=pts, y=bandpassFunc(z, 2, 8, st.stats.delta))            
            self.waveform_canvas.addItem(roi_waveform)

            
            sps = st.stats.sampling_rate
            dt = 1/st.stats.sampling_rate
            length = len(z_r)
            freq = np.linspace(1/len_of_region, (sps/2), length)*sps/length
            
            FAS = abs(fft(z_r))
            FAS_n = abs(fft(z_n))
            fas_data = pg.PlotDataItem()
            fas_data.setData(x=freq, y=FAS, pen=(255, 0, 0))

            fas_noise_data = pg.PlotDataItem()
            fas_noise_data.setData(x=freq, y=FAS_n, pen=(255, 255, 255))

            fas_diff_data = pg.PlotDataItem()
            fas_diff_data.setData(x=freq, y=np.abs(FAS/FAS_n), pen=(0, 125, 255))

            print("Red - Data", "White - Noise", "Blue - Data/Noise")
            self.waveform_fft_canvas.addItem(fas_data)
            self.waveform_fft_canvas.addItem(fas_noise_data)
            self.waveform_fft_canvas.addItem(fas_diff_data)

            # self.rescalePlot()
if __name__ == '__main__':

    pass