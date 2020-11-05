import os
import numpy as np
import pickle

from functools import partial

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *
from supra.Utils.TryObj import *

from supra.Stations.Filters import *


class BandpassWindow(QWidget):

    def __init__(self, bam, stn, channel, t_arrival=0):

        self.bam = bam
        self.stn = stn
        self.channel = channel
        self.t_arrival = t_arrival

        QWidget.__init__(self)
        self.buildGUI()


    def buildGUI(self):

        self.setWindowTitle('Bandpass Optimizer')

        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)
        
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)
     
        self.bandpass_view = pg.GraphicsLayoutWidget()
        self.bandpass_canvas = self.bandpass_view.addPlot()

        self.bp_button = createButton('Bandpass', layout, 2, 2, self.bandpass)

        layout.addWidget(self.bandpass_view, 1, 1, 1, 2)

        self.stream = self.stn.stream.select(channel="{:}".format(self.channel))

        st = self.stn.stream.select(channel="{:}".format(self.channel))[0]
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
        st = st.remove_response(inventory=resp, output="DISP")
        st.remove_sensitivity(resp) 

        waveform_data = st.data

        self.orig_data = np.copy(waveform_data)

        waveform_data = waveform_data[:len(time_data)]
        time_data = time_data[:len(waveform_data)] + stn.offset

        self.current_waveform_processed = waveform_data

        # Init the butterworth bandpass filter
        butter_b, butter_a = butterworthBandpassFilter(2, 8, \
            1.0/self.current_waveform_delta, order=6)

        # Filter the data
        waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(waveform_data))

        self.current_station_waveform = pg.PlotDataItem(x=time_data, y=waveform_data, pen='w')
        self.bandpass_canvas.addItem(self.current_station_waveform)
        self.bandpass_canvas.setXRange(self.t_arrival-100, self.t_arrival+100, padding=1)
        self.bandpass_canvas.setLabel('bottom', "Time after {:} s".format(self.bam.setup.fireball_datetime))
        self.bandpass_canvas.setLabel('left', "Signal Response")

        self.bandpass_canvas.plot(x=[-10000, 10000], y=[0, 0], pen=pg.mkPen(color=(100, 100, 100)))

        self.noise_selector = pg.LinearRegionItem(values=[0, 10], brush=(255, 0, 0, 100))

        self.signal_selector = pg.LinearRegionItem(values=[200, 210], brush=(0, 255, 0, 100))

        self.bandpass_canvas.addItem(self.noise_selector)
        self.bandpass_canvas.addItem(self.signal_selector)

    def determineROIidx(self, roi):
        
        len_of_region = roi[1] - roi[0]
        
        st = self.stn.stream.select(channel="{:}".format(self.channel))[0]
        
        number_of_pts_per_s = st.stats.sampling_rate
        num_of_pts_in_roi = len_of_region*number_of_pts_per_s

        num_of_pts_in_offset = np.abs(number_of_pts_per_s*self.stn.offset)

        num_of_pts_to_roi = roi[0]*number_of_pts_per_s

        pt_0 = int(num_of_pts_in_offset + num_of_pts_to_roi)
        pt_1 = int(pt_0 + num_of_pts_in_roi)

        return pt_0, pt_1

    def bandpass(self):

        st = self.stn.stream.select(channel="{:}".format(self.channel))[0]

        noise_roi = self.noise_selector.getRegion()
        signal_roi = self.signal_selector.getRegion()

        noise_a, noise_b = self.determineROIidx(noise_roi)
        signal_a, signal_b = self.determineROIidx(signal_roi)

        S = 1
        N = 1
        S_N = S/N
        result = [np.nan, np.nan, np.nan]

        fs = 1.0/self.current_waveform_delta
        # 0.5*fs is the maximum possible filtering
        # 0 is the lowest

        low_pass =  np.linspace(0.01, 0.5*fs, 100, endpoint=False)
        high_pass = np.linspace(0.01, 0.5*fs, 100, endpoint=False)
        SmN_last = 0
        for l in low_pass:
            for h in high_pass:

                if h > l + 1:
                
                    waveform_data = np.copy(self.orig_data)

                    # Init the butterworth bandpass filter
                    butter_b, butter_a = butterworthBandpassFilter(l, h, \
                        fs, order=6)

                    # Filter the data
                    waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(waveform_data))

                    S = np.mean(np.abs(waveform_data[signal_a:signal_b]))
                    N = np.mean(np.abs(waveform_data[noise_a:noise_b]))

                    S_N_last = S_N
                    S_N = S/N

                    if S_N > S_N_last and S - N > SmN_last:
                        result = [S_N, l, h]
                        SmN_last = S - N

        print("Signal/Noise:     {:.2f}".format(result[0]))
        print("Optimal Lowpass:  {:.2f}".format(result[1]))
        print("Optimal Highpass: {:.2f}".format(result[2]))