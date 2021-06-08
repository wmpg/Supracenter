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
from supra.Stations.ProcessStation import procTrace, procStream, findChn, findDominantPeriodPSD, genFFT
from supra.Files.SaveLoad import save
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

        self.bp_button = createButton('Bandpass', layout, 4, 2, self.bandpass)
        self.save_button = createButton('Save', layout, 4, 3, self.bandpassSave)

        layout.addWidget(self.bandpass_view, 1, 1, 1, 2)

        _, self.low_bandpass_edits = createLabelEditObj('Low Bandpass', layout, 2, width=1, h_shift=0, tool_tip='', validate='float', default_txt='2')
        _, self.high_bandpass_edits  = createLabelEditObj('High Bandpass', layout, 3, width=1, h_shift=0, tool_tip='', validate='float', default_txt='8')

        self.stream = self.stn.stream.select(channel="{:}".format(self.channel))

        st = self.stn.stream.select(channel="{:}".format(self.channel))[0]
        self.orig_trace = st.copy()
        stn = self.stn

        delta = self.orig_trace.stats.delta
        start_datetime = self.orig_trace.stats.starttime.datetime
        end_datetime = self.orig_trace.stats.endtime.datetime

        stn.offset = (start_datetime - self.bam.setup.fireball_datetime).total_seconds()

        self.current_waveform_delta = delta
        self.current_waveform_time = np.arange(0, self.orig_trace.stats.npts / self.orig_trace.stats.sampling_rate, \
             delta)

        time_data = np.copy(self.current_waveform_time)
        
        self.orig_trace.detrend()

        resp = stn.response
        self.orig_trace = self.orig_trace.remove_response(inventory=resp, output="DISP")
        # st.remove_sensitivity(resp) 

        waveform_data = self.orig_trace.data

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

    def bandpassSave(self):
        bandpass = [float(self.low_bandpass_edits.text()), float(self.high_bandpass_edits.text())]
        self.stn.bandpass = bandpass
        save(self.bam, file_check=False)

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


        waveform_data, t = procTrace(self.orig_trace, ref_datetime=self.bam.setup.fireball_datetime, \
                    resp=self.stn.response, bandpass=None, backup=False)

        S = waveform_data[0][signal_a:signal_b]
        N = waveform_data[0][noise_a:noise_b]

        # Make sure the windows are the same length
        if len(N) >= len(S):
            N = N[:len(S)]

        freq, FAS_S = genFFT(S, self.orig_trace.stats.sampling_rate)
        freq, FAS_N = genFFT(N, self.orig_trace.stats.sampling_rate)

        S_N_FAS = FAS_S/FAS_N

        plt.loglog(freq, FAS_S, label="Signal")
        plt.loglog(freq, FAS_N, label="Noise")
        plt.loglog(freq, S_N_FAS, label="Signal/Noise")
        plt.axhline(y=1)

        plt.legend()
        plt.show()