import os
import numpy as np
import pickle

from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from functools import partial
from obspy.core.utcdatetime import UTCDateTime
from obspy.signal.rotate import rotate_ne_rt
from supra.Stations.ProcessStation import procTrace, procStream, findChn

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *
from supra.Utils.Classes import *


class RotatePolWindow(QWidget):

    def __init__(self, stn, chn, lside, rside, ref_time):

        QWidget.__init__(self)
        
        self.ref_datetime = ref_time
        self.stn = stn
        self.chn = chn
        self.lside = lside
        self.rside = rside

        self.buildGUI()

    def buildGUI(self):
        self.setWindowTitle('Polarization from Rotation')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)

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

        self.backaz_view = pg.GraphicsLayoutWidget()
        self.backaz_canvas = self.backaz_view.addPlot()
        layout.addWidget(self.backaz_view, 2, 4, 1, 1)
        layout.addWidget(self.zne_plot_view, 1, 1, 3, 3)
        self.runsim = createButton("Run", layout, 1, 4, self.procGUI)


    def procGUI(self):

        # get stream
        st = self.stn.stream

        self.resp = self.stn.response

        # get selected channels
        chn_filter = self.chn[0:2]
        st = st.select(channel="{:}*".format(chn_filter))

        # Shorten time range of stream to times given
        st = st.trim(starttime=UTCDateTime(self.lside), endtime=UTCDateTime(self.rside))


        STEP_DEG = 1

        ba_list = np.arange(0, 360, STEP_DEG)
        r_max = []
        t_max = []
        ba_max = []
        for ba in ba_list:
            z = st.select(channel="**Z")[0]
            n = st.select(channel="**N")[0]
            e = st.select(channel="**E")[0]

            # filter traces
            z, times =  procTrace(z, ref_datetime=self.ref_datetime, resp=self.resp, bandpass=[2, 8])
            n, _ =      procTrace(n, ref_datetime=self.ref_datetime, resp=self.resp, bandpass=[2, 8])
            e, _ =      procTrace(e, ref_datetime=self.ref_datetime, resp=self.resp, bandpass=[2, 8])

            z = z[0]
            n = n[0]
            e = e[0]
            times = times[0]

            r, t = rotate_ne_rt(n, e, ba)

            # a = st.rotate('NE->RT', back_azimuth=ba)
            # print(a[0].stats.back_azimuth)
            print("Current Back-Azimuth = {:} deg | RAD = {:.2f} | TRN = {:.2f}".format(ba, np.max(r), np.max(t)))
            ba_max.append(ba)
            r_max.append(np.max(r))
            t_max.append(np.max(t))
            self.zne_canvas[0].plot(x=times, y=z, update=True, clear=True)
            self.zne_canvas[1].plot(x=times, y=r, update=True, clear=True)
            self.zne_canvas[2].plot(x=times, y=t, update=True, clear=True)
            rad = pg.PlotCurveItem(ba_max, r_max, pen=pg.mkPen(color=(255, 0, 0)))
            trn = pg.PlotCurveItem(ba_max, t_max, pen=pg.mkPen(color=(0, 0, 255)))           

            self.backaz_canvas.addItem(rad, update=True, clear=True, name="Radial")
            self.backaz_canvas.addItem(trn, update=True, clear=True, name="Transverse")
            self.backaz_canvas.setLabel('left', "Positive Maximum Amplitude (Radial)")
            self.backaz_canvas.setLabel('bottom', "Back-Azimuth [deg]")

            # self.backaz_canvas.plot(x=ba_max, y=t_max, update=True, clear=True, name="Transverse", pen=pg.mkPen(color=(0, 0, 255)))
            pg.QtGui.QApplication.processEvents()
