import obspy
import sys
import os
import datetime

import numpy as np
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('QT5Agg')

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT

from supra.GUI.Tools.Theme import *
from supra.GUI.Tools.GUITools import *
from supra.Stations.ProcessStation import *

N_ELEMENTS = 4

CHN = "BDF"

class ArrayStacker(QMainWindow):

    def __init__(self):
        super().__init__()

        self.setWindowTitle('mSeed Reader (not working)')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('images', 'bam.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        theme(self)

        self.buildGUI()

    def buildGUI(self):

        self._main = QWidget()
        self.setCentralWidget(self._main)
        layout = QGridLayout(self._main)

        self.raw_traces = [None]*N_ELEMENTS

        self.stat_graph = MatplotlibPyQT()
        self.sum_graph = MatplotlibPyQT()
        self.mseed_browser_label = [None]*N_ELEMENTS  
        self.mseed_browser_edits = [None]*N_ELEMENTS 
        self.mseed_browser_buton = [None]*N_ELEMENTS
        self.stat_shifter_label = [None]*N_ELEMENTS
        self.stat_shifter_edits = [None]*N_ELEMENTS

        self.stat_graph.ax = []
        for i in range(N_ELEMENTS):
            if i == 0:
                self.stat_graph.ax.append(self.stat_graph.figure.add_subplot(N_ELEMENTS, 1, i+1))
            else:

                self.stat_graph.ax.append(self.stat_graph.figure.add_subplot(N_ELEMENTS, 1, i+1, sharex=self.stat_graph.ax[0]))

        layout.addWidget(self.stat_graph, 0, 4, N_ELEMENTS*4, 1)


        self.sum_graph.ax = self.sum_graph.figure.add_subplot(1, 1, 1, sharex=self.stat_graph.ax[0])
        layout.addWidget(self.sum_graph, N_ELEMENTS*4 + 1, 4)


        for ii in range(N_ELEMENTS):
            self.mseed_browser_label[ii], self.mseed_browser_edits[ii], self.mseed_browser_buton[ii] = createFileSearchObj('mSeed File: {:}'.format(ii+1), layout, 2*ii, width=1, h_shift=0)
            _, self.stat_shifter_edits[ii] = createLabelEditObj('Shift [s]', layout, 2*ii+1, width=1, h_shift=0, tool_tip='', validate='float', default_txt='0')

            self.mseed_browser_buton[ii].clicked.connect(partial(fileSearch, ['mSeed File (*.mseed)'], self.mseed_browser_edits[ii]))


        self.stackRaw_button = createButton("Raw Stack", layout, 2*N_ELEMENTS+2, 2, self.stackRaw)        
        self.stack_button = createButton("Stack", layout, 2*N_ELEMENTS+2, 1, self.stack)
        self.read = createButton("Read", layout, 2*N_ELEMENTS+2, 0, self.read)

    def resetGraphs(self):
        for ii in range(N_ELEMENTS):
            self.stat_graph.ax[ii].clear()
        self.sum_graph.ax.clear()

    def read(self):


        self.resetGraphs()

        for ii in range(N_ELEMENTS):
            stat_text = self.mseed_browser_edits[ii].text()
            if stat_text is not None and len(stat_text) > 0: 
                

                self.raw_traces[ii] = obspy.read(stat_text)


                tr = self.raw_traces[ii].copy().select(channel=CHN)[0]

                if ii == 0:
                    ref_time = tr.stats.starttime.datetime
                # tr.stats.starttime = tr.stats.starttime - float(self.stat_shifter_edits.text())

                y, x = procTrace(tr, ref_datetime=ref_time - datetime.timedelta(seconds=float(self.stat_shifter_edits[ii].text())), \
                                        resp=None, bandpass=[3, 9.5], backup=False)
                x = x[0]
                y = y[0]

                if ii == 3:
                    self.startpoint = x[0]
                # Bandpass
                # Add shifter

                self.stat_graph.ax[ii].plot(x, y)

        # self.stat_graph.ax.legend()
        self.stat_graph.show()

    def findCommonEnds(self):

        common_start = None
        common_end = None

        for ii in range(N_ELEMENTS):

            stat_text = self.mseed_browser_edits[ii].text()

            if stat_text is not None and len(stat_text) > 0: 
                st = obspy.read(stat_text)[0]

                start = st.stats.starttime
                end = st.stats.endtime

                if common_start is None:
                    common_start = start

                if common_end is None:
                    common_end = end

                if start - common_start > 0:
                    common_start = start

                if end - common_end < 0:
                    common_end = end

        return common_start, common_end

    def findTotalRange(self):
        common_start = None
        common_end = None

        for ii in range(N_ELEMENTS):

            stat_text = self.mseed_browser_edits[ii].text()

            if stat_text is not None and len(stat_text) > 0: 
                st = obspy.read(stat_text)[0]

                start = st.stats.starttime
                end = st.stats.endtime

                if common_start is None:
                    common_start = start

                if common_end is None:
                    common_end = end

                if start - common_start < 0:
                    common_start = start

                if end - common_end > 0:
                    common_end = end

        return common_start, common_end

    def stack(self):

        

        for ii in range(N_ELEMENTS):

            stat_text = self.mseed_browser_edits[ii].text()

            if stat_text is not None and len(stat_text) > 0: 
                st = obspy.read(stat_text)
                # print("Unshifted", st)
                # a = input("Shift?")


                st[0].stats.starttime = st[0].stats.starttime - datetime.timedelta(seconds=float(self.stat_shifter_edits[ii].text()))
                # print("Shifted", st)
                common_start, common_end = self.findCommonEnds()
                total_start, total_end = self.findTotalRange()
                if ii == 0: 
                    master_st = st
                    # common_start, common_end = obspy.core.utcdatetime.UTCDateTime("2015-01-07T01:00:58.526038Z"),\
                    #                             obspy.core.utcdatetime.UTCDateTime("2015-01-07T02:05:57.676038Z")
                    master_st[0].trim(starttime=common_start, endtime=common_end)
                else:

                    tr = st.select(channel=CHN)[0]
                    tr.trim(starttime=common_start, endtime=common_end)
                    master_st.append(tr)


        
        
        stacked_st = master_st.stack(npts_tol=100)
        stacked_st[0].stats.starttime = common_start
        tr = stacked_st[0]

        y, x = procTrace(tr, ref_datetime=tr.stats.starttime.datetime , \
                                        resp=None, bandpass=[2, 8], backup=False)
        x = x[0]
        y = y[0]


        self.sum_graph.ax.plot(x + self.startpoint, y)
        self.sum_graph.show()
        # print("stacked", stacked_st)
        file_name = saveFile("mseed", note="")
        stacked_st.write(file_name)  

    def stackRaw(self):

        common_start, common_end = self.findCommonEnds()
        print(common_start, common_end)
        for ii in range(N_ELEMENTS):

            stat_text = self.mseed_browser_edits[ii].text()

            if stat_text is not None and len(stat_text) > 0: 
                st = obspy.read(stat_text)

                
                

                if ii == 0:
                    master_st = st
                    master_st[0].trim(starttime=common_start, endtime=common_end)
                else:
                    tr = st.select(channel=CHN)[0]
                    tr.trim(starttime=common_start, endtime=common_end)
                    master_st.append(tr)



        stacked_st = np.sum(master_st, axis=0)
        stacked_st[0].stats.starttime = obspy.read(self.mseed_browser_edits[-1].text())[0].stats.starttime

        tr = stacked_st[0]

        y, x = procTrace(tr, ref_datetime=tr.stats.starttime.datetime , \
                                        resp=None, bandpass=[2, 8], backup=False)
        x = x[0]
        y = y[0]


        self.sum_graph.ax.plot(x + self.startpoint, y)
        self.sum_graph.show()
        # print("stacked", stacked_st)
        file_name = saveFile("mseed", note="")
        stacked_st.write(file_name)  

if __name__ == "__main__":
    app = QApplication(sys.argv)

    splash_pix = QPixmap(os.path.join('images', 'wmpl.png'))
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()

    app.processEvents()

    gui = ArrayStacker()

    # gui.showFullScreen()
    # gui.showMaximized()
    gui.show()

    splash.finish(gui)

    sys.exit(app.exec_())