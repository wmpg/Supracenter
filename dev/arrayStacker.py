import obspy
import sys
import os

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

N_ELEMENTS = 8
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

        self.stat_graph.ax = self.stat_graph.figure.add_subplot(111)
        layout.addWidget(self.stat_graph, 0, 4, N_ELEMENTS, 1)

        self.sum_graph.ax = self.sum_graph.figure.add_subplot(111)
        layout.addWidget(self.sum_graph, N_ELEMENTS, 4)
        for ii in range(N_ELEMENTS):
            self.mseed_browser_label[ii], self.mseed_browser_edits[ii], self.mseed_browser_buton[ii] = createFileSearchObj('mSeed File: {:}'.format(ii+1), layout, 2*ii, width=1, h_shift=0)
            _, self.stat_shifter_edits = createLabelEditObj('Shift [s]', layout, 2*ii+1, width=1, h_shift=0, tool_tip='', validate='float', default_txt='0')

            self.mseed_browser_buton[ii].clicked.connect(partial(fileSearch, ['mSeed File (*.mseed)'], self.mseed_browser_edits[ii]))

            
        self.stack = createButton("Stack", layout, 2*N_ELEMENTS+2, 1, self.stack)
        self.read = createButton("Read", layout, 2*N_ELEMENTS+2, 0, self.read)

    def read(self):

        self.stat_graph.ax.clear()
        self.sum_graph.ax.clear()

        for ii in range(N_ELEMENTS):
            stat_text = self.mseed_browser_edits[ii].text()
            if stat_text is not None and len(stat_text) > 0: 
                

                self.raw_traces[ii] = obspy.read(stat_text)
                tr = self.raw_traces[ii].copy().select(channel=CHN)[0]

                if ii == 0:
                    ref_time = tr.stats.starttime
                # tr.stats.starttime = tr.stats.starttime - float(self.stat_shifter_edits.text())
                x = tr.times(reftime=ref_time - float(self.stat_shifter_edits.text()))
                y = tr.data
                # Bandpass
                # Add shifter
                print(tr.stats.starttime)
                print(x)

                self.stat_graph.ax.plot(x, y, label="Element {:}".format(ii + 1), alpha=0.2)

        self.stat_graph.ax.legend()
        self.stat_graph.show()


    def stack(self):

        for ii in range(N_ELEMENTS):

            stat_text = self.mseed_browser_edits[ii].text()

            if stat_text is not None and len(stat_text) > 0: 
                st = obspy.read(stat_text)
                # print("Unshifted", st)
                # a = input("Shift?")


                # st[0].stats.starttime = st[0].stats.starttime - float(a)
                # print("Shifted", st)

                if ii == 0: 
                    master_st = st
                    common_start, common_end = obspy.core.utcdatetime.UTCDateTime("2015-01-07T01:00:58.526038Z"),\
                                                obspy.core.utcdatetime.UTCDateTime("2015-01-07T02:05:57.676038Z")
                    master_st[0].trim(starttime=common_start, endtime=common_end)
                else:

                    tr = st.select(channel=CHN)[0]
                    tr.trim(starttime=common_start, endtime=common_end)
                    master_st.append(tr)

        print("Master", master_st)


        stacked_st = master_st.stack(stack_type=('pw', 2), npts_tol=0, time_tol=0)
        stacked_st[0].stats.starttime = obspy.core.utcdatetime.UTCDateTime("2015-01-07T01:00:58.526038Z")


        print("stacked", stacked_st)
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