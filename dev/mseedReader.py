import obspy
import sys
import os

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg

from supra.GUI.Tools.Theme import *
from supra.GUI.Tools.GUITools import *

class MseedReader(QMainWindow):

    def __init__(self):
        super().__init__()


        self.setWindowTitle('mSeed Reader')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('images', 'bam.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        theme(self)

        self.buildGUI()

    def buildGUI(self):

        self._main = QWidget()
        self.setCentralWidget(self._main)
        layout = QGridLayout(self._main)

        self.mseed_browser_label, self.mseed_browser_edits, self.mseed_browser_buton = createFileSearchObj('mSeed File: ', layout, 1, width=1, h_shift=0)
        self.mseed_browser_buton.clicked.connect(partial(fileSearch, ['mSeed File (*.mseed)'], self.mseed_browser_edits))

        self.mresp_browser_label, self.mresp_browser_edits, self.mresp_browser_buton = createFileSearchObj('Resp File: ', layout, 2, width=1, h_shift=0)
        self.mresp_browser_buton.clicked.connect(partial(fileSearch, ['Resp File (*.XML)'], self.mresp_browser_edits))

        self.print_stream_button = createButton("Print Stream", layout, 3, 1, self.printStream)
        self.print_resp_button = createButton("Print Response", layout, 3, 2, self.printResp)


    def printResp(self):

        station_resp_name = self.mresp_browser_edits.text()

        resp = obspy.read_inventory(station_resp_name)

        print(resp)

    def printStream(self):

        station_file_name = self.mseed_browser_edits.text()

        st = obspy.read(station_file_name)

        print(st)




if __name__ == "__main__":
    app = QApplication(sys.argv)

    splash_pix = QPixmap(os.path.join('images', 'wmpl.png'))
    splash = QSplashScreen(splash_pix, Qt.WindowStaysOnTopHint)
    splash.setMask(splash_pix.mask())
    splash.show()

    app.processEvents()

    gui = MseedReader()

    # gui.showFullScreen()
    # gui.showMaximized()
    gui.show()

    splash.finish(gui)

    sys.exit(app.exec_())