import os
import numpy as np
import pickle
import obspy

from functools import partial

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *

from supra.Files.SaveObjs import Prefs

from supra.Utils.Classes import Position


class StationWindow(QWidget):

    def __init__(self, bam):

        QWidget.__init__(self)
        
        self.bam = bam

        self.buildGUI()

    def buildGUI(self):
        self.setWindowTitle('Add Station From File')
        
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)
     
        _, self.edit, self.btn = createFileSearchObj('Source File', layout, 1, width=1, h_shift=0, tool_tip='')
        self.btn.clicked.connect(partial(fileSearch, ['Mini SEED (*.mseed)'], self.edit))
        self.edit.textChanged.connect(self.loadStn)

        _, self.network = createLabelEditObj('Network', layout, 2)
        _, self.code = createLabelEditObj('Code', layout, 3)
        _, self.lat = createLabelEditObj('Latitude', layout, 4)
        _, self.lon = createLabelEditObj('Longitude', layout, 5)
        _, self.elev = createLabelEditObj('Elevation', layout, 6)
        _, self.name = createLabelEditObj('Name', layout, 7)

        self.save_btn = createButton('Save Station', layout, 8, 2, self.saveStn)

    def loadStn(self):

        st = obspy.read(self.edit.text())
        network = st[0].stats.network
        code = st[0].stats.station

        self.network.setText(network)
        self.code.setText(code)

    def saveStn(self):

        position = Position(float(self.lat.text()), float(self.lon.text()), float(self.elev.text()))

        st = obspy.read(self.edit.text())

        stn = Station(self.network.text(), self.code.text(), position, st, name=self.name.text())

        self.bam.stn_list.append(stn)

