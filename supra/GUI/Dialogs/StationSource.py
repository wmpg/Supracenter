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

from supra.Stations.StationObj import *

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
        self.btn.clicked.connect(partial(fileSearch, ['Mini SEED (*.mseed)', 'SAC File (*.sac)', "GCF (*.gcf)"], self.edit))
        self.edit.textChanged.connect(self.loadStn)

        _, self.mresp_browser_edits, self.mresp_browser_buton = createFileSearchObj('Response File: ', layout, 2, width=1, h_shift=0)
        self.mresp_browser_buton.clicked.connect(partial(fileSearch, ['Response File (*.XML)'], self.mresp_browser_edits))


        _, self.network = createLabelEditObj('Network', layout, 3)
        _, self.code = createLabelEditObj('Code', layout, 4)
        _, self.lat = createLabelEditObj('Latitude', layout, 5)
        _, self.lon = createLabelEditObj('Longitude', layout, 6)
        _, self.elev = createLabelEditObj('Elevation', layout, 7)
        _, self.name = createLabelEditObj('Name', layout, 8)

        self.save_btn = createButton('Save Station', layout, 9, 2, self.saveStn)

    def loadStn(self):


        st = obspy.read(self.edit.text())

        network = st[0].stats.network
        code = st[0].stats.station


        self.network.setText(network)
        self.code.setText(code)

    def saveStn(self):

        position = Position(float(self.lat.text()), float(self.lon.text()), float(self.elev.text()))

        st = obspy.read(self.edit.text())

        meta = Metadata(self.network.text(), self.code.text(), position, self.name.text())

        try:
            response = obspy.read_inventory(self.mresp_browser_edits.text())
        except FileNotFoundError:
            response = None

        stn = Station(meta, st, response=response)

        self.bam.stn_list.append(stn)

        self.close()
