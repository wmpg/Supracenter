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

class StationList(QWidget):

    def __init__(self):

        QWidget.__init__(self)
        
        self.buildGUI()


    def buildGUI(self):
        self.setWindowTitle('Station Download Sources')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)
        
        self.s = [None]*20
        self.s[1]= createToggle('IRIS Data Management Center', layout, 1, width=1, h_shift=0, tool_tip='')
        self.s[2]= createToggle('Northern California Earthquake Data Center', layout, 2, width=1, h_shift=0, tool_tip='')
        self.s[3]= createToggle('Southern California Earthquake Data Center', layout, 3, width=1, h_shift=0, tool_tip='')
        self.s[4]= createToggle('TexNet – Texas Earthquake Data Center', layout, 4, width=1, h_shift=0, tool_tip='')
        self.s[5]= createToggle('BGR Hannover, Germany', layout, 5, width=1, h_shift=0, tool_tip='')
        self.s[6]= createToggle('Boğaziçi University, Kandilli Observatory', layout, 6, width=1, h_shift=0, tool_tip='')
        self.s[7]= createToggle('ETHZ', layout, 7, width=1, h_shift=0, tool_tip='')
        self.s[8]= createToggle('GEOFON Program, GFZ', layout, 8, width=1, h_shift=0, tool_tip='')
        self.s[9]= createToggle('ICGC', layout, 9, width=1, h_shift=0, tool_tip='')
        self.s[10] = createToggle('IPGP Data Center', layout, 10, width=1, h_shift=0, tool_tip='')
        self.s[11] = createToggle('INGV', layout, 11, width=1, h_shift=0, tool_tip='')
        self.s[12] = createToggle('LMU Munich, Germany', layout, 12, width=1, h_shift=0, tool_tip='')
        self.s[13] = createToggle('NIEP, Romania', layout, 13, width=1, h_shift=0, tool_tip='')
        self.s[14] = createToggle('NOA, Greece', layout, 14, width=1, h_shift=0, tool_tip='')
        self.s[15] = createToggle('ORFEUS Data Center', layout, 15, width=1, h_shift=0, tool_tip='')
        self.s[16] = createToggle('RESIF', layout, 16, width=1, h_shift=0, tool_tip='')
        self.s[17] = createToggle('USP Seismological Center, Brazil', layout, 17, width=1, h_shift=0, tool_tip='')
        self.s[18] = createToggle('Raspberry Shake, S.A', layout, 18, width=1, h_shift=0, tool_tip='')
        self.s[19] = createToggle('AusPass, Australia', layout, 19, width=1, h_shift=0, tool_tip='')

        try:
            with open(os.path.join('supra', 'Misc', 'BAMStationprefs.bam'), 'rb') as f:
                sources = pickle.load(f)
            for i in range(len(sources)):
                if sources[i] == 1:
                    self.s[i+1].setChecked(True)
        except:
            self.selectAll()   

        createButton('Select All', layout, 20, 1, self.selectAll)
        createButton('Select None', layout, 20, 2, self.selectNone)
        createButton('Save', layout, 20, 3, self.saveStats)

    def selectAll(self):
        for i in range(20):
            if i > 0:
                self.s[i].setChecked(True)

    def selectNone(self):
        for i in range(20):
            if i > 0:
                self.s[i].setChecked(False)

    def saveStats(self):
        a = []

        for i in range(20):
            if i > 0:
                if self.s[i].isChecked():
                    a.append(1)
                else:
                    a.append(0)

        with open(os.path.join('supra', 'Misc', 'BAMStationprefs.bam') , 'wb') as f:
            pickle.dump(a, f)

        self.close()