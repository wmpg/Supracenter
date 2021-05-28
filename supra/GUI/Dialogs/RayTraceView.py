import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from functools import partial

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *


from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT
from supra.Stations.ProcessStation import procTrace, procStream, findChn, findDominantPeriodPSD, genFFT
from supra.Utils.Formatting import *
from supra.Geminus.geminusSearch import periodSearch, presSearch



class rtvWindowDialog(QWidget):

    def __init__(self, bam):

        QWidget.__init__(self)
        
        self.bam = bam
        self.buildGUI()
       



    def buildGUI(self):
        self.setWindowTitle('Ray-Trace Viewer')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)

        self.rtv_graph = MatplotlibPyQT()
      
        self.rtv_graph.ax = self.rtv_graph.figure.add_subplot(111)

        layout.addWidget(self.rtv_graph, 1, 1, 100, 1)