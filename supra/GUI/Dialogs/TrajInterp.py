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
from supra.Files.SaveLoad import save, loadSourcesIntoBam
from supra.Utils.Classes import Position, Supracenter

from supra.GUI.Dialogs.AddSource import Source 

class TrajInterpWindow(QWidget):

    def __init__(self, bam, obj):

        QWidget.__init__(self)
        
        self.obj = obj
        self.bam = bam
        self.setup = bam.setup
        self.buildGUI()


    def buildGUI(self):
        self.setWindowTitle('Trajectory Interpolator')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'preferences.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)
     
        _, self.high_point = createLabelEditObj('Maximum Height', layout, 0, width=1, h_shift=0, tool_tip='', validate='float', default_txt='50000')
        _, self.low_point  = createLabelEditObj('Minimum Height', layout, 1, width=1, h_shift=0, tool_tip='', validate='float', default_txt='17000')
        _, self.divisions  = createLabelEditObj('Divisions', layout, 2, width=1, h_shift=0, tool_tip='', validate='int', default_txt='250')

        save_button = QPushButton('Save to CSV')
        layout.addWidget(save_button, 3, 1)
        save_button.clicked.connect(self.trajInterp)

        save_to_file_button = QPushButton('Save to Event')
        layout.addWidget(save_to_file_button, 3, 2)
        save_to_file_button.clicked.connect(self.saveEvent)

    def trajInterp(self):

        traj = self.setup.trajectory

        
        try:
            points = traj.trajInterp2(div=tryInt(self.divisions.text()),\
                                      min_p=tryFloat(self.low_point.text()),\
                                      max_p=tryFloat(self.high_point.text()))
        except AttributeError as e:
            errorMessage("Unable to interpolate trajectory!", 2, info='Please define a trajectory source in the "sources" tab, save and then restart BAM', title='Yikes!', detail='{:}'.format(e))
            return None

        file_name = saveFile('csv')

        with open(file_name, 'w+') as f:
            for pt in points:
                f.write('{:}, {:}, {:}, {:} \n'.format(pt[0], pt[1], pt[2], pt[3]))


        self.close()

    def saveEvent(self):
        traj = self.setup.trajectory
        
        try:
            points = traj.trajInterp2(div=tryInt(self.divisions.text()),\
                                      min_p=tryFloat(self.low_point.text()),\
                                      max_p=tryFloat(self.high_point.text()))
        except AttributeError as e:
            errorMessage("Unable to interpolate trajectory!", 2, info='Please define a trajectory source in the "sources" tab, save and then restart BAM', title='Yikes!', detail='{:}'.format(e))
            return None

        for pt in points:
            source = Supracenter(Position(pt[0], pt[1], pt[2]), pt[3])

            S = Source('Traj point - {:.2f} km'.format(pt[2]/1000), 'Fragmentation', source)

            self.bam.source_list.append(S)

        save(self.obj, True)
        loadSourcesIntoBam(self.obj.bam)

        self.close()