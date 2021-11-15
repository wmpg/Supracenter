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

from supra.Utils.Classes import Position, Supracenter, Trajectory, Angle



class Source:

    def __init__(self, title, source_type, source, notes='', color=(255, 255, 255)):

        self.title = title
        self.color = color
        self.notes = notes
        self.source_type = source_type
        self.source = source

    def __str__(self):
        return "Source Obj: {:} type with name {:}".format(self.source_type, self.title)

class SourceWindow(QWidget):

    def __init__(self, bam):

        QWidget.__init__(self)
        
        self.bam = bam
        self.annote_color = (255, 255, 255)

        self.buildGUI()

    def color_picker(self):
        color = QColorDialog.getColor()
        
        self.color_button.setStyleSheet("background-color: %s" % color.name())
        self.annote_color = color

    def buildGUI(self):
        self.setWindowTitle('Add Source')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'Trajectory.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)
        
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        self.layout = QGridLayout()
        self.setLayout(self.layout)

        _, self.source_type = createComboBoxObj("Source Type", self.layout, 3, items=["Fragmentation", "Ballistic"], width=1, h_shift=0, tool_tip='')
        self.source_type.currentIndexChanged.connect(self.sourceChanged)

        self.clearLayout()
        self.sourceChanged()

    def clearLayout(self):

        for i in reversed(range(self.layout.count())):
            self.layout.itemAt(i).widget().setParent(None)


        title_label = QLabel("Title: ")
        self.layout.addWidget(title_label, 1, 1, 1, 1)

        self.title_edits = QLineEdit()
        self.layout.addWidget(self.title_edits, 1, 2, 1, 1)

        self.color_button = QPushButton()
        self.layout.addWidget(self.color_button, 1, 3, 1, 1)
        self.color_button.clicked.connect(self.color_picker)

        notes_label = QLabel("Description")
        self.layout.addWidget(notes_label, 2, 1)

        self.notes_box = QPlainTextEdit()
        self.layout.addWidget(self.notes_box, 2, 2, 1, 2)

        self.add_source_button = createButton("Add Source", self.layout, 14, 2, self.addSources)


    def sourceChanged(self):

        source_type = self.source_type.currentText()

        self.clearLayout()

        if source_type == "Fragmentation":

            _, self.lat_edit = createLabelEditObj('Latitude: ', self.layout, 4)
            _, self.lon_edit = createLabelEditObj('Longitude: ', self.layout, 5)
            _, self.elev_edit = createLabelEditObj('Elevation [m]: ', self.layout, 6)
            _, self.time_edit = createLabelEditObj('Time after reference [s]: ', self.layout, 7)

        elif source_type == "Ballistic":

            _, self.lat_i_edit = createLabelEditObj('Latitude (initial): ', self.layout, 4)
            _, self.lat_f_edit = createLabelEditObj('Latitude (final): ', self.layout, 5)
            _, self.lon_i_edit = createLabelEditObj('Longitude (initial): ', self.layout, 6)
            _, self.lon_f_edit = createLabelEditObj('Longitude (final): ', self.layout, 7)
            _, self.elev_i_edit = createLabelEditObj('Elevation (initial) [m]: ', self.layout, 8)
            _, self.elev_f_edit = createLabelEditObj('Elevation (final) [m]: ', self.layout, 9)
            _, self.azimuth_edit = createLabelEditObj('Azimuth: ', self.layout, 10)
            _, self.zenith_edit = createLabelEditObj('Zenith: ', self.layout, 11)
            _, self.time_edit = createLabelEditObj('Time after reference [s]: ', self.layout, 12)
            _, self.velocity_edit = createLabelEditObj('Velocity [m/s]: ', self.layout, 13)

        else:
            print("Unknown source")

        _, self.source_type = createComboBoxObj("Source Type", self.layout, 3, items=["Fragmentation", "Ballistic"], width=1, h_shift=0, tool_tip='')
        self.source_type.setCurrentText(source_type)
        self.source_type.currentIndexChanged.connect(self.sourceChanged)

    def addSources(self):
        
        title = self.title_edits.text()
        color = self.annote_color
        notes = self.notes_box.toPlainText()
        source_type = self.source_type.currentText()

        if source_type == "Fragmentation":
            source = Supracenter(Position(float(self.lat_edit.text()),\
                                                float(self.lon_edit.text()), \
                                                float(self.elev_edit.text())), \
                                                float(self.time_edit.text()))

        elif source_type == "Ballistic":

            try:
                initial_pos = Position(float(self.lat_i_edit.text()), float(self.lon_i_edit.text()), float(self.elev_i_edit.text()))
            except ValueError:
                initial_pos = Position(None, None, None)
            try:    
                final_pos =   Position(float(self.lat_f_edit.text()), float(self.lon_f_edit.text()), float(self.elev_f_edit.text()))
            except ValueError:
                final_pos = Position(None, None, None)

            t = float(self.time_edit.text())
            v = float(self.velocity_edit.text())
            
            try:
                az = Angle(float(self.azimuth_edit.text()))
            except ValueError:
                az = None

            try:    
                ze = Angle(float(self.zenith_edit.text()))
            except ValueError:
                ze = None

            source = Trajectory(t, v, zenith=ze, azimuth=az, pos_i=initial_pos, pos_f=final_pos, v_f=None)

            print(source)
        else:
            print("Error - unknown source")

        S = Source(title, source_type, source, notes=notes, color=color)
        

        if not hasattr(self.bam, "source_list"):
            self.bam.source_list = []

        self.bam.source_list.append(S)

        self.close()
if __name__ == '__main__':

    pass