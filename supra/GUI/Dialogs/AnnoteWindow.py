
import numpy as np

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.WidgetBuilder import center
from supra.Utils.Classes import Annote
from supra.Stations.ProcessStation import *


class AnnoteWindow(QWidget):

    def __init__(self, time, stn, bam, mode="new", an=None):

        QWidget.__init__(self)
        
        self.buildGUI(time)

        self.annote_color = QColor(0, 0, 255)

        self.stn = stn

        self.an = an

        self.bam = bam

        self.mode = mode
        
        if mode == "edit":
            trace = self.stn.stream[0]

            start_time = self.an.time
            end_time = self.an.time + self.an.length

            cut_waveform, cut_time = subTrace(trace, start_time, end_time, self.bam.setup.fireball_datetime)
            station_waveform = pg.PlotDataItem(x=cut_time, y=cut_waveform, pen='w')
            self.annote_waveform_canvas.addItem(station_waveform)

            self.loadAnnote(an)

    def color_picker(self):
        color = QColorDialog.getColor()
        
        self.color_button.setStyleSheet("background-color: %s" % color.name())
        self.annote_color = color

    def loadAnnote(self, an):
        self.title_edits.setText(an.title)
        self.time_edits.setText(str(an.time))
        self.length_edits.setText(str(an.length))
        # group = self.group_edits.currentText()
        # source = self.source_edits.currentText()
        # height = float(self.height_edits.text())
        # notes = self.notes_box.toPlainText()
        # color = self.annote_color

    def addAnnote(self):
        
        title = self.title_edits.text()
        time = float(self.time_edits.text())
        length = float(self.length_edits.text())
        group = self.group_edits.currentText()
        source = self.source_edits.currentText()
        height = float(self.height_edits.text())
        notes = self.notes_box.toPlainText()
        color = self.annote_color

        an = Annote(title, time, length, group, source, height, notes, color)

        if self.mode == 'new':
            self.stn.annotation.add(an)
        elif self.mode == 'edit':
            self.stn.annotation.overwrite(an)
        self.close()

    def delAnnote(self):
        self.stn.annotation.remove(self.an)
        self.close()

    def buildGUI(self, time):
        self.setWindowTitle('Edit Annotation')
        
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)
     
        title_label = QLabel("Title: ")
        layout.addWidget(title_label, 1, 1, 1, 1)

        self.title_edits = QLineEdit()
        layout.addWidget(self.title_edits, 1, 2, 1, 1)

        self.color_button = QPushButton()
        layout.addWidget(self.color_button, 1, 3, 1, 1)
        self.color_button.clicked.connect(self.color_picker)

        group_label = QLabel("Group: ")
        layout.addWidget(group_label, 2, 1, 1, 1)

        self.group_edits = QComboBox()
        layout.addWidget(self.group_edits, 2, 2, 1, 1)

        time_label = QLabel("Time: ")
        layout.addWidget(time_label, 3, 1, 1, 1)

        self.time_edits = QLineEdit("{:.2f}".format(time))
        layout.addWidget(self.time_edits, 3, 2, 1, 1)
        self.time_edits.setValidator(QDoubleValidator())

        length_label = QLabel("Length: ")
        layout.addWidget(length_label, 4, 1, 1, 1)

        self.length_edits = QLineEdit('0')
        layout.addWidget(self.length_edits, 4, 2, 1, 1)
        self.length_edits.setValidator(QDoubleValidator())

        source_label = QLabel("Source: ")
        layout.addWidget(source_label, 5, 1, 1, 1)

        self.source_edits = QComboBox()
        layout.addWidget(self.source_edits, 5, 2, 1, 1)

        height_label = QLabel("Height: ")
        layout.addWidget(height_label, 6, 1, 1, 1)

        self.height_edits = QLineEdit('0')
        layout.addWidget(self.height_edits, 6, 2, 1, 1)
        self.height_edits.setValidator(QDoubleValidator())

        notes_label = QLabel("Notes: ")
        layout.addWidget(notes_label, 7, 1, 1, 1)

        self.notes_box = QPlainTextEdit()
        layout.addWidget(self.notes_box, 7, 2, 1, 2)

        self.annote_waveform_view = pg.GraphicsLayoutWidget()
        self.annote_waveform_canvas = self.annote_waveform_view.addPlot()
        layout.addWidget(self.annote_waveform_view, 8, 1, 1, 2)

        save_button = QPushButton('Save')
        layout.addWidget(save_button, 10, 2, 1, 1)
        save_button.clicked.connect(self.addAnnote)

        del_button = QPushButton('Delete')
        layout.addWidget(del_button, 10, 1, 1, 1)
        del_button.clicked.connect(self.delAnnote)