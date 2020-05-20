
import numpy as np

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.WidgetBuilder import center
from supra.Utils.Classes import Annote


class AnnoteWindow(QWidget):

    def __init__(self, time, stn):

        QWidget.__init__(self)
        
        self.buildGUI(time)

        self.annote_color = QColor(0, 0, 255)

        self.stn = stn
        

    def color_picker(self):
        color = QColorDialog.getColor()
        
        self.color_button.setStyleSheet("background-color: %s" % color.name())
        self.annote_color = color

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
        self.stn.annotations.append(an)
        self.close()

    def delAnnote(self):
        pass

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

        save_button = QPushButton('Save')
        layout.addWidget(save_button, 8, 2, 1, 1)
        save_button.clicked.connect(self.addAnnote)

        del_button = QPushButton('Delete')
        layout.addWidget(del_button, 8, 1, 1, 1)
        del_button.clicked.connect(self.delAnnote)