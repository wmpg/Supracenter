import os

from functools import partial

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

import pyqtgraph as pg

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure



class SourceEx(QGroupBox):

    def __init__(self, S):
        super(SourceEx, self).__init__()

        self.source = S


        # self.setTitle('')
        main_layout = QGridLayout()

        left_side = QGridLayout()
        self.layout = QGridLayout()

        main_layout.addLayout(left_side, 0, 0, 2, 1)
        main_layout.addLayout(self.layout, 0, 1, 2, 100)

        self.toggle = ToggleButton(True, 0)
        left_side.addWidget(self.toggle, 1, 1)
        self.toggle.clicked.connect(self.toggle.clickedEvt)
        
        self.icon = QLabel()

        FRAG_ICON = os.path.join('supra', 'GUI', 'Images', "Fragmentation.png")
        TRAJ_ICON = os.path.join('supra', 'GUI', 'Images', "Trajectory.png")

        if self.source.source_type == "Fragmentation": 
            pixmap = QPixmap(FRAG_ICON)
            self.buildFragSource(self.layout)
        elif self.source.source_type == "Ballistic":
            pixmap = QPixmap(TRAJ_ICON)
            self.buildTrajSource(self.layout)
        else:
            print("Unknown source type")
            pixmap = QPixmap()

        self.icon.setPixmap(pixmap)
        left_side.addWidget(self.icon, 1, 0)

        self.setMinimumHeight(126)

        self.setLayout(main_layout)

        self.setStyleSheet("""
                QGroupBox {
                    border: 3px solid grey;
                    border-radius: 10px;
                    }
                """)


    def buildFragSource(self, layout):

        self.srctype_label = QLabel(self.source.title)
        # self.srctype_label.setStyleSheet("background-color: %s" % self.source.color.name())
        layout.addWidget(self.srctype_label, 1, 1)

        self.notes_label = QLabel(self.source.notes)
        layout.addWidget(self.notes_label, 1, 2)

        lat = self.source.source.position.lat
        lon = self.source.source.position.lon
        elev = self.source.source.position.elev
        time = self.source.source.time

        self.lat_label = QLabel("Latitude {:.4f}°N".format(lat))
        layout.addWidget(self.lat_label, 2, 1)

        self.lon_label = QLabel("Longitude {:.4f}°E".format(lon))
        layout.addWidget(self.lon_label, 2, 2)

        self.elev_label = QLabel("Elevation {:.2f} km".format(elev/1000))
        layout.addWidget(self.elev_label, 2, 3)

        self.time_label = QLabel("Time {:.2f} s".format(time))
        layout.addWidget(self.time_label, 2, 4)

    def buildTrajSource(self, layout):
        self.srctype_label = QLabel(self.source.title)
        # self.srctype_label.setStyleSheet("background-color: %s" % self.source.color.name())
        layout.addWidget(self.srctype_label, 1, 1)

        self.notes_label = QLabel(self.source.notes)
        layout.addWidget(self.notes_label, 1, 2)

        t = self.source.source.t
        v = self.source.source.v
        az = self.source.source.azimuth.deg
        ze = self.source.source.zenith.deg
        lat = self.source.source.pos_i.lat
        lon = self.source.source.pos_i.lon
        elev = self.source.source.pos_i.elev

        self.lat_label = QLabel("Latitude {:.4f}°N".format(lat))
        layout.addWidget(self.lat_label, 2, 1)

        self.lon_label = QLabel("Longitude {:.4f}°E".format(lon))
        layout.addWidget(self.lon_label, 2, 2)

        self.elev_label = QLabel("Elevation {:.2f} km".format(elev/1000))
        layout.addWidget(self.elev_label, 2, 3)

        self.time_label = QLabel("Time {:.2f} s".format(t))
        layout.addWidget(self.time_label, 2, 4)

        self.v_label = QLabel("Speed {:.2f} km/s".format(v/1000))
        layout.addWidget(self.v_label, 3, 1)

        self.az_label = QLabel("Azimuth {:.2f}°".format(az))
        layout.addWidget(self.az_label, 3, 2)

        self.ze_label = QLabel("Zenith {:.2f}°".format(ze))
        layout.addWidget(self.ze_label, 3, 3)



class StationEx(QGroupBox):

    def __init__(self, infrasound=False):
        super(StationEx, self).__init__()

        self.infrasound = infrasound

        # self.setTitle('')
        main_layout = QHBoxLayout()

        left_side = QGridLayout()
        layout = QGridLayout()

        main_layout.addLayout(left_side)
        main_layout.addLayout(layout)

        self.toggle = ToggleButton(True, 0)
        left_side.addWidget(self.toggle, 1, 1)
        self.toggle.clicked.connect(self.toggle.clickedEvt)
        self.toggle.clicked.connect(partial(self.func, self.toggle))

        self.network = QLabel()
        layout.addWidget(self.network, 1, 1)

        self.code = QLabel()
        layout.addWidget(self.code, 1, 2)

        self.name = QLabel()
        layout.addWidget(self.name, 1, 3)

        self.position = PositionEx()
        layout.addWidget(self.position, 2, 1, 1, 4)

        self.stream = []
        self.response = []

        self.setMinimumHeight(126)

        if self.toggle.getState():
            self.setStyle('enabled')
        else:
            self.setStyle('disabled')

        self.setLayout(main_layout)


    def setStyle(self, state):
        if state == 'enabled':
            if self.infrasound:
                self.setStyleSheet("""
                    QGroupBox {
                        border: 3px solid red;
                        border-radius: 10px;
                        background-color: rgb(10, 10, 10);
                        }
                    """)
            else:
                self.setStyleSheet("""
                    QGroupBox {
                        border: 3px solid rgb(0, 100, 200);
                        border-radius: 10px;
                        background-color: rgb(10, 10, 10);
                        }
                    """)
            self.network.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.code.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.name.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.name.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.position.poslabel.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.position.latlabel.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.position.lat.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.position.lonlabel.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.position.lon.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.position.elevlabel.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")
            self.position.elev.setStyleSheet("color: white; background-color: rgb(10, 10, 10)")

        elif state == 'disabled':
            self.setStyleSheet("""
                QGroupBox {
                    border: 3px solid grey;
                    border-radius: 10px;
                    background-color: rgb(125, 125, 125);
                    }
                """)
            self.network.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.code.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.name.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.name.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.position.poslabel.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.position.latlabel.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.position.lat.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.position.lonlabel.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.position.lon.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.position.elevlabel.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")
            self.position.elev.setStyleSheet("color: black; background-color: rgb(125, 125, 125)")

    def func(self, toggle):
        
        if toggle.getState():
            self.setStyle('enabled')
        else:
            self.setStyle('disabled')

class PositionEx(QWidget):

    def __init__(self):
        super(PositionEx, self).__init__()

        layout = QHBoxLayout()

        self.poslabel = QLabel('Position:')
        self.latlabel = QLabel('Lat: ')
        self.lat = QLabel()
        self.lonlabel = QLabel('Lon: ')
        self.lon = QLabel()
        self.elevlabel = QLabel('Elev: ')
        self.elev = QLabel()

        layout.addWidget(self.poslabel)
        layout.addWidget(self.latlabel)
        layout.addWidget(self.lat)
        layout.addWidget(self.lonlabel)
        layout.addWidget(self.lon)
        layout.addWidget(self.elevlabel)
        layout.addWidget(self.elev)

        self.setLayout(layout)

class ToggleButton(QAbstractButton):

    """ Creates a basic toggle button with a picture as one of the following.
    A more general version

    """

    def __init__(self, status, dict_item):
        super(ToggleButton, self).__init__()

        # Preloads images here to make it faster
        on_pix, off_pix = self.detImg(dict_item)

        self.on_pix = os.path.join('supra', 'GUI', 'Images', on_pix)
        self.off_pix = os.path.join('supra', 'GUI', 'Images', off_pix)

        self.status = status
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.update()

    def detImg(self, dict_item):
        # Preloads images here to make it faster
        if dict_item == 0:
            on_pix = 'ButtonStateOn.png'
            off_pix = 'ButtonStateOff.png'
        elif dict_item == 1:
            on_pix = 'Pick_down.png'
            off_pix = 'Pick_up.png'
        elif dict_item == 2:
            on_pix = 'Annote_down.png'
            off_pix = 'Annote_up.png'
        elif dict_item == 3:
            on_pix = 'Ground_Motion_down.png'
            off_pix = 'Ground_Motion_up.png'
        elif dict_item == 4:
            on_pix = 'RM_Pick_down.png'
            off_pix = 'RM_Pick_up.png'
        elif dict_item == 5:
            on_pix = 'Pick_down_b.png'
            off_pix = 'Pick_up_b.png'
        elif dict_item == 6:
            on_pix = 'Bandpass_down.png'
            off_pix = 'Bandpass_up.png'
        elif dict_item == 7:
            on_pix = 'Polmap_down.png'
            off_pix = 'Polmap_up.png'


        return on_pix, off_pix

    def changeImg(self, dict_item):

        on_pix, off_pix = self.detImg(dict_item)
        self.on_pix = os.path.join('supra', 'GUI', 'Images', on_pix)
        self.off_pix = os.path.join('supra', 'GUI', 'Images', off_pix)
        self.update()

    def setState(self, state):
        self.status = state
        qApp.processEvents()
        self.update()

    def getState(self):
        return self.status

    def isChecked(self):
        return self.status

    def paintEvent(self, event):
        painter = QPainter(self)

        if self.status:
            painter.drawPixmap(event.rect(), QPixmap(self.on_pix))
        else:
            painter.drawPixmap(event.rect(), QPixmap(self.off_pix))

    def sizeHint(self):
        return QPixmap(self.on_pix).size()

    def clickedEvt(self):
        self.status = not self.status
        self.update()


class MatplotlibPyQT(QWidget):

    def __init__(self):
        super(MatplotlibPyQT, self).__init__() 


        # a figure instance to plot on
        self.figure = Figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        # this is the Navigation widget
        # it takes the Canvas widget and a parent
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setStyleSheet("background-color:Gray;")
        # create an axis
        self.ax = self.figure.add_subplot(111)
      
        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def show(self):
        self.canvas.draw()