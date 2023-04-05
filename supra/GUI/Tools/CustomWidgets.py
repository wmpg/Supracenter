import os
import time

from functools import partial

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

import pyqtgraph as pg

from supra.GUI.Tools.GUITools import createLabelEditObj, createButton

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import matplotlib.pyplot as plt



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

class TauEx(QWidget):
    def __init__(self, energy, heights):
        super(TauEx, self).__init__()

        self.mark_for_deletion = False
        self.energy = energy

        self.heights = heights

        main_layout = QGridLayout()

        self.toggle = ToggleButton(True, 0)
        main_layout.addWidget(self.toggle, 1, 1)
        self.toggle.clicked.connect(self.toggle.clickedEvt)
        
        self.icon = QLabel()

        self.energy_label = QLabel("{:} @ {:} km".format(self.energy.source_type, self.energy.height/1000))
        main_layout.addWidget(self.energy_label, 1, 2)

        self.energy_amount_label = QLabel("{:.2E} J/{:.2f} kT TNT".format(self.energy.chem_pres, self.energy.chem_pres/4.184e12))
        main_layout.addWidget(self.energy_amount_label, 2, 2)

        try:
            station_name = "{:}-{:}".format(self.energy.station.metadata.network, self.energy.station.metadata.code)
        except AttributeError:
            station_name = "Unknown Station"

        self.energy_station_label = QLabel("Station: {:}".format(station_name))
        main_layout.addWidget(self.energy_station_label, 3, 2)


        # tau_label, self.tau_edits = createLabelEditObj("Tau: ", main_layout, 1, width=1, h_shift=2, tool_tip='', validate='float', default_txt='{:.2f}'.format(tau))
        self.delete_self = createButton("Delete", main_layout, 2, 5, self.deleteSelf, args=[])

        FRAG_ICON = os.path.join('supra', 'GUI', 'Images', "Fragmentation.png")
        TRAJ_ICON = os.path.join('supra', 'GUI', 'Images', "Trajectory.png")

        if heights[0] is None:
            min_h_label, self.min_h_edits = createLabelEditObj("height_min [km]: ", main_layout, 3, width=1, h_shift=2, tool_tip='', validate='float', default_txt='{:.2f}'.format(self.energy.height/1000))
            max_h_label, self.max_h_edits = createLabelEditObj("height_max [km]: ", main_layout, 4, width=1, h_shift=2, tool_tip='', validate='float', default_txt='{:.2f}'.format(self.energy.height/1000))
        else:
            min_h_label, self.min_h_edits = createLabelEditObj("height_min [km]: ", main_layout, 3, width=1, h_shift=2, tool_tip='', validate='float', default_txt='{:.2f}'.format(heights[0]))
            max_h_label, self.max_h_edits = createLabelEditObj("height_max [km]: ", main_layout, 4, width=1, h_shift=2, tool_tip='', validate='float', default_txt='{:.2f}'.format(heights[1]))
        
        if self.energy.source_type.lower() == "fragmentation": 
            pixmap = QPixmap(FRAG_ICON)
        elif self.energy.source_type.lower() == "ballistic":
            pixmap = QPixmap(TRAJ_ICON)
        else:
            print("Unknown source type")
            pixmap = QPixmap()

        self.icon.setPixmap(pixmap)
        main_layout.addWidget(self.icon, 1, 0)

        self.setMinimumHeight(126)

        self.setLayout(main_layout)

        self.setStyleSheet("""
                QGroupBox {
                    border: 3px solid grey;
                    border-radius: 10px;
                    }

                QPushButton{color: white; background-color: rgb(0, 100, 200);}
                """)
    # def getTau(self):
    #     return float(self.tau_edits.text())

    def getHeights(self):
        min_height = float(self.min_h_edits.text())
        max_height = float(self.max_h_edits.text())

        return [min_height, max_height]

    def deleteSelf(self):
        self.mark_for_deletion = True


class ToggleButton(QAbstractButton):

    """ Creates a basic toggle button with a picture as one of the following.
    A more general version.

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
        elif dict_item == 8:
            on_pix = 'save_icon_down.png'
            off_pix = 'save_icon_up.png'
        elif dict_item == 9:
            on_pix = 'Invert_down.png'
            off_pix = 'Invert_up.png'
        elif dict_item == 10:
            on_pix = 'show_frags_down.png'
            off_pix = 'show_frags_up.png'
        elif dict_item == 11:
            on_pix = 'show_traj_down.png'
            off_pix = 'show_traj_up.png'
        elif dict_item == 12:
            on_pix = 'perts_down.png'
            off_pix = 'perts_up.png'
        elif dict_item == 13:
            on_pix = 'right_arrow_down.png'
            off_pix = 'right_arrow_up.png'
        elif dict_item == 14:
            on_pix = 'left_arrow_down.png'
            off_pix = 'left_arrow_up.png'
        elif dict_item == 15:
            on_pix = 'folder_down.png'
            off_pix = 'folder_up.png'
        elif dict_item == 16:
            on_pix = 'save_trace_down.png'
            off_pix = 'save_trace_up.png'
        elif dict_item == 17:
            on_pix = 'frag_contour_down.png'
            off_pix = 'frag_contour_up.png'   
        elif dict_item == 18:
            on_pix = 'ball_contour_down.png'
            off_pix = 'ball_contour_up.png'           
        elif dict_item == 19:
            on_pix = 'pick_save_down.png'
            off_pix = 'pick_save_up.png'
        elif dict_item == 20:
            on_pix = 'Export_times_down.png'
            off_pix = 'Export_times_up.png'   
        elif dict_item == 21:
            on_pix = 'calc_yield_down.png'
            off_pix = 'calc_yield_up.png'
        elif dict_item == 22:
            on_pix = 'tau_down.png'
            off_pix = 'tau_up.png'   
        elif dict_item == 23:
            on_pix = 'staff_down.png'
            off_pix = 'staff_up.png'   


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

    def switchState(self):
        self.setState(not self.getState())
        self.update()


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

        plt.style.use('dark_background')

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
        # self.ax = self.figure.add_subplot(111)

        # set the layout
        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def show(self):
        self.canvas.draw()

    # def clear(self):
    #     self.ax.clear()