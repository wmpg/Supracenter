import os

from functools import partial

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

import pyqtgraph as pg

class StationEx(QGroupBox):

    def __init__(self):
        super(StationEx, self).__init__()

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
            self.setStyleSheet("""
                QGroupBox {
                    border: 3px solid grey;
                    border-radius: 10px;
                    background-color: rgb(10, 10, 10);
                    }
                """)
        elif state == 'disabled':
            self.setStyleSheet("""
                QGroupBox {
                    border: 3px solid grey;
                    border-radius: 10px;
                    background-color: rgb(125, 125, 125);
                    }
                """)

    def func(self, toggle):
        
        if toggle.getState():
            self.setStyle('enabled')
        else:
            self.setStyle('disabled')

class PositionEx(QWidget):

    def __init__(self):
        super(PositionEx, self).__init__()

        layout = QHBoxLayout()

        pos = QLabel('Position:')
        lat = QLabel('Lat: ')
        self.lat = QLabel()
        lon = QLabel('Lon: ')
        self.lon = QLabel()
        elev = QLabel('Elev: ')
        self.elev = QLabel()

        layout.addWidget(pos)
        layout.addWidget(lat)
        layout.addWidget(self.lat)
        layout.addWidget(lon)
        layout.addWidget(self.lon)
        layout.addWidget(elev)
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

