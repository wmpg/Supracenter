from functools import partial

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

class StationEx(QGroupBox):

    def __init__(self):
        super(StationEx, self).__init__()

        # self.setTitle('')
        main_layout = QHBoxLayout()

        left_side = QGridLayout()
        layout = QGridLayout()

        main_layout.addLayout(left_side)
        main_layout.addLayout(layout)

        self.toggle = ToggleButton(True)
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
    def __init__(self, status):
        super(ToggleButton, self).__init__()

        self.status = status
        self.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.update()

    def setState(self, state):
        self.status = state

        self.update()

    def getState(self):
        return self.status

    def paintEvent(self, event):
        painter = QPainter(self)

        if self.status:
            painter.drawPixmap(event.rect(), QPixmap('supra\\GUI\\Images\\ButtonStateOn.png'))
        else:
            painter.drawPixmap(event.rect(), QPixmap('supra\\GUI\\Images\\ButtonStateOff.png'))

    def sizeHint(self):
        return QPixmap('supra\\GUI\\Images\\ButtonStateOff.png').size()

    def clickedEvt(self):
        self.status = not self.status
        self.update()