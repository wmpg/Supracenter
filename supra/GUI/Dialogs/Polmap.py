import os
import numpy as np
import pickle
import time

from functools import partial

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *

from supra.Files.SaveObjs import Prefs
from supra.Utils.AngleConv import roundToNearest

from supra.Utils.Classes import Position, RectangleItem

class Polmap(QWidget):

    def __init__(self, bam, points):

        QWidget.__init__(self)
        
        self.bam = bam
        self.points = points

        self.buildGUI()
        self.createGrid()

    def buildGUI(self):

        self.setWindowTitle('Polmap')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)
     
        self.polmap_view = pg.GraphicsLayoutWidget()
        self.polmap_canvas = self.polmap_view.addPlot()

        self.polmap_height = QSlider(Qt.Horizontal)
        self.polmap_grid = QSlider(Qt.Horizontal)

        self.polmap_height.setMinimum(0)
        self.polmap_height.setMaximum(500)
        self.polmap_height.setValue(340)
        self.polmap_height.setSingleStep(1)
        self.polmap_height.sliderMoved.connect(self.heightChange)
        self.polmap_height.valueChanged.connect(self.heightChange)

        self.polmap_grid.setMinimum(1)
        self.polmap_grid.setMaximum(100)
        self.polmap_grid.setValue(10)
        self.polmap_grid.setSingleStep(1)
        self.polmap_grid.sliderMoved.connect(self.gridChange)
        self.polmap_grid.valueChanged.connect(self.gridChange)

        self.polmap_height_l = QLabel("34000")
        self.polmap_grid_l = QLabel("0.01")

        layout.addWidget(self.polmap_view, 0, 0)
        layout.addWidget(self.polmap_height, 1, 0)
        layout.addWidget(self.polmap_grid, 2, 0)
        layout.addWidget(self.polmap_height_l, 1, 1)
        layout.addWidget(self.polmap_grid_l, 2, 1)

        self.station_marker = [None]*len(self.bam.stn_list)

        for ii, stn in enumerate(self.bam.stn_list):
            self.station_marker[ii] = pg.ScatterPlotItem()

            if stn.polarization.azimuth is not None:
                color = (255,   0,   0)
            else:
                color = (255, 255, 255)

            self.station_marker[ii].setPoints(x=[stn.metadata.position.lon], y=[stn.metadata.position.lat], pen=color, brush=color, symbol='t')
            self.polmap_canvas.addItem(self.station_marker[ii], update=True)
            txt = pg.TextItem("{:}".format(stn.metadata.code))
            txt.setPos(stn.metadata.position.lon, stn.metadata.position.lat)
            self.polmap_canvas.addItem(txt)

    def heightChange(self):
        self.polmap_height_l.setText(str(self.polmap_height.value()*100))
        self.drawGrid()

    def gridChange(self):
        self.polmap_grid_l.setText(str(self.polmap_grid.value()/1000))
        self.createGrid()

    def createGrid(self):

        tol = 500

        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

        radius = self.bam.setup.deg_radius

        min_lat = ref_pos.lat - radius
        max_lat = ref_pos.lat + radius
        min_lon = ref_pos.lon - radius
        max_lon = ref_pos.lon + radius
        min_elev = 0
        max_elev = 50000

        # height = float(self.polmap_height_l.text())

        self.grid = float(self.polmap_grid_l.text())

        self.grid_lat = np.arange(min_lat, max_lat, self.grid)
        self.grid_lon = np.arange(min_lon, max_lon, self.grid)
        self.grid_elev = np.arange(min_elev, max_elev, tol)

        self.grid_array = np.zeros((len(self.grid_lat), len(self.grid_lon), len(self.grid_elev)))

        pts_at_height = []

        for pt in self.points:

            # if roundToNearest(pt.elev, 5*tol) == roundToNearest(height, 5*tol):
            elev_idx = np.nanargmin(np.abs(self.grid_elev - pt.elev))
            lat_idx = np.nanargmin(np.abs(self.grid_lat - pt.lat))
            lon_idx = np.nanargmin(np.abs(self.grid_lon - pt.lon))

            self.grid_array[lat_idx, lon_idx, elev_idx] += 1


        self.drawGrid()
        
        


    def drawGrid(self):
        
        try:
            self.clearGrid()
        except:
            pass

        height = float(self.polmap_height_l.text())
        
        data = []

        h_idx = np.nanargmin(np.abs(height - self.grid_elev))

        indicies = np.argwhere(self.grid_array[:, :, h_idx]>0.0)

        for ii in indicies:
            data.append((self.grid_lon[ii[1]], self.grid_lat[ii[0]], self.grid, \
                            self.grid, self.grid_array[ii[0], ii[1], h_idx]))


        self.heat_map_data_squares = RectangleItem(data, c_map="reverse")
        self.polmap_canvas.addItem(self.heat_map_data_squares)

        self.drawTrajPoint(height)

    def clearGrid(self):
        self.polmap_canvas.removeItem(self.heat_map_data_squares)
        self.polmap_canvas.removeItem(self.traj_marker)

    def drawTrajPoint(self, height):
    
        traj_pos = self.bam.setup.trajectory.findGeo(height)

        self.traj_marker = pg.ScatterPlotItem()
        self.traj_marker.setPoints(x=[traj_pos.lon], y=[traj_pos.lat], pen=(255, 255, 255), brush=(255, 255, 255), symbol='+')
        self.polmap_canvas.addItem(self.traj_marker)