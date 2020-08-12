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
        self.polmap_grid.setValue(25)
        self.polmap_grid.setSingleStep(1)
        self.polmap_grid.sliderMoved.connect(self.gridChange)
        self.polmap_grid.valueChanged.connect(self.gridChange)

        self.polmap_height_l = QLabel("34000")
        self.polmap_grid_l = QLabel("0.25")

        layout.addWidget(self.polmap_view, 0, 0)
        layout.addWidget(self.polmap_height, 1, 0)
        layout.addWidget(self.polmap_grid, 2, 0)
        layout.addWidget(self.polmap_height_l, 1, 1)
        layout.addWidget(self.polmap_grid_l, 2, 1)

        self.station_marker = [None]*len(self.bam.stn_list)

        for ii, stn in enumerate(self.bam.stn_list):
            self.station_marker[ii] = pg.ScatterPlotItem()
            self.station_marker[ii].setPoints(x=[stn.metadata.position.lon], y=[stn.metadata.position.lat], pen=(255, 255, 255), brush=(255, 255, 255), symbol='t')
            self.polmap_canvas.addItem(self.station_marker[ii], update=True)
            txt = pg.TextItem("{:}".format(stn.metadata.code))
            txt.setPos(stn.metadata.position.lon, stn.metadata.position.lat)
            self.polmap_canvas.addItem(txt)

    def heightChange(self):
        self.polmap_height_l.setText(str(self.polmap_height.value()*100))
        self.createGrid()

    def gridChange(self):
        self.polmap_grid_l.setText(str(self.polmap_grid.value()/100))
        self.createGrid()

    def createGrid(self):

        try:
            self.clearGrid()
        except:
            pass

        tol = 1000

        ref_pos = Position(self.bam.setup.lat_centre, self.bam.setup.lon_centre, 0)

        radius = self.bam.setup.deg_radius

        min_lat = ref_pos.lat - radius
        max_lat = ref_pos.lat + radius
        min_lon = ref_pos.lon - radius
        max_lon = ref_pos.lon + radius

        height = float(self.polmap_height_l.text())

        grid = float(self.polmap_grid_l.text())

        grid_lat = np.arange(min_lat, max_lat, grid)
        grid_lon = np.arange(min_lon, max_lon, grid)

        grid_array = np.zeros((len(grid_lat), len(grid_lon)))

        pts_at_height = []

        for pt in self.points:
            if roundToNearest(pt.elev, tol) == roundToNearest(height, tol):
                result_lat = roundToNearest(pt.lat, grid)
                result_lon = roundToNearest(pt.lon, grid)

                lat_idx = np.nanargmin(np.abs(grid_lat - result_lat))
                lon_idx = np.nanargmin(np.abs(grid_lon - result_lon))

                grid_array[lat_idx, lon_idx] += 1

        self.drawGrid(grid_array, grid_lat, grid_lon, grid)
        self.drawTrajPoint(height)

    def drawGrid(self, grid_array, grid_lat, grid_lon, grid):
        
        data = []

        for ii in range(len(grid_lat)):
            for jj in range(len(grid_lon)):
                if grid_array[ii][jj] > 0:
                    data.append((grid_lon[jj], grid_lat[ii], grid, grid, grid_array[ii][jj]))

        self.heat_map_data_squares = RectangleItem(data, c_map="reverse")
        self.polmap_canvas.addItem(self.heat_map_data_squares)

    def clearGrid(self):
        self.polmap_canvas.removeItem(self.heat_map_data_squares)
        self.polmap_canvas.removeItem(self.traj_marker)

    def drawTrajPoint(self, height):
    
        traj_pos = self.bam.setup.trajectory.findGeo(height)

        self.traj_marker = pg.ScatterPlotItem()
        self.traj_marker.setPoints(x=[traj_pos.lon], y=[traj_pos.lat], pen=(255, 255, 255), brush=(255, 255, 255), symbol='+')
        self.polmap_canvas.addItem(self.traj_marker)