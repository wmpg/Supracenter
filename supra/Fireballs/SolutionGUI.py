
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time
import datetime
import copy
import obspy
import scipy.signal

from netCDF4 import Dataset

from PyQt5.QtWidgets import *
import sys
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import matplotlib.gridspec as gridspec

from functools import partial

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from mpl_toolkits.mplot3d import Axes3D

from supra.Supracenter.stationDat import convStationDat
from supra.Supracenter.psoSearch import psoSearch
from supra.Fireballs.Program import position, configRead, configWrite, station
from supra.Fireballs.SeismicTrajectory import Constants, parseWeather, waveReleasePoint, timeOfArrival
from supra.Supracenter.angleConv import loc2Geo, geo2Loc, angle2NDE
from supra.Supracenter.convLevels import convLevels
from supra.Supracenter.netCDFconv import findECMWFSound
from supra.Supracenter.SPPT import perturb as perturbation_method
from supra.Supracenter.fetchCopernicus import copernicusAPI
from supra.Fireballs.MakeIRISPicks import WaveformPicker
from supra.Fireballs.GetIRISData import readStationAndWaveformsListFile, plotAllWaveforms, \
butterworthBandpassFilter, convolutionDifferenceFilter
from wmpl.Utils.Earth import greatCircleDistance
from wmpl.Utils.PlotMap import GroundMap
from supra.Supracenter.cyscan import cyscan
from supra.Supracenter.cyweatherInterp import getWeather
from supra.Supracenter.SPPT import perturb

global arrTimes 
global sounding

DATA_FILE = 'data.txt'
OUTPUT_CSV = 'data_picks.csv'

class SolutionGUI(QMainWindow):
    def __init__(self):
        super().__init__()

        self._main = QWidget()
        self.setCentralWidget(self._main)
        layout = QGridLayout(self._main)

        self.setWindowTitle('WMPG Seismic and Infrasound Meteor Program')
        app_icon = QtGui.QIcon()
        app_icon.addFile('wmpl.png', QtCore.QSize(16,16))
        self.setWindowIcon(app_icon)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        self.slider_scale = 0.25
        self.bandpass_scale = 0.1

        self.tab_widget = QTabWidget()
        self.tab_widget.blockSignals(True)
        self.tab_widget.currentChanged.connect(self.tabChange)
        # self.tab_widget.blockSignals(False)

        self.addIniWidgets()
        self.addSupraWidgets()
        self.addSupWidgets()
        self.addMakePicksWidgets()
        self.addFetchATMWidgets()
        self.addProfileWidgets()
        self.addDocsWidgets()

        self.var_typ = 't'

        self.tab_widget.blockSignals(False)
        layout.addWidget(self.tab_widget, 1, 1)

        menu_bar = self.menuBar() 
        layout.addWidget(menu_bar, 0, 1)
        file_menu = menu_bar.addMenu('&File')
        about_menu = menu_bar.addMenu('&About')

        file_exit = QAction("Exit", self)
        file_exit.triggered.connect(qApp.quit)
        file_menu.addAction(file_exit)

        stylesheet = """ 
        QTabWidget>QWidget>QWidget{background: gray;}
        QLabel{color: white;}
        QCheckBox(color: white;)
        """

        self.setStyleSheet(stylesheet)

        # timer = QtCore.QTimer()
        # timer.timeout.connect(self.update)
        # timer.start(0)
    

        pg.setConfigOptions(antialias=True)

    def toolTime(self, var):

        tool_time_dict = {}

        with open('supra/Fireballs/tool_tips.csv') as f:
            for line in f:
                line = line.split(':')
                tool_time_dict[line[0]] = line[1]

        return tool_time_dict[var]

    def errorMessage(self, message, level, info='', title='Yikes!', detail=''):
        """
        """
        msg = QMessageBox()

        if level == 0:
            msg.setIcon(QMessageBox.Information)
        elif level == 1:
            msg.setIcon(QMessageBox.Warning)
        else:
            msg.setIcon(QMessageBox.Critical)

        msg.setText(message)
        
        msg.setInformativeText(info)
        msg.setWindowTitle(title)
        msg.setDetailedText(detail)

        msg.setStandardButtons(QMessageBox.Ok)

        msg.exec_()

    def supSaveChanges(self):
        pass

    def supraSaveChanges(self):

        self.sounding_file_edits.setText(self.atmospheric_file_edit.text())
        self.station_picks_file_edits.setText(self.picks_file_edit.text())
        self.start_datetime_edits.setDateTime(self.ref_edit.dateTime().toPyDateTime())
        self.lat_frag_edits.setText(str(self.lat_edit.text()))
        self.lon_frag_edits.setText(str(self.lon_edit.text()))
        self.elev_frag_edits.setText(str(self.elev_edit.text()))
        self.time_frag_edits.setText(str(self.time_edit.text()))

        self.errorMessage('File Copied!', 0, title='Ini Information Copied')
   
    def addSupWidgets(self):

        sup_tab = QWidget()
        self.master_sup = QHBoxLayout()
        self.sup_tab_content = QGridLayout()
        self.sup_plots = QVBoxLayout()

        self.sup_two_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.sup_two_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.sup_plots.addWidget(self.sup_two_canvas)
        
        self.sup_three_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.sup_three_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.sup_plots.addWidget(self.sup_three_canvas)

        self.sup_sub_tab = QTabWidget()
        self.sup_sub_tab.blockSignals(True)
        self.sup_tab_content.addWidget(self.sup_sub_tab, 1, 1, 1, 2)

        bounds = QWidget()
        bounds_content = QGridLayout()
        bounds.setLayout(bounds_content)
        self.sup_sub_tab.addTab(bounds, "Bounds")

        self.sup_range_label = QLabel("Range: ")
        bounds_content.addWidget(self.sup_range_label, 1, 1)

        self.sup_north_label = QLabel("North: ")
        bounds_content.addWidget(self.sup_north_label, 1, 3)

        self.sup_north_edits = QLineEdit("")
        bounds_content.addWidget(self.sup_north_edits, 2, 3)
        
        self.sup_west_label = QLabel("West: ")
        bounds_content.addWidget(self.sup_west_label, 3, 1)

        self.sup_west_edits = QLineEdit("")
        bounds_content.addWidget(self.sup_west_edits, 3, 2)
        
        self.sup_east_label = QLabel("East: ")
        bounds_content.addWidget(self.sup_east_label, 3, 5)

        self.sup_east_edits = QLineEdit("")
        bounds_content.addWidget(self.sup_east_edits, 3, 4)
        
        self.sup_south_label = QLabel("South: ")
        bounds_content.addWidget(self.sup_south_label, 5, 3)

        self.sup_south_edits = QLineEdit("")
        bounds_content.addWidget(self.sup_south_edits, 4, 3) 

        self.sup_height_label = QLabel("Height: ")
        bounds_content.addWidget(self.sup_height_label, 2, 6)

        self.sup_height_max_edits = QLineEdit("")
        bounds_content.addWidget(self.sup_height_max_edits, 3, 6)

        self.sup_height_min_edits = QLineEdit("")
        bounds_content.addWidget(self.sup_height_min_edits, 4, 6)         

        self.sup_ref_time_label = QLabel("Reference Time: ")
        bounds_content.addWidget(self.sup_ref_time_label, 6, 1)

        self.sup_ref_time_edits = QDateTimeEdit()
        bounds_content.addWidget(self.sup_ref_time_edits, 6, 2, 1, 4)
        self.sup_ref_time_edits.setCalendarPopup(True)
        
        self.sup_time_label = QLabel("Time: ")
        bounds_content.addWidget(self.sup_time_label, 7, 1)
        
        self.sup_max_time_label = QLabel("Max: ")
        bounds_content.addWidget(self.sup_max_time_label, 7, 4)

        self.sup_max_time_edits = QLineEdit("")
        bounds_content.addWidget(self.sup_max_time_edits, 7, 5)
        
        self.sup_min_time_label = QLabel("Min: ")
        bounds_content.addWidget(self.sup_min_time_label, 7, 2)

        self.sup_min_time_edits = QLineEdit("")
        bounds_content.addWidget(self.sup_min_time_edits, 7, 3)

        self.sup_picks_search_label = QLabel("Picks File: ")
        bounds_content.addWidget(self.sup_picks_search_label, 8, 1, 1, 1)
        
        self.sup_picks_file_edit = QLineEdit("")
        bounds_content.addWidget(self.sup_picks_file_edit, 8, 2, 1, 3)

        self.sup_picks_file_button = QPushButton('Browse')
        bounds_content.addWidget(self.sup_picks_file_button, 8, 5, 1, 1)
        self.sup_picks_file_button.clicked.connect(partial(self.fileSearch, ["CSV Picks File (*.csv)"], self.sup_picks_file_edit))

        atm = QWidget()
        atm_content = QGridLayout()
        atm.setLayout(atm_content)
        self.sup_sub_tab.addTab(atm, "Atmosphere")

        self.sup_copy_button = QPushButton('Copy to INI')
        self.sup_tab_content.addWidget(self.sup_copy_button, 2, 1, 1, 1)
        self.sup_copy_button.clicked.connect(self.supSaveChanges)

        self.sup_search_button = QPushButton('Search')
        self.sup_tab_content.addWidget(self.sup_search_button, 2, 2, 1, 1)
        self.sup_search_button.clicked.connect(self.supSearch)

        self.sup_results_label = QLabel("Results: ")
        self.sup_tab_content.addWidget(self.sup_results_label, 3, 1, 1, 1)

        self.sup_results_table = QTableWidget(0, 0)
        self.sup_tab_content.addWidget(self.sup_results_table, 4, 1, 1, 2)

        self.master_sup.addLayout(self.sup_tab_content)
        self.master_sup.addLayout(self.sup_plots)

        sup_tab.setLayout(self.master_sup)
        self.tab_widget.addTab(sup_tab, "Supracenter PSO Search")
  
    def addSupraWidgets(self):

        #self.tab_widget.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        supra_tab = QWidget()
        self.master_supra = QHBoxLayout()
        self.supra_tab_content = QGridLayout()
        self.plots = QVBoxLayout()

        self.two_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.two_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.plots.addWidget(self.two_canvas)
        #self.two_ax = self.two_canvas.figure.subplots()
        
        # self.three_canvas = gl.GLViewWidget()
        # self.three_canvas.opts['distance'] = 20
        # g = gl.GLGridItem()
        # self.three_canvas.addItem(g)
        # self.three_ax = gl.GLScatterPlotItem(pxMode=False)
        # self.three_ax.translate(5,5,0)
        # self.three_canvas.addItem(self.three_ax)
        # self.three_canvas.show()
        # self.w = gl.GLViewWidget()
        # self.w.opts['distance'] = 40
        # self.w.setWindowTitle('graph')
        # self.w.setGeometry(0, 110, 1920, 1080)
        # self.w.show()
        
        # self.gx = gl.GLGridItem()
        # self.gx.scale(10000, 10000, 10000)
        # # gx.rotate(90, 0, 1, 0)
        # # gx.translate(-10, 0, 0)

        # pos = np.empty((53, 3))
        # size = np.empty((53))
        # color = np.empty((53, 4))
        # pos[0] = (1,0,0); size[0] = 0.5;   color[0] = (1.0, 0.0, 0.0, 0.5)
        # pos[1] = (0,1,0); size[1] = 0.2;   color[1] = (0.0, 0.0, 1.0, 0.5)
        # pos[2] = (0,0,1); size[2] = 2./3.; color[2] = (0.0, 1.0, 0.0, 0.5)

        # z = 0.5
        # d = 6.0
        # for i in range(3,53):
        #     pos[i] = (0,0,z)
        #     size[i] = 2./d
        #     color[i] = (0.0, 1.0, 0.0, 0.5)
        #     z *= 0.5
        #     d *= 2.0
            
        # self.p13d = gl.GLScatterPlotItem(pos=pos, size=size, color=color, pxMode=False)
        # self.p13d = gl.GLScatterPlotItem(pxMode=False)
        # self.w.addItem(self.p13d)
        # self.w.addItem(self.gx)


        # self.plots.addWidget(self.view)
        # self.plots.addWidget(self.three_canvas)

        self.three_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.three_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.plots.addWidget(self.three_canvas)
        #self.three_ax = self.three_canvas.figure.subplots()

        # self.two_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        # self.two_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        # self.plots.addWidget(self.two_canvas)
        #self.two_ax = self.two_canvas.figure.subplots()

        self.search_button = QPushButton('Search')
        self.supra_tab_content.addWidget(self.search_button, 11, 3, 1, 3)
        self.search_button.clicked.connect(self.supraSearch)

        self.lat_label = QLabel("Latitude: ")
        self.supra_tab_content.addWidget(self.lat_label, 1, 3, 1, 1)

        self.lon_label = QLabel("Longitude: ")
        self.supra_tab_content.addWidget(self.lon_label, 2, 3, 1, 1)

        self.elev_label = QLabel("Elevation: ")
        self.supra_tab_content.addWidget(self.elev_label, 3, 3, 1, 1)

        self.time_label = QLabel("Time: ")
        self.supra_tab_content.addWidget(self.time_label, 4, 3, 1, 1)

        self.ref_label = QLabel("Reference Datetime: ")
        self.supra_tab_content.addWidget(self.ref_label, 5, 3, 1, 1)

        self.picks_search_label = QLabel("Picks File: ")
        self.supra_tab_content.addWidget(self.picks_search_label, 6, 3, 1, 1)

        self.atmospheric_search_label = QLabel("Atmospheric File: ")
        self.supra_tab_content.addWidget(self.atmospheric_search_label, 7, 3, 1, 1)

        self.supra_save_changes_button = QPushButton('Copy to INI Builder')
        self.supra_tab_content.addWidget(self.supra_save_changes_button, 8, 3, 1, 1)
        self.supra_save_changes_button.clicked.connect(self.supraSaveChanges)

        self.results_label = QLabel("Results: ")
        self.supra_tab_content.addWidget(self.results_label, 9, 1, 1, 1)

        self.tableWidget = QTableWidget(0, 0)
        self.supra_tab_content.addWidget(self.tableWidget, 10, 1, 1, 10)

        self.lat_edit = QLineEdit("0")
        self.supra_tab_content.addWidget(self.lat_edit, 1, 4, 1, 2)

        self.lon_edit = QLineEdit("0")
        self.supra_tab_content.addWidget(self.lon_edit, 2, 4, 1, 2)

        self.elev_edit = QLineEdit("0")
        self.supra_tab_content.addWidget(self.elev_edit, 3, 4, 1, 2)

        self.time_edit = QLineEdit("0")
        self.supra_tab_content.addWidget(self.time_edit, 4, 4, 1, 2)

        self.ref_edit = QDateTimeEdit()
        self.supra_tab_content.addWidget(self.ref_edit, 5, 4, 1, 2)
        self.ref_edit.setCalendarPopup(True)

        self.picks_file_edit = QLineEdit("")
        self.supra_tab_content.addWidget(self.picks_file_edit, 6, 4, 1, 1)

        self.atmospheric_file_edit = QLineEdit("")
        self.supra_tab_content.addWidget(self.atmospheric_file_edit, 7, 4, 1, 1)

        self.picks_file_button = QPushButton('Browse')
        self.supra_tab_content.addWidget(self.picks_file_button, 6, 5, 1, 1)
        self.picks_file_button.clicked.connect(partial(self.fileSearch, ["CSV Picks File (*.csv)"], self.picks_file_edit))

        self.atmospheric_file_button = QPushButton('Browse')
        self.supra_tab_content.addWidget(self.atmospheric_file_button, 7, 5, 1, 1)
        self.atmospheric_file_button.clicked.connect(partial(self.fileSearch, ["NetCDF (*.nc)"], self.atmospheric_file_edit))
        self.master_supra.addLayout(self.supra_tab_content)
        self.master_supra.addLayout(self.plots)

        supra_tab.setLayout(self.master_supra)
        self.tab_widget.addTab(supra_tab, "Supracenter Manual Search")

    def fatmPlot(self):

        consts = Constants()
        setup = self.saveINI(False)

        setup.lat_centre = self.fatm_lat_slide.value()*self.slider_scale
        setup.lon_centre = self.fatm_lon_slide.value()*self.slider_scale
        setup.sounding_file = self.fatm_name_edits.text()
        setup.weather_type = 'ecmwf'

        try:
            dataset = parseWeather(setup, consts)
            sounding = findECMWFSound(setup.lat_centre, setup.lon_centre, dataset)
        except:
            self.errorMessage('Error reading weather profile in fatmPlotProfile', 2)
            return None

        self.atm_canvas.setLabel('left', "Height", units='m')
        if self.fatm_variable_combo.currentText() == 'Temperature':
            X = sounding[:, 1]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Temperature", units='K')
        elif self.fatm_variable_combo.currentText() == 'Wind Magnitude':
            X = sounding[:, 2]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Wind Magnitude", units='m/s')
        elif self.fatm_variable_combo.currentText() == 'Wind Direction':
            X = sounding[:, 3]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Wind Direction", units='deg from N')
        elif self.fatm_variable_combo.currentText() == 'U-Component of Wind':
            dirs = angle2NDE(np.degrees(sounding[:, 3]))
            mags = sounding[:, 2]
            X = mags*np.cos(np.radians(dirs))
            Y = sounding[:, 0]
        elif self.fatm_variable_combo.currentText() == 'V-Component of Wind':
            dirs = angle2NDE(np.degrees(sounding[:, 3]))
            mags = sounding[:, 2]
            X = mags*np.sin(np.radians(dirs))
            Y = sounding[:, 0]
        else:
            self.errorMessage('Error reading fatmPlotProfile combo box', 2, detail=self.fatm_variable_combo.currentText())
            X = []
            Y = []

        self.fatm_canvas.clear()
        self.fatm_canvas.plot(x=X, y=Y, pen='w')
        SolutionGUI.update(self)

    def fatmFetch(self, download):

        variables = []
        self.fatm_variable_combo.clear() 

        if self.fatm_temp.isChecked():
            variables.append('temperature')
            self.fatm_variable_combo.addItem('Temperature')
        if self.fatm_u_wind.isChecked():
            variables.append('u_component_of_wind')
            self.fatm_variable_combo.addItem('U-Component of Wind')
        if self.fatm_v_wind.isChecked():
            variables.append('v_component_of_wind')
            self.fatm_variable_combo.addItem('V-Component of Wind')
        if self.fatm_u_wind.isChecked() and self.fatm_v_wind.isChecked():
            self.fatm_variable_combo.addItem('Wind Magnitude')
            self.fatm_variable_combo.addItem('Wind Direction')

        year = str(self.fatm_datetime_edits.dateTime().date().year())
        month = str(self.fatm_datetime_edits.dateTime().date().month())
        day = str(self.fatm_datetime_edits.dateTime().date().day())

        time_of = str("{:02d}".format(self.fatm_datetime_edits.dateTime().time().hour())) + \
        ':' + str("{:02d}".format(self.fatm_datetime_edits.dateTime().time().minute()))

        loc = self.fatm_name_edits.text()

        if download:
            self.errorMessage('Please wait until your file has downloaded...', 0)
            copernicusAPI(variables, year, month, day, time_of, loc)
            return None
        else:
            self.fatmPlot()

    def fatmValueChange(self, obj, slider):
        
        if obj == self.fatm_lat_label:
            obj.setText('Latitude: {:8.2f}'.format(slider.value()*self.slider_scale))
        elif obj == self.fatm_lon_label:
            obj.setText('Longitude: {:8.2f}'.format(slider.value()*self.slider_scale))
        else:
            self.errorMessage('Bad atm slider pass in fatmValueChange', 2)

        #self.atmPlotProfile(self.atm_lat_slide.value(), self.atm_lon_slide.value(), self.var_typ)
    
    def parseGeneralECMWF(self, file_name, lat, lon, time_of, variables):

        try:
        # Read the file
            dataset = Dataset(file_name, "r+", format="NETCDF4")
        except:
            self.errorMessage("Unable to read weather file: {:}".format(file_name), 2)
            return None

        lon_index = int(np.around((lon%360)*4))
        lat_index = int(np.around(-(lat+90)*4))
        time_index = int(time_of)

        sounding = []
        pressures = np.array(dataset.variables['level'])
        # level = convLevels()
        # # level = np.flipud(np.array(level))

        # sounding.append(level)
        sounding.append(pressures)

        for var in variables:
            if (set([var]) < set(dataset.variables.keys())):
                sounding.append(np.array(dataset.variables[var][time_index, :, lat_index, lon_index]))
            else:
                print("WARNING: Variable {:} not found in dataset!".format(var))

        sounding = np.array(sounding).transpose()

        dataset.close()
        return sounding

    def fatmPrint(self):
        
        filename = QFileDialog.getSaveFileName(self, 'Save File', '', 'Text File (*.txt)')

        setup = self.saveINI(False)

        setup.lat_centre = self.fatm_lat_slide.value()*self.slider_scale
        setup.lon_centre = self.fatm_lon_slide.value()*self.slider_scale
        setup.sounding_file = self.fatm_name_edits.text()
        atm_time = self.fatm_datetime_edits.dateTime().time().hour()
        setup.weather_type = 'ecmwf'

        variables = []

        if self.fatm_temp.isChecked():
            variables.append('t')
        if self.fatm_u_wind.isChecked():
            variables.append('u')
        if self.fatm_v_wind.isChecked():
            variables.append('v')
            
        # try:
        sounding = self.parseGeneralECMWF(setup.sounding_file, setup.lat_centre, setup.lon_centre, atm_time, variables)
        #     # dataset = parseWeather(setup, consts)
        #     # sounding = findECMWFSound(setup.lat_centre, setup.lon_centre, dataset)
        # except:
        #     self.errorMessage('Error reading weather profile in fatmPrint', 2)
        #     return None

        if '.txt' not in filename[0]:
            filename[0] = filename[0] + '.txt'


        header = 'Pressures (hPa)'

        for element in variables:

            if element == variables[-1]:
                header = header + ', ' + element + '\n'
            else:
                header = header + ', ' + element



        with open(str(filename[0]), 'w') as f:
            
            f.write(header)

            for line in sounding:
                info = ''
                for element in line:
                    if element == line[-1]:
                        info = info + str(element) + '\n'
                    else:
                        info = info + str(element) + ','
                f.write(info)

        self.errorMessage('Printed out sounding data', 0, title="Print Done")

    def addFetchATMWidgets(self):
        fetch = QWidget()
        fetch_master = QHBoxLayout()
        fetch_content = QGridLayout()
        fetch_plots = QVBoxLayout()
        fetch_master.addLayout(fetch_plots)
        fetch_master.addLayout(fetch_content)
        fetch.setLayout(fetch_master)
        self.tab_widget.addTab(fetch, "Fetch Atmosphere")

        self.fatm_view = pg.GraphicsLayoutWidget()
        self.fatm_canvas = self.fatm_view.addPlot()
        fetch_plots.addWidget(self.fatm_view)
        self.fatm_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        self.fatm_variable_combo = QComboBox()
        fetch_plots.addWidget(self.fatm_variable_combo)
        self.fatm_variable_combo.currentTextChanged.connect(self.fatmPlot)

        self.fatm_name_label = QLabel("Name: ")
        fetch_content.addWidget(self.fatm_name_label, 1, 1)

        self.fatm_name_edits = QLineEdit("")
        fetch_content.addWidget(self.fatm_name_edits, 1, 2)

        self.fatm_button = QPushButton('Browse')
        fetch_content.addWidget(self.fatm_button, 1, 3)
        self.fatm_button.clicked.connect(partial(self.folderSearch, self.fatm_name_edits))

        self.fatm_datetime_label = QLabel("Time of Profile:")
        self.fatm_datetime_edits = QDateTimeEdit()
        fetch_content.addWidget(self.fatm_datetime_label, 2, 1)
        fetch_content.addWidget(self.fatm_datetime_edits, 2, 2)
        self.fatm_datetime_edits.setCalendarPopup(True)

        self.fatm_variable_group = QGroupBox("Variables")
        fetch_content.addWidget(self.fatm_variable_group, 4, 1, 1, 2)

        #############################
        group_box = QGridLayout()
        self.fatm_variable_group.setLayout(group_box)

        self.fatm_temp = QCheckBox('Temperature')
        group_box.addWidget(self.fatm_temp, 1, 1)

        self.fatm_u_wind = QCheckBox('U Wind')
        group_box.addWidget(self.fatm_u_wind, 1, 2)

        self.fatm_v_wind = QCheckBox('V Wind')
        group_box.addWidget(self.fatm_v_wind, 1, 3)
        ##############################

        self.fatm_fetch = QPushButton("Download")
        fetch_content.addWidget(self.fatm_fetch, 5, 1, 1, 2)
        self.fatm_fetch.clicked.connect(partial(self.fatmFetch, True))

        self.fatm_open = QPushButton("Open")
        fetch_content.addWidget(self.fatm_open, 6, 1, 1, 2)
        self.fatm_open.clicked.connect(partial(self.fatmFetch, False))

        self.fatm_print = QPushButton("Print")
        fetch_content.addWidget(self.fatm_print, 9, 1, 1, 2)
        self.fatm_print.clicked.connect(self.fatmPrint)

        self.fatm_lat_label = QLabel("Latitude: 0")
        self.fatm_lat_slide = QSlider(Qt.Horizontal)
        fetch_content.addWidget(self.fatm_lat_label, 7, 1)
        fetch_content.addWidget(self.fatm_lat_slide, 7, 2)
        self.fatm_lat_slide.setMinimum(-90/self.slider_scale)
        self.fatm_lat_slide.setMaximum(90/self.slider_scale)
        self.fatm_lat_slide.setValue(0)
        self.fatm_lat_slide.setTickInterval(0.5)
        self.fatm_lat_slide.valueChanged.connect(partial(self.fatmValueChange, self.fatm_lat_label, self.fatm_lat_slide))

        self.fatm_lon_label = QLabel("Longitude: 0")
        self.fatm_lon_slide = QSlider(Qt.Horizontal)
        fetch_content.addWidget(self.fatm_lon_label, 8, 1)
        fetch_content.addWidget(self.fatm_lon_slide, 8, 2,)
        self.fatm_lon_slide.setMinimum(-180/self.slider_scale)
        self.fatm_lon_slide.setMaximum(180/self.slider_scale)
        self.fatm_lon_slide.setValue(0)
        self.fatm_lon_slide.setTickInterval(0.5)
        self.fatm_lon_slide.valueChanged.connect(partial(self.fatmValueChange, self.fatm_lon_label, self.fatm_lon_slide))

    def addDocsWidgets(self):
        docs_tab = QWidget()
        docs_tab_content = QVBoxLayout()
        docs_tab.setLayout(docs_tab_content)

        self.tab_widget.addTab(docs_tab, "Documentation")

        self.docTree = QTreeWidget()

        root = QTreeWidgetItem(self.docTree, ["Documentation"])
        root.setData(2, Qt.EditRole, '')

        folder1 = QTreeWidgetItem(root, ["Usage"])
        folder1.setData(2, Qt.EditRole, '')

        folder1_1 = QTreeWidgetItem(folder1, ["INI Builder"])
        folder1_1.setData(2, Qt.EditRole, '') 

        folder1_1_1 = QTreeWidgetItem(folder1_1, ['''Used to create configuration files for this program. Configurations
can be loaded or saved through this tab. The other tabs in this program 
are linked to these ones, and by default, will run off of these 
configurations, unless specified on the respective tab.'''])
        folder1_1_1.setData(2, Qt.EditRole, '') 

        folder1_2 = QTreeWidgetItem(folder1, ["Supracenter Manual Search"])
        folder1_2.setData(2, Qt.EditRole, '') 

        folder1_2_1 = QTreeWidgetItem(folder1_2, ['''Used to search for an individual Supracenter point, and return the 
residuals from each station. The Supracenter can be adjusted using the edit boxes on the tab.'''])
        folder1_2_1.setData(2, Qt.EditRole, '') 


        folder1_3 = QTreeWidgetItem(folder1, ["Atmospheric Profile"])
        folder1_3.setData(2, Qt.EditRole, '') 

        folder1_4 = QTreeWidgetItem(folder1, ["Documentation"])
        folder1_4.setData(2, Qt.EditRole, '') 


        folder2 = QTreeWidgetItem(root, ["Extra Stations"])
        folder2.setData(2, Qt.EditRole, '') 

        folder3 = QTreeWidgetItem(root, ["Output Files"])
        folder3.setData(2, Qt.EditRole, '') 

        folder4 = QTreeWidgetItem(root, ["Atmospheric"])
        folder4.setData(2, Qt.EditRole, '') 

        folder5 = QTreeWidgetItem(root, ["Weather Perturbations"])
        folder5.setData(2, Qt.EditRole, '') 


        docs_tab_content.addWidget(self.docTree)

    def initGeneralTab(self):
        general = QWidget()
        general_content = QGridLayout()
        general.setLayout(general_content)
        self.section_widget.addTab(general, "General")

        self.fireball_name_label = QLabel("Fireball Name:")
        self.fireball_name_edits = QLineEdit("")
        general_content.addWidget(self.fireball_name_label, 1, 1)
        general_content.addWidget(self.fireball_name_edits, 1, 2)
        self.fireball_name_label.setToolTip(self.toolTime('fireball_name'))

        self.difference_filter_label = QLabel("Difference Filter:")
        self.difference_filter_edits = QComboBox()
        general_content.addWidget(self.difference_filter_label, 2, 1)
        general_content.addWidget(self.difference_filter_edits, 2, 2)
        self.difference_filter_edits.addItem("True")
        self.difference_filter_edits.addItem("False")
        self.difference_filter_label.setToolTip(self.toolTime('difference_filter'))

        self.get_data_label = QLabel("Get Data:")
        self.get_data_edits = QComboBox()
        general_content.addWidget(self.get_data_label, 3, 1)
        general_content.addWidget(self.get_data_edits, 3, 2)
        self.get_data_edits.addItem("True")
        self.get_data_edits.addItem("False")
        self.get_data_label.setToolTip(self.toolTime('get_data'))

        self.run_mode_label = QLabel("Run Mode:")
        self.run_mode_edits = QComboBox()
        general_content.addWidget(self.run_mode_label, 4, 1)
        general_content.addWidget(self.run_mode_edits, 4, 2)
        self.run_mode_edits.addItem("Search")
        self.run_mode_edits.addItem("Replot")
        self.run_mode_edits.addItem("Manual")
        self.run_mode_label.setToolTip(self.toolTime('run_mode'))

        self.debug_label = QLabel("Debug:")
        self.debug_edits = QComboBox()
        general_content.addWidget(self.debug_label, 5, 1)
        general_content.addWidget(self.debug_edits, 5, 2)
        self.debug_edits.addItem("True")
        self.debug_edits.addItem("False")
        self.debug_label.setToolTip(self.toolTime('debug'))

    def initFilesTab(self):
        files = QWidget()
        files_content = QGridLayout()
        files.setLayout(files_content)
        self.section_widget.addTab(files, "Files")

        self.working_directory_label = QLabel("Working Directory:")
        self.working_directory_edits = QLineEdit("")
        self.working_directory_buton = QPushButton("Browse")
        files_content.addWidget(self.working_directory_label, 1, 1)
        files_content.addWidget(self.working_directory_edits, 1, 2)
        files_content.addWidget(self.working_directory_buton, 1, 3)
        self.working_directory_buton.clicked.connect(partial(self.folderSearch, self.working_directory_edits))
        self.working_directory_label.setToolTip(self.toolTime('working_directory'))

        self.arrival_times_label = QLabel("Arrival Times:")
        self.arrival_times_edits = QLineEdit("")
        self.arrival_times_buton = QPushButton("Browse")
        files_content.addWidget(self.arrival_times_label, 2, 1)
        files_content.addWidget(self.arrival_times_edits, 2, 2)
        files_content.addWidget(self.arrival_times_buton, 2, 3)
        self.arrival_times_buton.clicked.connect(partial(self.fileSearch, ['Numpy Array (*.npy)'], self.arrival_times_edits))
        self.arrival_times_label.setToolTip(self.toolTime('arrival_times_file'))

        self.sounding_file_label = QLabel("Sounding File:")
        self.sounding_file_edits = QLineEdit("")
        self.sounding_file_buton = QPushButton("Browse")
        files_content.addWidget(self.sounding_file_label, 3, 1)
        files_content.addWidget(self.sounding_file_edits, 3, 2)
        files_content.addWidget(self.sounding_file_buton, 3, 3)
        self.sounding_file_buton.clicked.connect(partial(self.fileSearch, ['NetCDF (*.nc)', 'HDF (*.HDF)'], self.sounding_file_edits))
        self.sounding_file_label.setToolTip(self.toolTime('sounding_file'))

        self.perturbation_file_label = QLabel("Perturbation File:")
        self.perturbation_file_edits = QLineEdit("")
        self.perturbation_file_buton = QPushButton("Browse")
        files_content.addWidget(self.perturbation_file_label, 4, 1)
        files_content.addWidget(self.perturbation_file_edits, 4, 2)
        files_content.addWidget(self.perturbation_file_buton, 4, 3)
        self.perturbation_file_buton.clicked.connect(partial(self.fileSearch, ['NetCDF (*.nc)'], self.perturbation_file_edits))
        self.perturbation_file_label.setToolTip(self.toolTime('perturbation_spread_file'))

        self.station_picks_file_label = QLabel("Station Picks File:")
        self.station_picks_file_edits = QLineEdit("")
        self.station_picks_file_buton = QPushButton("Browse")
        files_content.addWidget(self.station_picks_file_label, 5, 1)
        files_content.addWidget(self.station_picks_file_edits, 5, 2)
        files_content.addWidget(self.station_picks_file_buton, 5, 3)
        self.station_picks_file_buton.clicked.connect(partial(self.fileSearch, ['CSV (*.csv)', 'Text File (*.txt)'], self.station_picks_file_edits))
        self.station_picks_file_label.setToolTip(self.toolTime('station_picks_file'))

        self.points_name_label = QLabel("Replot Points File:")
        self.points_name_edits = QLineEdit("")
        self.points_name_buton = QPushButton("Browse")
        files_content.addWidget(self.points_name_label, 6, 1)
        files_content.addWidget(self.points_name_edits, 6, 2)
        files_content.addWidget(self.points_name_buton, 6, 3)
        self.points_name_buton.clicked.connect(partial(self.fileSearch, ['CSV (*.csv)'], self.points_name_edits))
        self.points_name_label.setToolTip(self.toolTime('points_name'))

    def initParametersTab(self):
        params = QWidget()
        params_content = QGridLayout()
        params.setLayout(params_content)
        self.section_widget.addTab(params, "Parameters")

        self.lat_centre_label = QLabel("Latitude Center:")
        self.lat_centre_edits = QLineEdit("")
        params_content.addWidget(self.lat_centre_label, 1, 1)
        params_content.addWidget(self.lat_centre_edits, 1, 2)
        self.lat_centre_label.setToolTip(self.toolTime('lat_centre'))

        self.lon_centre_label = QLabel("Longitude Center:")
        self.lon_centre_edits = QLineEdit("")
        params_content.addWidget(self.lon_centre_label, 2, 1)
        params_content.addWidget(self.lon_centre_edits, 2, 2)
        self.lon_centre_label.setToolTip(self.toolTime('lon_centre'))

        self.deg_radius_label = QLabel("Degrees in Search Radius:")
        self.deg_radius_edits = QLineEdit("")
        params_content.addWidget(self.deg_radius_label, 3, 1)
        params_content.addWidget(self.deg_radius_edits, 3, 2)
        self.deg_radius_label.setToolTip(self.toolTime('deg_radius'))

        self.start_datetime_label = QLabel("Start Datetime:")
        self.start_datetime_edits = QDateTimeEdit()
        params_content.addWidget(self.start_datetime_label, 4, 1)
        params_content.addWidget(self.start_datetime_edits, 4, 2)
        self.start_datetime_edits.setCalendarPopup(True)
        self.start_datetime_label.setToolTip(self.toolTime('start_datetime'))

        self.end_datetime_label = QLabel("End Datetime:")
        self.end_datetime_edits = QDateTimeEdit()
        params_content.addWidget(self.end_datetime_label, 5, 1)
        params_content.addWidget(self.end_datetime_edits, 5, 2)
        self.end_datetime_edits.setCalendarPopup(True)
        self.end_datetime_label.setToolTip(self.toolTime('end_datetime'))

        self.v_sound_label = QLabel("Average Speed of Sound:")
        self.v_sound_edits = QLineEdit("")
        params_content.addWidget(self.v_sound_label, 6, 1)
        params_content.addWidget(self.v_sound_edits, 6, 2)
        self.v_sound_label.setToolTip(self.toolTime('v_sound'))

    def initBallisticTab(self):

        ballistic = QWidget()
        ballistic_content = QGridLayout()
        ballistic.setLayout(ballistic_content)
        self.section_widget.addTab(ballistic, "Ballistic")

        self.t0_label = QLabel("t0:")
        self.t0_edits = QLineEdit("")
        ballistic_content.addWidget(self.t0_label, 1, 1)
        ballistic_content.addWidget(self.t0_edits, 1, 2, 1, 3)
        self.t0_label.setToolTip(self.toolTime('t0'))

        self.v_label = QLabel("v:")
        self.v_edits = QLineEdit("")
        ballistic_content.addWidget(self.v_label, 2, 1)
        ballistic_content.addWidget(self.v_edits, 2, 2, 1, 3)
        self.v_label.setToolTip(self.toolTime('v'))

        self.azim_label = QLabel("azim:")
        self.azim_edits = QLineEdit("")
        ballistic_content.addWidget(self.azim_label, 3, 1)
        ballistic_content.addWidget(self.azim_edits, 3, 2, 1, 3)
        self.azim_label.setToolTip(self.toolTime('azim'))

        self.zangle_label = QLabel("zangle:")
        self.zangle_edits = QLineEdit("")
        ballistic_content.addWidget(self.zangle_label, 4, 1)
        ballistic_content.addWidget(self.zangle_edits, 4, 2, 1, 3)
        self.zangle_label.setToolTip(self.toolTime('zangle'))

        self.lat_i_label = QLabel("lat_i:")
        self.lat_i_edits = QLineEdit("")
        ballistic_content.addWidget(self.lat_i_label, 5, 1)
        ballistic_content.addWidget(self.lat_i_edits, 5, 2)
        self.lat_i_label.setToolTip(self.toolTime('lat_i'))

        self.lon_i_label = QLabel("lon_i:")
        self.lon_i_edits = QLineEdit("")
        ballistic_content.addWidget(self.lon_i_label, 6, 1)
        ballistic_content.addWidget(self.lon_i_edits, 6, 2)
        self.lon_i_label.setToolTip(self.toolTime('lon_i'))

        self.elev_i_label = QLabel("elev_i:")
        self.elev_i_edits = QLineEdit("")
        ballistic_content.addWidget(self.elev_i_label, 7, 1)
        ballistic_content.addWidget(self.elev_i_edits, 7, 2)
        self.elev_i_label.setToolTip(self.toolTime('elev_i'))

        self.lat_f_label = QLabel("lat_f:")
        self.lat_f_edits = QLineEdit("")
        ballistic_content.addWidget(self.lat_f_label, 5, 3)
        ballistic_content.addWidget(self.lat_f_edits, 5, 4)
        self.lat_f_label.setToolTip(self.toolTime('lat_f'))

        self.lon_f_label = QLabel("lon_f:")
        self.lon_f_edits = QLineEdit("")
        ballistic_content.addWidget(self.lon_f_label, 6, 3)
        ballistic_content.addWidget(self.lon_f_edits, 6, 4)
        self.lon_f_label.setToolTip(self.toolTime('lon_f'))

        self.elev_f_label = QLabel("elev_f:")
        self.elev_f_edits = QLineEdit("")
        ballistic_content.addWidget(self.elev_f_label, 7, 3)
        ballistic_content.addWidget(self.elev_f_edits, 7, 4)
        self.elev_f_label.setToolTip(self.toolTime('elev_f'))

        self.show_ballistic_waveform_label = QLabel("Show Ballistic Waveform:")
        self.show_ballistic_waveform_edits = QComboBox()
        ballistic_content.addWidget(self.show_ballistic_waveform_label, 8, 1, 1, 2)
        ballistic_content.addWidget(self.show_ballistic_waveform_edits, 8, 3, 1, 2)
        self.show_ballistic_waveform_edits.addItem("True")
        self.show_ballistic_waveform_edits.addItem("False")
        self.show_ballistic_waveform_label.setToolTip(self.toolTime('show_ballistic_waveform'))

    def initFragmentationTab(self):
        fragmentation = QWidget()
        fragmentation_content = QGridLayout()
        fragmentation.setLayout(fragmentation_content)
        self.section_widget.addTab(fragmentation, "Fragmentation")

        
        self.fragmentation_point_label = QLabel("Fragmentation Point(s)")
        self.fragmentation_point = QTableWidget(0, 4)
        fragmentation_content.addWidget(self.fragmentation_point, 1, 1, 1, 4)
        fragmentation_content.addWidget(self.fragmentation_point_label, 0, 1, 1, 4)
        self.fragmentation_point.setHorizontalHeaderLabels(['Latitude', 'Longitude', 'Elevation', 'Time'])
        header = self.fragmentation_point.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        self.fragmentation_point_label.setToolTip(self.toolTime('fragmentation_point'))

        self.fragmentation_point_add = QPushButton("+")
        fragmentation_content.addWidget(self.fragmentation_point_add, 2, 3, 1, 2)
        self.fragmentation_point_add.clicked.connect(partial(self.changeRows, self.fragmentation_point, 1))
        self.fragmentation_point_add.setToolTip("Add row")

        self.fragmentation_point_min = QPushButton("-")
        fragmentation_content.addWidget(self.fragmentation_point_min, 2, 1, 1, 2)
        self.fragmentation_point_min.clicked.connect(partial(self.changeRows, self.fragmentation_point, -1))
        self.fragmentation_point_min.setToolTip("Remove row")

        self.show_fragmentation_waveform_label = QLabel("Show Fragmentation Waveform:")
        self.show_fragmentation_waveform_edits = QComboBox()
        fragmentation_content.addWidget(self.show_fragmentation_waveform_label, 3, 1, 1, 2)
        fragmentation_content.addWidget(self.show_fragmentation_waveform_edits, 3, 3, 1, 2)
        self.show_fragmentation_waveform_edits.addItem("True")
        self.show_fragmentation_waveform_edits.addItem("False")
        self.show_fragmentation_waveform_label.setToolTip("show_fragmentation_waveform")

        self.manual_label = QLabel("Manual Fragmentation Search:")
        fragmentation_content.addWidget(self.manual_label, 4, 1, 1, 4)
        self.manual_label.setToolTip("manual_fragmentation_search")

        self.lat_frag_label = QLabel("Latitude:")
        self.lat_frag_edits = QLineEdit("")
        fragmentation_content.addWidget(self.lat_frag_label, 5, 1)
        fragmentation_content.addWidget(self.lat_frag_edits, 6, 1)

        self.lon_frag_label = QLabel("Longitude:")
        self.lon_frag_edits = QLineEdit("")
        fragmentation_content.addWidget(self.lon_frag_label, 5, 2)
        fragmentation_content.addWidget(self.lon_frag_edits, 6, 2)

        self.elev_frag_label = QLabel("Elevation:")
        self.elev_frag_edits = QLineEdit("")
        fragmentation_content.addWidget(self.elev_frag_label, 5, 3)
        fragmentation_content.addWidget(self.elev_frag_edits, 6, 3)

        self.time_frag_label = QLabel("Time:")
        self.time_frag_edits = QLineEdit("")
        fragmentation_content.addWidget(self.time_frag_label, 5, 4)
        fragmentation_content.addWidget(self.time_frag_edits, 6, 4)

    def initRestrictionTab(self):
        restriction = QWidget()
        restriction_content = QGridLayout()
        restriction.setLayout(restriction_content)
        self.section_widget.addTab(restriction, "Restriction")

        self.v_fixed_label = QLabel("v_fixed:")
        self.v_fixed_edits = QLineEdit("")
        restriction_content.addWidget(self.v_fixed_label, 1, 1, 1, 2)
        restriction_content.addWidget(self.v_fixed_edits, 1, 3, 1, 2)
        self.v_fixed_label.setToolTip(self.toolTime('v_fixed'))

        self.max_error_label = QLabel("max_error:")
        self.max_error_edits = QLineEdit("")
        restriction_content.addWidget(self.max_error_label, 2, 1, 1, 2)
        restriction_content.addWidget(self.max_error_edits, 2, 3, 1, 2)
        self.max_error_label.setToolTip(self.toolTime('max_error'))

        self.restricted_time_label = QLabel("restricted_time:")
        self.restricted_time_edits = QDateTimeEdit()
        self.restricted_time_check = QCheckBox("Enable Restricted Time: ")
        restriction_content.addWidget(self.restricted_time_label, 3, 1, 1, 2)
        restriction_content.addWidget(self.restricted_time_edits, 3, 3, 1, 1)
        restriction_content.addWidget(self.restricted_time_check, 3, 4, 1, 1)
        self.restricted_time_edits.setCalendarPopup(True)
        self.restricted_time_label.setToolTip(self.toolTime('restricted_time'))

        self.traj_tol_label = QLabel("traj_tol:")
        self.traj_tol_edits = QLineEdit("")
        restriction_content.addWidget(self.traj_tol_label, 4, 1, 1, 2)
        restriction_content.addWidget(self.traj_tol_edits, 4, 3, 1, 2)
        self.traj_tol_label.setToolTip(self.toolTime('traj_tol'))

        self.restrict_to_trajectory_label = QLabel("restrict_to_trajectory:")
        self.restrict_to_trajectory_edits = QComboBox()
        restriction_content.addWidget(self.restrict_to_trajectory_label, 5, 1, 1, 2)
        restriction_content.addWidget(self.restrict_to_trajectory_edits, 5, 3, 1, 2)
        self.restrict_to_trajectory_edits.addItem("True")
        self.restrict_to_trajectory_edits.addItem("False")
        self.restrict_to_trajectory_label.setToolTip(self.toolTime('restrict_to_trajectory'))

        self.azimuth_min_label = QLabel("azimuth_min:")
        self.azimuth_min_edits = QLineEdit("")
        restriction_content.addWidget(self.azimuth_min_label, 6, 1, 1, 1)
        restriction_content.addWidget(self.azimuth_min_edits, 6, 2, 1, 1)
        self.azimuth_min_label.setToolTip(self.toolTime('azimuth_min'))

        self.azimuth_max_label = QLabel("azimuth_max:")
        self.azimuth_max_edits = QLineEdit("")
        restriction_content.addWidget(self.azimuth_max_label, 6, 3, 1, 1)
        restriction_content.addWidget(self.azimuth_max_edits, 6, 4, 1, 1)
        self.azimuth_max_label.setToolTip(self.toolTime('azimuth_max'))

        self.zangle_min_label = QLabel("zangle_min:")
        self.zangle_min_edits = QLineEdit("")
        restriction_content.addWidget(self.zangle_min_label, 7, 1, 1, 1)
        restriction_content.addWidget(self.zangle_min_edits, 7, 2, 1, 1)
        self.zangle_min_label.setToolTip(self.toolTime('zenith_min'))

        self.zangle_max_label = QLabel("zangle_max:")
        self.zangle_max_edits = QLineEdit("")
        restriction_content.addWidget(self.zangle_max_label, 7, 3, 1, 1)
        restriction_content.addWidget(self.zangle_max_edits, 7, 4, 1, 1)
        self.zangle_max_label.setToolTip(self.toolTime('zenith_max'))

        self.x_min_label = QLabel("x_min:")
        self.x_min_edits = QLineEdit("")
        restriction_content.addWidget(self.x_min_label, 8, 1, 1, 1)
        restriction_content.addWidget(self.x_min_edits, 8, 2, 1, 1)
        self.x_min_label.setToolTip(self.toolTime('x_min'))

        self.x_max_label = QLabel("x_max:")
        self.x_max_edits = QLineEdit("")
        restriction_content.addWidget(self.x_max_label, 8, 3, 1, 1)
        restriction_content.addWidget(self.x_max_edits, 8, 4, 1, 1)
        self.x_max_label.setToolTip(self.toolTime('x_max'))

        self.y_min_label = QLabel("y_min:")
        self.y_min_edits = QLineEdit("")
        restriction_content.addWidget(self.y_min_label, 9, 1, 1, 1)
        restriction_content.addWidget(self.y_min_edits, 9, 2, 1, 1)
        self.y_min_label.setToolTip(self.toolTime('y_min'))

        self.y_max_label = QLabel("y_max:")
        self.y_max_edits = QLineEdit("")
        restriction_content.addWidget(self.y_max_label, 9, 3, 1, 1)
        restriction_content.addWidget(self.y_max_edits, 9, 4, 1, 1)
        self.y_max_label.setToolTip(self.toolTime('y_max'))

        self.t_min_label = QLabel("t_min:")
        self.t_min_edits = QLineEdit("")
        restriction_content.addWidget(self.t_min_label, 10, 1, 1, 1)
        restriction_content.addWidget(self.t_min_edits, 10, 2, 1, 1)
        self.t_min_label.setToolTip(self.toolTime('t_min'))

        self.t_max_label = QLabel("t_max:")
        self.t_max_edits = QLineEdit("")
        restriction_content.addWidget(self.t_max_label, 10, 3, 1, 1)
        restriction_content.addWidget(self.t_max_edits, 10, 4, 1, 1)
        self.t_max_label.setToolTip(self.toolTime('t_max'))

        self.v_min_label = QLabel("v_min:")
        self.v_min_edits = QLineEdit("")
        restriction_content.addWidget(self.v_min_label, 11, 1, 1, 1)
        restriction_content.addWidget(self.v_min_edits, 11, 2, 1, 1)
        self.v_min_label.setToolTip(self.toolTime('v_min'))

        self.v_max_label = QLabel("v_max:")
        self.v_max_edits = QLineEdit("")
        restriction_content.addWidget(self.v_max_label, 11, 3, 1, 1)
        restriction_content.addWidget(self.v_max_edits, 11, 4, 1, 1)
        self.v_max_label.setToolTip(self.toolTime('v_max'))

        self.weight_distance_min_label = QLabel("weight_distance_min:")
        self.weight_distance_min_edits = QLineEdit("")
        restriction_content.addWidget(self.weight_distance_min_label, 12, 1, 1, 1)
        restriction_content.addWidget(self.weight_distance_min_edits, 12, 2, 1, 1)
        self.weight_distance_min_label.setToolTip(self.toolTime('weight_distance_min'))

        self.weight_distance_max_label = QLabel("weight_distance_max:")
        self.weight_distance_max_edits = QLineEdit("")
        restriction_content.addWidget(self.weight_distance_max_label, 12, 3, 1, 1)
        restriction_content.addWidget(self.weight_distance_max_edits, 12, 4, 1, 1)
        self.weight_distance_max_label.setToolTip(self.toolTime('weight_distance_max'))

        self.search_time_min_label = QLabel("search_time_min:")
        self.search_time_min_edits = QLineEdit("")
        restriction_content.addWidget(self.search_time_min_label, 13, 1, 1, 1)
        restriction_content.addWidget(self.search_time_min_edits, 13, 2, 1, 1)
        self.search_time_min_label.setToolTip(self.toolTime('min_time'))

        self.search_time_max_label = QLabel("search_time_max:")
        self.search_time_max_edits = QLineEdit("")
        restriction_content.addWidget(self.search_time_max_label, 13, 3, 1, 1)
        restriction_content.addWidget(self.search_time_max_edits, 13, 4, 1, 1)
        self.search_time_max_label.setToolTip(self.toolTime('max_time'))

        self.search_lat_min_label = QLabel("search_lat_min:")
        self.search_lat_min_edits = QLineEdit("")
        restriction_content.addWidget(self.search_lat_min_label, 14, 1, 1, 1)
        restriction_content.addWidget(self.search_lat_min_edits, 14, 2, 1, 1)
        self.search_lat_min_label.setToolTip(self.toolTime('search_area'))

        self.search_lat_max_label = QLabel("search_lat_max:")
        self.search_lat_max_edits = QLineEdit("")
        restriction_content.addWidget(self.search_lat_max_label, 14, 3, 1, 1)
        restriction_content.addWidget(self.search_lat_max_edits, 14, 4, 1, 1)
        self.search_lat_max_label.setToolTip(self.toolTime('search_area'))

        self.search_lon_min_label = QLabel("search_lon_min:")
        self.search_lon_min_edits = QLineEdit("")
        restriction_content.addWidget(self.search_lon_min_label, 15, 1, 1, 1)
        restriction_content.addWidget(self.search_lon_min_edits, 15, 2, 1, 1)
        self.search_lon_min_label.setToolTip(self.toolTime('search_area'))

        self.search_lon_max_label = QLabel("search_lon_max:")
        self.search_lon_max_edits = QLineEdit("")
        restriction_content.addWidget(self.search_lon_max_label, 15, 3, 1, 1)
        restriction_content.addWidget(self.search_lon_max_edits, 15, 4, 1, 1)
        self.search_lon_max_label.setToolTip(self.toolTime('search_area'))

        self.search_elev_min_label = QLabel("search_elev_min:")
        self.search_elev_min_edits = QLineEdit("")
        restriction_content.addWidget(self.search_elev_min_label, 16, 1, 1, 1)
        restriction_content.addWidget(self.search_elev_min_edits, 16, 2, 1, 1)
        self.search_elev_min_label.setToolTip(self.toolTime('search_area'))

        self.search_elev_max_label = QLabel("search_elev_max:")
        self.search_elev_max_edits = QLineEdit("")
        restriction_content.addWidget(self.search_elev_max_label, 16, 3, 1, 1)
        restriction_content.addWidget(self.search_elev_max_edits, 16, 4, 1, 1)
        self.search_elev_max_label.setToolTip(self.toolTime('search_area'))

    def initAtmosphereTab(self):

        atmosphere = QWidget()
        atmosphere_content = QGridLayout()
        atmosphere.setLayout(atmosphere_content)
        self.section_widget.addTab(atmosphere, "Atmosphere")

        self.enable_winds_label = QLabel("Enable Winds:")
        self.enable_winds_edits = QComboBox()
        atmosphere_content.addWidget(self.enable_winds_label, 1, 1)
        atmosphere_content.addWidget(self.enable_winds_edits, 1, 2)
        self.enable_winds_edits.addItem("True")
        self.enable_winds_edits.addItem("False")
        self.enable_winds_label.setToolTip(self.toolTime('enable_winds'))

        self.weather_type_label = QLabel("Weather Type:")
        self.weather_type_edits = QComboBox()
        atmosphere_content.addWidget(self.weather_type_label, 2, 1)
        atmosphere_content.addWidget(self.weather_type_edits, 2, 2)
        self.weather_type_edits.addItem("none")
        self.weather_type_edits.addItem("ecmwf")
        self.weather_type_edits.addItem("ukmo")
        self.weather_type_edits.addItem("merra")
        self.weather_type_edits.addItem("custom")
        self.weather_type_label.setToolTip(self.toolTime('weather_type'))

        self.grid_size_label = QLabel("Grid Size:")
        self.grid_size_edits = QLineEdit("")
        atmosphere_content.addWidget(self.grid_size_label, 3, 1)
        atmosphere_content.addWidget(self.grid_size_edits, 3, 2)
        self.grid_size_label.setToolTip(self.toolTime('grid_size'))

    def initPerturbationsTab(self):

        perturb = QWidget()
        perturb_content = QGridLayout()
        perturb.setLayout(perturb_content)
        self.section_widget.addTab(perturb, "Perturbations")

        self.perturb_times_label = QLabel("Perturbation Times:")
        self.perturb_times_edits = QLineEdit("")
        perturb_content.addWidget(self.perturb_times_label, 1, 1)
        perturb_content.addWidget(self.perturb_times_edits, 1, 2)
        self.perturb_times_label.setToolTip(self.toolTime('perturb_times'))

        self.frag_no_label = QLabel("Fragmentation Number:")
        self.frag_no_edits = QLineEdit("")
        perturb_content.addWidget(self.frag_no_label, 2, 1)
        perturb_content.addWidget(self.frag_no_edits, 2, 2)
        self.frag_no_label.setToolTip(self.toolTime('fragno'))

        self.perturb_label = QLabel("Perturb:")
        self.perturb_edits = QComboBox()
        perturb_content.addWidget(self.perturb_label, 3, 1)
        perturb_content.addWidget(self.perturb_edits, 3, 2)
        self.perturb_edits.addItem("True")
        self.perturb_edits.addItem("False")
        self.perturb_label.setToolTip(self.toolTime('perturb'))

        self.perturb_method_label = QLabel("Perturb Method:")
        self.perturb_method_edits = QComboBox()
        perturb_content.addWidget(self.perturb_method_label, 4, 1)
        perturb_content.addWidget(self.perturb_method_edits, 4, 2)
        self.perturb_method_edits.addItem("none")
        self.perturb_method_edits.addItem("bmp")
        self.perturb_method_edits.addItem("sppt")
        self.perturb_method_edits.addItem("temporal")
        self.perturb_method_edits.addItem("spread")
        self.perturb_method_edits.addItem("spread_r")
        self.perturb_method_label.setToolTip(self.toolTime('perturb_method'))

    def initSpeedTab(self):

        speed = QWidget()
        speed_content = QGridLayout()
        speed.setLayout(speed_content)
        self.section_widget.addTab(speed, "Speed")

        self.fast_ballistic_label = QLabel("Fast Ballistic:")
        self.fast_ballistic_edits = QComboBox()
        speed_content.addWidget(self.fast_ballistic_label, 1, 1)
        speed_content.addWidget(self.fast_ballistic_edits, 1, 2)
        self.fast_ballistic_edits.addItem("True")
        self.fast_ballistic_edits.addItem("False")
        self.fast_ballistic_label.setToolTip(self.toolTime('fast_ballistic')) 

        self.fit_type_label = QLabel("Fit Type:")
        self.fit_type_edits = QLineEdit("")
        speed_content.addWidget(self.fit_type_label, 2, 1)
        speed_content.addWidget(self.fit_type_edits, 2, 2)
        self.fit_type_label.setToolTip(self.toolTime('fit_type'))

        self.n_theta_label = QLabel("Theta Resolution:")
        self.n_theta_edits = QLineEdit("")
        speed_content.addWidget(self.n_theta_label, 3, 1)
        speed_content.addWidget(self.n_theta_edits, 3, 2)
        self.n_theta_label.setToolTip(self.toolTime('n_theta'))

        self.n_phi_label = QLabel("Phi Resolution:")
        self.n_phi_edits = QLineEdit("")
        speed_content.addWidget(self.n_phi_label, 4, 1)
        speed_content.addWidget(self.n_phi_edits, 4, 2)
        self.n_phi_label.setToolTip(self.toolTime('n_phi'))

        self.angle_precision_label = QLabel("Angle Precision:")
        self.angle_precision_edits = QLineEdit("")
        speed_content.addWidget(self.angle_precision_label, 5, 1)
        speed_content.addWidget(self.angle_precision_edits, 5, 2)
        self.angle_precision_label.setToolTip(self.toolTime('angle_precision'))

        self.angle_error_tol_label = QLabel("Angle Error Tolerance:")
        self.angle_error_tol_edits = QLineEdit("")
        speed_content.addWidget(self.angle_error_tol_label, 6, 1)
        speed_content.addWidget(self.angle_error_tol_edits, 6, 2)
        self.angle_error_tol_label.setToolTip(self.toolTime('angle_error_tol'))

    def initPSOTab(self):
        pso = QWidget()
        pso_content = QGridLayout()
        pso.setLayout(pso_content)
        self.section_widget.addTab(pso, "PSO")

        self.maxiter_label = QLabel("Max Iterations:")
        self.maxiter_edits = QLineEdit("50")
        pso_content.addWidget(self.maxiter_label, 1, 1)
        pso_content.addWidget(self.maxiter_edits, 1, 2)
        self.maxiter_label.setToolTip(self.toolTime('maxiter'))

        self.swarmsize_label = QLabel("Swarm Size:")
        self.swarmsize_edits = QLineEdit("250")
        pso_content.addWidget(self.swarmsize_label, 2, 1)
        pso_content.addWidget(self.swarmsize_edits, 2, 2)
        self.swarmsize_label.setToolTip(self.toolTime('swarmsize'))

        self.run_times_label = QLabel("Run Times:")
        self.run_times_edits = QLineEdit("1")
        pso_content.addWidget(self.run_times_label, 3, 1)
        pso_content.addWidget(self.run_times_edits, 3, 2)
        self.run_times_label.setToolTip(self.toolTime('run_times'))

        self.minfunc_label = QLabel("minfunc:")
        self.minfunc_edits = QLineEdit("1e-8")
        pso_content.addWidget(self.minfunc_label, 4, 1)
        pso_content.addWidget(self.minfunc_edits, 4, 2)
        self.minfunc_label.setToolTip(self.toolTime('minfunc'))

        self.minstep_label = QLabel("minstep:")
        self.minstep_edits = QLineEdit("1e-8")
        pso_content.addWidget(self.minstep_label, 5, 1)
        pso_content.addWidget(self.minstep_edits, 5, 2)
        self.minstep_label.setToolTip(self.toolTime('minstep'))

        self.phip_label = QLabel("phip:")
        self.phip_edits = QLineEdit("0.5")
        pso_content.addWidget(self.phip_label, 6, 1)
        pso_content.addWidget(self.phip_edits, 6, 2)
        self.phip_label.setToolTip(self.toolTime('phip'))

        self.phig_label = QLabel("phig:")
        self.phig_edits = QLineEdit("0.5")
        pso_content.addWidget(self.phig_label, 7, 1)
        pso_content.addWidget(self.phig_edits, 7, 2)
        self.phig_label.setToolTip(self.toolTime('phig'))

        self.omega_label = QLabel("omega:")
        self.omega_edits = QLineEdit("0.5")
        pso_content.addWidget(self.omega_label, 8, 1)
        pso_content.addWidget(self.omega_edits, 8, 2)
        self.omega_label.setToolTip(self.toolTime('omega'))

        self.pso_debug_label = QLabel("PSO Debug:")
        self.pso_debug_edits = QComboBox()
        pso_content.addWidget(self.pso_debug_label, 9, 1)
        pso_content.addWidget(self.pso_debug_edits, 9, 2)
        self.pso_debug_edits.addItem("True")
        self.pso_debug_edits.addItem("False")
        self.pso_debug_label.setToolTip(self.toolTime('pso_debug'))

    def initGraphingTab(self):

        graphing = QWidget()
        graphing_content = QGridLayout()
        graphing.setLayout(graphing_content)
        self.section_widget.addTab(graphing, "Graphing")

        self.plot_all_stations_label = QLabel("Plot All Stations:")
        self.plot_all_stations_edits = QComboBox()
        graphing_content.addWidget(self.plot_all_stations_label, 1, 1)
        graphing_content.addWidget(self.plot_all_stations_edits, 1, 2)
        self.plot_all_stations_edits.addItem("True")
        self.plot_all_stations_edits.addItem("False")
        self.plot_all_stations_label.setToolTip(self.toolTime('plot_all_stations'))

        self.color_toggle_label = QLabel("Toggle Color:")
        self.color_toggle_edits = QComboBox()
        graphing_content.addWidget(self.color_toggle_label, 2, 1)
        graphing_content.addWidget(self.color_toggle_edits, 2, 2)
        self.color_toggle_edits.addItem("True")
        self.color_toggle_edits.addItem("False")
        self.color_toggle_label.setToolTip(self.toolTime('colortoggle'))

        self.dot_tol_label = QLabel("Dot Product Tolerance:")
        self.dot_tol_edits = QLineEdit("")
        graphing_content.addWidget(self.dot_tol_label, 3, 1)
        graphing_content.addWidget(self.dot_tol_edits, 3, 2)
        self.dot_tol_label.setToolTip(self.toolTime('dot_tol'))

        self.contour_res_label = QLabel("Contour Resolution:")
        self.contour_res_edits = QLineEdit("")
        graphing_content.addWidget(self.contour_res_label, 4, 1)
        graphing_content.addWidget(self.contour_res_edits, 4, 2)
        self.contour_res_label.setToolTip(self.toolTime('contour_res'))

        self.high_f_label = QLabel("Highlight Fragmentation:")
        self.high_f_edits = QLineEdit("")
        graphing_content.addWidget(self.high_f_label, 5, 1)
        graphing_content.addWidget(self.high_f_edits, 5, 2)
        self.high_f_label.setToolTip(self.toolTime('high_f'))

        self.high_b_label = QLabel("Highlight Ballistic:")
        self.high_b_edits = QLineEdit("")
        graphing_content.addWidget(self.high_b_label, 6, 1)
        graphing_content.addWidget(self.high_b_edits, 6, 2)
        self.high_b_label.setToolTip(self.toolTime('high_b'))

        self.rm_stat_label = QLabel("Remove Stations:")
        self.rm_stat_edits = QLineEdit("")
        graphing_content.addWidget(self.rm_stat_label, 7, 1)
        graphing_content.addWidget(self.rm_stat_edits, 7, 2)
        self.rm_stat_label.setToolTip(self.toolTime('rm_stat'))

        self.img_dim_label = QLabel("Image Dimensions:")
        self.img_dim_edits = QLineEdit("")
        graphing_content.addWidget(self.img_dim_label, 8, 1)
        graphing_content.addWidget(self.img_dim_edits, 8, 2)
        self.img_dim_label.setToolTip(self.toolTime('img_dim'))

        self.reported_points_label = QLabel("Reported Points")
        self.reported_points = QTableWidget(0, 4)
        graphing_content.addWidget(self.reported_points, 10, 1, 1, 2)
        graphing_content.addWidget(self.reported_points_label, 9, 1, 1, 4)
        self.reported_points.setHorizontalHeaderLabels(['Latitude', 'Longitude', 'Elevation', 'Time'])
        header = self.reported_points.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        self.reported_points_label.setToolTip(self.toolTime('reported_points'))

        self.reported_points_add = QPushButton("+")
        graphing_content.addWidget(self.reported_points_add, 11, 2)
        self.reported_points_add.clicked.connect(partial(self.changeRows, self.reported_points, 1))

        self.reported_points_min = QPushButton("-")
        graphing_content.addWidget(self.reported_points_min, 11, 1)
        self.reported_points_min.clicked.connect(partial(self.changeRows, self.reported_points, -1))

    def initExtraStationsTab(self):
        extra = QWidget()
        extra_content = QGridLayout()
        extra.setLayout(extra_content)
        self.section_widget.addTab(extra, "Extra Stations")

        self.extra_label = QLabel("Manually Added Stations")
        self.extra_point = QTableWidget(0, 8)
        extra_content.addWidget(self.extra_point, 1, 1, 1, 4)
        extra_content.addWidget(self.extra_label, 0, 1, 1, 4)
        self.extra_point.setHorizontalHeaderLabels(\
            ['Station Network', 'Station Code', 'Latitude', 'Longitude', \
            'Elevation', 'Station Display Name', 'Channel', 'File Name'])
        header = self.extra_point.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)
        self.extra_label.setToolTip(self.toolTime('stations'))

        self.extra_point_add = QPushButton("+")
        extra_content.addWidget(self.extra_point_add, 2, 3, 1, 2)
        self.extra_point_add.clicked.connect(partial(self.changeRows, self.extra_point, 1))

        self.extra_point_min = QPushButton("-")
        extra_content.addWidget(self.extra_point_min, 2, 1, 1, 2)
        self.extra_point_min.clicked.connect(partial(self.changeRows, self.extra_point, -1))

    def addIniWidgets(self):

        ini_builder_tab = QWidget()
        ini_builder_tab_content = QVBoxLayout()
        ini_builder_tab.setLayout(ini_builder_tab_content)
        button_panel = QHBoxLayout()

        self.load_button = QPushButton('Load')
        button_panel.addWidget(self.load_button)
        self.load_button.clicked.connect(self.loadINI)

        self.save_button = QPushButton('Save')
        button_panel.addWidget(self.save_button)
        self.save_button.clicked.connect(partial(self.saveINI, True))

        self.tab_widget.addTab(ini_builder_tab, "Ini Builder")

        self.section_widget = QTabWidget()
        
        self.initGeneralTab()
        self.initFilesTab()
        self.initParametersTab()
        self.initBallisticTab()
        self.initFragmentationTab()
        self.initRestrictionTab()
        self.initAtmosphereTab()
        self.initPerturbationsTab()
        self.initSpeedTab()
        self.initPSOTab()
        self.initGraphingTab()
        self.initExtraStationsTab()

        ini_builder_tab_content.addWidget(self.section_widget)
        ini_builder_tab_content.addLayout(button_panel)

    def atmPlotProfile(self, lat, lon, var_typ='t', perturb='none'):
        
        consts = Constants()
        setup = self.saveINI(False)

        setup.sounding_file = self.atm_atm_file_edits.text()
        setup.weather_type = self.atm_weather_type_edits.currentText()
        setup.perturbation_spread_file = self.atm_perturbation_file_edits.text()
        try:
            setup.perturb_times = int(self.atm_perturb_times_edits.text())
        except:
            setup.perturb_times = 0
        setup.perturb_method = self.atm_perturb_method_edits.currentText()

        if setup.weather_type == 'none':
            self.errorMessage('Weather type is set to "none", no weather can be displayed', 1)
            return None

        try:
            dataset = parseWeather(setup, consts)
            sounding = findECMWFSound(lat, lon, dataset)
        except:
            self.errorMessage('Error reading weather profile in atmPlotProfile', 2)
            return None

        self.var_typ = var_typ
        self.atm_canvas.setLabel('left', "Height", units='m')
        if self.var_typ == 't':
            X = sounding[:, 1]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Temperature", units='K')
        elif self.var_typ == 'm':
            X = sounding[:, 2]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Wind Magnitude", units='m/s')
        elif self.var_typ == 'd':
            X = sounding[:, 3]
            Y = sounding[:, 0]
            self.atm_canvas.setLabel('bottom', "Wind Direction", units='deg from N')
        else:
            self.errorMessage('Error reading var_typ in atmPlotProfile', 2)
            return None

        self.atm_canvas.clear()
        self.atm_canvas.plot(x=X, y=Y, pen='w')
        SolutionGUI.update(self)

        if setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(setup, consts, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(setup, consts, time=-1)

        else:
            sounding_u = []
            sounding_l = []

        if setup.perturb_method != 'none':
            for ptb_n in range(setup.perturb_times):

                if ptb_n > 0:
                    
                    if setup.debug:
                        print("STATUS: Perturbation {:}".format(ptb_n))

                    # generate a perturbed sounding profile
                    sounding_p = perturbation_method(setup, dataset, setup.perturb_method, \
                        sounding_u=sounding_u, sounding_l=sounding_l, \
                        spread_file=setup.perturbation_spread_file, lat=setup.lat_centre, lon=setup.lon_centre)
                    sounding_p = findECMWFSound(lat, lon, sounding_p)
                    
                    if self.var_typ == 't':
                        X = sounding_p[:, 1]
                        Y = sounding_p[:, 0]
                    elif self.var_typ == 'm':
                        X = sounding_p[:, 2]
                        Y = sounding_p[:, 0]
                    elif self.var_typ == 'd':
                        X = sounding_p[:, 3]
                        Y = sounding_p[:, 0]
                    else:
                        print('error, atmPlotProfile')

                    self.atm_canvas.plot(x=X, y=Y, pen='g')
                    SolutionGUI.update(self)

    def atmValueChange(self, obj, slider):
        
        if obj == self.atm_lat_label:
            obj.setText('Latitude: {:8.2f}'.format(slider.value()*self.slider_scale))
        elif obj == self.atm_lon_label:
            obj.setText('Longitude: {:8.2f}'.format(slider.value()*self.slider_scale))
        else:
            self.errorMessage(self, 'Bad atm slider pass in atmValueChange', 2)

        self.atmPlotProfile(self.atm_lat_slide.value()*self.slider_scale, self.atm_lon_slide.value()*self.slider_scale, self.var_typ)
    
    def addProfileWidgets(self):
        profile_tab = QWidget()
        profile_master = QHBoxLayout()
        profile_tab_content = QGridLayout()
        profile_tab_content_graph = QVBoxLayout()
        profile_tab.setLayout(profile_master)
        profile_master.addLayout(profile_tab_content_graph)
        profile_master.addLayout(profile_tab_content)

        self.atm_lat_label = QLabel("Latitude: 0")
        self.atm_lat_slide = QSlider(Qt.Horizontal)
        profile_tab_content.addWidget(self.atm_lat_label, 6, 1)
        profile_tab_content.addWidget(self.atm_lat_slide, 6, 2, 1, 3)
        self.atm_lat_slide.setMinimum(-90/self.slider_scale)
        self.atm_lat_slide.setMaximum(90/self.slider_scale)
        self.atm_lat_slide.setValue(0)
        self.atm_lat_slide.setTickInterval(0.5)
        self.atm_lat_slide.valueChanged.connect(partial(self.atmValueChange, self.atm_lat_label, self.atm_lat_slide))

        self.atm_lon_label = QLabel("Longitude: 0")
        self.atm_lon_slide = QSlider(Qt.Horizontal)
        profile_tab_content.addWidget(self.atm_lon_label, 7, 1)
        profile_tab_content.addWidget(self.atm_lon_slide, 7, 2, 1, 3)
        self.atm_lon_slide.setMinimum(-180/self.slider_scale)
        self.atm_lon_slide.setMaximum(180/self.slider_scale)
        self.atm_lon_slide.setValue(0)
        self.atm_lon_slide.setTickInterval(0.5)
        self.atm_lon_slide.valueChanged.connect(partial(self.atmValueChange, self.atm_lon_label, self.atm_lon_slide))

        self.atm_view = pg.GraphicsLayoutWidget()
        self.atm_canvas = self.atm_view.addPlot()
        #self.atm_ax = pg.PlotItem(size=10, title='Atmosphere')
        #self.atm_canvas.addItem(self.atm_ax)

        profile_tab_content_graph.addWidget(self.atm_view)
        self.atm_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        self.atm_T_button = QPushButton('Temperature')
        profile_tab_content_graph.addWidget(self.atm_T_button)
        self.atm_T_button.clicked.connect(partial(self.atmPlotProfile, self.atm_lat_slide.value()*self.slider_scale, self.atm_lon_slide.value()*self.slider_scale, var_typ='t'))

        self.atm_mag_button = QPushButton('Wind Magnitude')
        profile_tab_content_graph.addWidget(self.atm_mag_button)
        self.atm_mag_button.clicked.connect(partial(self.atmPlotProfile, self.atm_lat_slide.value()*self.slider_scale, self.atm_lon_slide.value()*self.slider_scale, var_typ='m'))

        self.atm_dir_button = QPushButton('Wind Direction')
        profile_tab_content_graph.addWidget(self.atm_dir_button)
        self.atm_dir_button.clicked.connect(partial(self.atmPlotProfile, self.atm_lat_slide.value()*self.slider_scale, self.atm_lon_slide.value()*self.slider_scale, var_typ='d'))

        self.atm_atm_file_label = QLabel("Atmospheric File:")
        self.atm_atm_file_edits = QLineEdit("")
        self.atm_atm_file_buton = QPushButton("Browse")
        profile_tab_content.addWidget(self.atm_atm_file_label, 1, 1)
        profile_tab_content.addWidget(self.atm_atm_file_edits, 1, 2, 1, 2)
        profile_tab_content.addWidget(self.atm_atm_file_buton, 1, 4)
        self.atm_atm_file_buton.clicked.connect(partial(self.fileSearch, ['NetCDF (*.nc)', 'HDF (*.HDF)'], self.atm_atm_file_edits))

        self.atm_weather_type_label = QLabel("Weather Type:")
        self.atm_weather_type_edits = QComboBox()
        profile_tab_content.addWidget(self.atm_weather_type_label, 2, 1)
        profile_tab_content.addWidget(self.atm_weather_type_edits, 2, 2, 1, 3)
        self.atm_weather_type_edits.addItem("none")
        self.atm_weather_type_edits.addItem("ecmwf")
        self.atm_weather_type_edits.addItem("ukmo")
        self.atm_weather_type_edits.addItem("merra")
        self.atm_weather_type_edits.addItem("custom")

        self.atm_perturbation_file_label = QLabel("Perturbation File:")
        self.atm_perturbation_file_edits = QLineEdit("")
        self.atm_perturbation_file_buton = QPushButton("Browse")
        profile_tab_content.addWidget(self.atm_perturbation_file_label, 3, 1)
        profile_tab_content.addWidget(self.atm_perturbation_file_edits, 3, 2, 1, 2)
        profile_tab_content.addWidget(self.atm_perturbation_file_buton, 3, 4)
        self.atm_perturbation_file_buton.clicked.connect(partial(self.fileSearch, ['NetCDF (*.nc)'], self.atm_perturbation_file_edits))

        self.atm_perturb_times_label = QLabel("Perturbation Times:")
        self.atm_perturb_times_edits = QLineEdit("")
        profile_tab_content.addWidget(self.atm_perturb_times_label, 4, 1)
        profile_tab_content.addWidget(self.atm_perturb_times_edits, 4, 2, 1, 3)

        self.atm_perturb_method_label = QLabel("Perturb Method:")
        self.atm_perturb_method_edits = QComboBox()
        profile_tab_content.addWidget(self.atm_perturb_method_label, 5, 1)
        profile_tab_content.addWidget(self.atm_perturb_method_edits, 5, 2, 1, 3)
        self.atm_perturb_method_edits.addItem("none")
        self.atm_perturb_method_edits.addItem("bmp")
        self.atm_perturb_method_edits.addItem("sppt")
        self.atm_perturb_method_edits.addItem("temporal")
        self.atm_perturb_method_edits.addItem("spread")
        self.atm_perturb_method_edits.addItem("spread_r")

        self.tab_widget.addTab(profile_tab, "Atmospheric Profile")

    def makeStationObj(self, lst):

        new_lst = []

        for line in lst:

            pos = position(line[2], line[3], line[4])
            stn = station(line[0], line[1], pos, line[5], line[6], line[7])
            new_lst.append(stn)

        return new_lst
    
    def makePicks(self):

        setup  = self.saveINI(write=False)
        
        if not os.path.exists(setup.working_directory):
            os.makedirs(setup.working_directory)

            #Build seismic data path
        dir_path = os.path.join(setup.working_directory, setup.fireball_name)

        # Load the station and waveform files list
        data_file_path = os.path.join(dir_path, DATA_FILE)
        if os.path.isfile(data_file_path):
            
            stn_list = readStationAndWaveformsListFile(data_file_path, rm_stat=setup.rm_stat)

        else:
            print('Station and waveform data file not found! Download the waveform files first!')
            sys.exit()

        if setup.stations is not None:
            stn_list = stn_list + self.makeStationObj(setup.stations)

        # Init the constants
        consts = Constants()
        setup.search_area = [0, 0, 0, 0]
        setup.search_area[0] = setup.lat_centre - setup.deg_radius 
        setup.search_area[1] = setup.lat_centre + setup.deg_radius
        setup.search_area[2] = setup.lon_centre - setup.deg_radius
        setup.search_area[3] = setup.lon_centre + setup.deg_radius

        sounding = parseWeather(setup, consts)
        arrTimes = np.array([])

        if len(stn_list) == 0:
            self.errorMessage('No Stations to load', 2)
            return None 

        try:
            #turn coordinates into position objects
            setup.traj_i = position(setup.lat_i, setup.lon_i, setup.elev_i)
            setup.traj_f = position(setup.lat_f, setup.lon_f, setup.elev_f)
        except:
            setup.traj_i = position(0, 0, 0)
            setup.traj_f = position(0, 0, 0)
            print("Warning: Unable to build trajectory points")

        self.waveformPicker(dir_path, setup, sounding, stn_list, difference_filter_all=setup.difference_filter_all, stn_list=stn_list)
    
    def trajCalc(setup):
        """ Creates trajectory between point A and the ground (B) based off of the initial position and the angle of travel

        Arguments:
        setup: [Object] ini file parameters

        Returns: 
        A [list] lat/lon/elev of the tail of the trajectory
        B [list] lat/lon/elev of the head of the trajectory
        """
        print("Trajectory Calculation...")

        ref = position(setup.lat_centre, setup.lon_centre, 0)

        # Tail of the trajectory
        A = geo2Loc(ref.lat, ref.lon, ref.elev, setup.lat_i, setup.lon_i, setup.elev_i*1000)

        # convert angles to radians
        ze = np.radians(setup.zangle)
        az = np.radians(setup.azim)

        # Create trajectory vector
        traj = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

        # How far along the trajectory until it reaches the ground
        n = -A[2]/traj[2]

        # B is the intersection between the trajectory vector and the ground
        B = A + n*traj

        # Convert back to geo coordinates
        B = np.array(loc2Geo(ref.lat, ref.lon, ref.elev, B))
        A = np.array(loc2Geo(ref.lat, ref.lon, ref.elev, A))

        print("Created Trajectory between A and B:")
        print("     A = {:10.4f}N {:10.4f}E {:10.2f}m".format(A[0], A[1], A[2]))
        print("     B = {:10.4f}N {:10.4f}E {:10.2f}m".format(B[0], B[1], B[2]))

        return A, B


    def calcAllTimes(self, stn_list, setup, sounding):
        """ 
        Method to calculate all arrival times and place them into an array, so that they are not recalculated every
        time a new waveform is opened
                
        Arguments:
        data_list: [ndarray] array containing all station names and locations
        setup: [object] ini file parameters
        sounding: [ndarray] atmospheric profile

        Returns:
        allTimes: [ndarray] a ndarray that stores the arrival times for all stations over every perturbation
        [perturbation, station, 0 - ballistic/ 1 - fragmentation, frag number (0 for ballistic)]
        """

        #All perturbation happens here
        allTimes = [0]*setup.perturb_times

        # Ballistic Prediction
        ref_pos = position(setup.lat_centre, setup.lon_centre, 0)

        no_of_frags = len(setup.fragmentation_point)

        # array of frags and ballistic arrivals have to be the same size. So, minimum can be 1
        if no_of_frags == 0:
            no_of_frags = 1

        # Initialize variables
        b_time = 0

        consts = Constants()

        # For temporal perturbations, fetch the soudning data for the hour before and after the event
        if setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(setup, consts, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(setup, consts, time=-1)

        else:
            sounding_u = []
            sounding_l = []

        d_time = 2*(setup.perturb_times*len(stn_list)*no_of_frags)
        count = 0

        #number of perturbations
        for ptb_n in range(setup.perturb_times):

            if ptb_n > 0:
                
                if setup.debug:
                    print("STATUS: Perturbation {:}".format(ptb_n))

                # generate a perturbed sounding profile
                sounding_p = perturb(setup, sounding, setup.perturb_method, \
                    sounding_u=sounding_u, sounding_l=sounding_l, \
                    spread_file=setup.perturbation_spread_file, lat=setup.lat_centre, lon=setup.lon_centre)
            else:

                # if not using perturbations on this current step, then return the original sounding profile
                sounding_p = sounding

            # Initialize station times array
            stnTimes = [0]*len(stn_list)

            #number of stations
            for n, stn in enumerate(stn_list):

                # For ballistic arrivals
                if setup.show_ballistic_waveform:
                    bTimes = [0]*no_of_frags
                    for i in range(no_of_frags):
                        count += 1
                        sys.stdout.write("\rCalculating all times: {:5.2f} % ".format(count/d_time * 100))
                        sys.stdout.flush()
                        time.sleep(0.001)
                        #need filler values to make this a numpy array with fragmentation
                        if i == 0:
                            
                            stn.position.pos_loc(ref_pos)
                            setup.traj_f.pos_loc(ref_pos)

                            # Time to travel from trajectory to station
                            b_time = timeOfArrival([stn.position.x, stn.position.y, stn.position.z], setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, \
                                                        np.radians(setup.azim), np.radians(setup.zangle), setup, sounding=sounding_p, travel=False, fast=False, ref_loc=[ref_pos.lat, ref_pos.lon, ref_pos.elev])# + setup.t 


                            bTimes[i] = b_time
                        else:
                            bTimes[i] = np.nan
                else:
                    bTimes = [np.nan]*no_of_frags

                # Fragmentation Prediction
                f_time = np.array([0]*no_of_frags)

                # If manual fragmentation search is on
                if setup.show_fragmentation_waveform:
                    fTimes = [0]*no_of_frags
                    for i, line in enumerate(setup.fragmentation_point):
                        count += 1
                        sys.stdout.write("\rCalculating all times: {:5.2f} % ".format(count/d_time * 100))

                        # location of supracenter
                        supra = position(float(line[0]), float(line[1]), float(line[2]))
                        
                        # convert to local coordinates based off of the ref_pos
                        supra.pos_loc(ref_pos)

                        # convert station coordinates to local coordinates based on the ref_pos
                        stn.position.pos_loc(ref_pos)

                        # Cut down atmospheric profile to the correct heights, and interp
                        zProfile, _ = getWeather(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), setup.weather_type, \
                                [ref_pos.lat, ref_pos.lon, ref_pos.elev], copy.copy(sounding_p), convert=True)

                        # Travel time of the fragmentation wave
                        f_time, _, _ = cyscan(np.array([supra.x, supra.y, supra.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True, \
                            n_theta=setup.n_theta, n_phi=setup.n_phi, precision=setup.angle_precision, tol=setup.angle_error_tol)

                        fTimes[i] = f_time + line[3]

                else:

                    # Repack all arrays into allTimes array
                    fTimes = [np.nan]*no_of_frags

                stnTimes[n] = ([np.array(bTimes), np.array(fTimes)])

            allTimes[ptb_n] = np.array(stnTimes)

        allTimes = np.array(allTimes)

        # Save as .npy file to be reused in SeismicTrajectory and Supracenter
        np.save(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times'), allTimes)
        print("All Times File saved as {:}".format(os.path.join(setup.working_directory, setup.fireball_name, 'all_pick_times.npy')))

        return allTimes


    def waveformPicker(self, dir_path, setup, sounding, data_list, waveform_window=600, \
        difference_filter_all=False, stn_list=[]):
        """

        Arguments:
            data_list: [list]

        Keyword arguments:
            waveform_window: [int] Number of seconds for the wavefrom window.
            difference_filter_all: [bool] If True, the Kalenda et al. (2014) difference filter will be applied
                on the data plotted in the overview plot of all waveforms.
        """
        self.setup = setup
        self.dir_path = dir_path

        self.v_sound = setup.v_sound
        self.t0 = setup.t0

        self.stn_list = stn_list

        # Filter out all stations for which the mseed file does not exist
        filtered_stn_list = []

        names = []
        lats = []
        lons = []

        for stn in self.stn_list:
            
            mseed_file = stn.file_name
            mseed_file_path = os.path.join(self.dir_path, mseed_file)

            if os.path.isfile(mseed_file_path):
                filtered_stn_list.append(stn)

            else:
                print('mseed file does not exist:', mseed_file_path)

        self.stn_list = filtered_stn_list


        self.lat_centre = setup.lat_centre
        self.lon_centre = setup.lon_centre

        self.waveform_window = waveform_window


        self.current_station = 0
        self.current_wavefrom_raw = None
        self.current_wavefrom_delta = None
        self.current_waveform_processed = None

        # List of picks
        self.pick_list = []

        self.pick_group = 0

        # Define a list of colors for groups
        self.pick_group_colors = ['r', 'g', 'm', 'k', 'y']

        # Current station map handle
        self.current_station_scat = None

        # Station waveform marker handle
        self.current_station_all_markers = None

        # Picks on all waveform plot handle
        self.all_waves_picks_handle = None

        # Handle for pick text
        self.pick_text_handle = None
        self.pick_markers_handles = []

        # handle for pick marker on the wavefrom
        self.pick_wavefrom_handle = None


        # Default bandpass values
        self.bandpass_low_default = 2.0
        self.bandpass_high_default = 8.0

        # Flag indicating whether CTRL is pressed or not
        self.ctrl_pressed = False


        ### Sort stations by distance from source ###

        # Calculate distances of station from source
        self.source_dists = []


        for stn in self.stn_list:

            stat_name, stat_lat, stat_lon = stn.code, stn.position.lat, stn.position.lon

            names.append(stat_name)
            lats.append(stat_lat)
            lons.append(stat_lon)

            # Calculate the distance in kilometers
            dist = greatCircleDistance(np.radians(setup.lat_centre), np.radians(setup.lon_centre), \
                np.radians(stat_lat), np.radians(stat_lon))

            self.source_dists.append(dist)

        # Get sorted arguments
        dist_sorted_args = np.argsort(self.source_dists)

        # Sort the stations by distance
        self.stn_list = [self.stn_list[i] for i in dist_sorted_args]
        self.source_dists = [self.source_dists[i] for i in dist_sorted_args]
        #############################################
        
        # Init the plot framework
        self.initPlot(setup, sounding)

        # Extract the list of station locations
        self.lat_list = [stn.position.lat_r for stn in stn_list]
        self.lon_list = [stn.position.lon_r for stn in stn_list]

        plt.style.use('dark_background')
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(8, 5)
        self.map_ax = fig.add_subplot(1, 1, 1)

        # Init ground map
        self.m = GroundMap(self.lat_list, self.lon_list, ax=self.map_ax, color_scheme='light')

        for stn in self.stn_list:

            # Plot stations
            if stn.code in setup.high_f:
                self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='g', s=2)
            elif stn.code in setup.high_b:
                self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='b', s=2)
            else:
                self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='k', s=2)

        # Manual Supracenter search
        if setup.show_fragmentation_waveform:
            
            # Fragmentation plot
            for i, line in enumerate(setup.fragmentation_point):
                self.m.scatter([np.radians(float(line[0]))], [np.radians(float(line[1]))], c=self.pick_group_colors[(i+1)%4], marker='x')

        # Extract coordinates of the reference station
        ref_pos = position(setup.lat_centre, setup.lon_centre, 0)

        # Plot source location
        self.m.scatter([np.radians(setup.lat_centre)], [np.radians(setup.lon_centre)], marker='*', c='yellow')

        # Manual trajectory search
        if setup.show_ballistic_waveform:

            # Plot the trajectory with the bottom point known
            self.m.plot([setup.traj_i.lat_r, setup.traj_f.lat_r], [setup.traj_i.lon_r, setup.traj_f.lon_r], c='b')
            # Plot intersection with the ground
            self.m.scatter(setup.traj_f.lat_r, setup.traj_f.lon_r, s=10, marker='x', c='b')

            ### CONTOUR ###

            # Get the limits of the plot
            x_min = setup.lat_f - 100000*setup.deg_radius
            x_max = setup.lat_f + 100000*setup.deg_radius
            y_min = setup.lon_f - 100000*setup.deg_radius
            y_max = setup.lon_f + 100000*setup.deg_radius

            img_dim = int(setup.contour_res)
            x_data = np.linspace(x_min, x_max, img_dim)
            y_data = np.linspace(y_min, y_max, img_dim)
            xx, yy = np.meshgrid(x_data, y_data)


            # # Make an array of all plane coordinates
            plane_coordinates = np.c_[xx.ravel(), yy.ravel(), np.zeros_like(xx.ravel())]

            times_of_arrival = np.zeros_like(xx.ravel())

            az = np.radians(setup.azim)
            ze = np.radians(setup.zangle)

            # vector of the trajectory of the fireball
            traj_vect = np.array([np.sin(az)*np.sin(ze), np.cos(az)*np.sin(ze), -np.cos(ze)])

            for i, plane_coords in enumerate(plane_coordinates):

                # Print out percentage complete
                if (i + 1) % 10 == 0:
                    sys.stdout.write("\rDrawing Contour: {:.2f} %".format(100*(i + 1)/img_dim**2))
                    sys.stdout.flush()
                    time.sleep(0.001)

                setup.traj_f.pos_loc(ref_pos)
                # Point on the trajectory where the sound wave that will hit the plane_coord originated from

                p = waveReleasePoint(plane_coords, setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, az, \
                                          ze, setup.v_sound)

                # # vector between the wave release point and the plane coordinate
                d_vect = plane_coords - p

                # Since the arrivals are always perpendicular to the fireball trajectory, only take arrivals where the dot product
                # of the vectors are small.

                if abs(np.dot(d_vect/1000, traj_vect)) < setup.dot_tol:

                    # time of arrival from the trajectory
                    ti = timeOfArrival(plane_coords, setup.traj_f.x, setup.traj_f.y, setup.t0, 1000*setup.v, \
                                       az, ze, setup, sounding=sounding, ref_loc=[ref_pos.lat, ref_pos.lon, 0], travel=True, fast=True)# - setup.t + setup.t0

                # escape value for if sound never reaches the plane_coord
                else:
                   ti = np.nan

                times_of_arrival[i] = ti + setup.t0

            print('')

            # if sound never reaches the plane_coord, set to maximum value of the contour
            max_time = np.nanmax(times_of_arrival)
            for i in range(len(times_of_arrival)):
                if np.isnan(times_of_arrival[i]):
                    times_of_arrival[i] = max_time

            times_of_arrival = times_of_arrival.reshape(img_dim, img_dim)

            # Determine range and number of contour levels, so they are always centred around 0
            toa_abs_max = np.max([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
            #  toa_abs_min = np.min([np.abs(np.min(times_of_arrival)), np.max(times_of_arrival)])
            levels = np.linspace(0, toa_abs_max, 25)

            # Convert contour local coordinated to geo coordinates
            lat_cont = []
            lon_cont = []

            for x_cont, y_cont in zip(xx.ravel(), yy.ravel()):
                
                lat_c, lon_c, _ = loc2Geo(ref_pos.lat, ref_pos.lon, ref_pos.elev, np.array([x_cont, y_cont, 0]))

                lat_cont.append(lat_c)
                lon_cont.append(lon_c)

            lat_cont = np.array(lat_cont).reshape(img_dim, img_dim)
            lon_cont = np.array(lon_cont).reshape(img_dim, img_dim)

            # Plot the time of arrival contours
            toa_conture = self.m.m.contourf(lon_cont, lat_cont, times_of_arrival, levels, zorder=3, \
                latlon=True, cmap='viridis_r', alpha=0.5)

            # # Add a color bar which maps values to colors
            self.m.m.colorbar(toa_conture, label='Time of arrival (s)')

        if setup.arrival_times_file != '':
            try:
                self.arrTimes = np.load(setup.arrival_times_file)
                print("Reading in arrival times file...")
            except:
                print("WARNING: Unable to load allTimes_file {:} . Please check that file exists".format(setup.arrival_times_file))
                self.arrTimes = self.calcAllTimes(self.stn_list, setup, sounding)
        else:  
            # Calculate all arrival times
            self.arrTimes = self.calcAllTimes(self.stn_list, setup, sounding)
        

        self.make_picks_top_graphs.removeWidget(self.make_picks_map_graph_canvas)
        self.make_picks_map_graph_canvas = FigureCanvas(Figure(figsize=(1, 1)))
        self.make_picks_map_graph_canvas = FigureCanvas(fig)
        self.make_picks_map_graph_canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.make_picks_top_graphs.addWidget(self.make_picks_map_graph_canvas)    
        self.make_picks_map_graph_canvas.draw()
        SolutionGUI.update(self)

        self.updatePlot(setup)

    def makeValueChange(self, obj, slider):
        
        if obj == self.low_bandpass_label:
            obj.setText('Low: {:8.2f} Hz'.format(slider.value()*self.bandpass_scale))
        elif obj == self.high_bandpass_label:
            obj.setText('High: {:8.2f} Hz'.format(slider.value()*self.bandpass_scale))
        else:
            self.errorMessage('Bad atm slider pass in makeValueChange', 2)

        self.updatePlot()
        #self.atmPlotProfile(self.atm_lat_slide.value()*self.slider_scale, self.atm_lon_slide.value()*self.slider_scale, self.var_typ)
    
    def chooseFilter(self, obj):

        if obj.currentText() == 'Bandpass Filter':
            self.filterBandpass()
        elif obj.currentText() == 'Spectrogram of Raw Data':
            self.showSpectrogram()
        elif obj.currentText() == 'Difference Filter':
            self.filterConvolution()
        else:
            pass

    def navStats(self):

        self.current_station = self.make_picks_station_choice.currentIndex()
        self.updatePlot()


    def initPlot(self, setup, sounding):
        """ Initializes the plot framework. """
        
               ### Init the basic grid ###


        # Register a mouse press event on the waveform axis
        # plt.gca().figure.canvas.mpl_connect('button_press_event', self.onWaveMousePress)


        # Register window resize
        #plt.gca().figure.canvas.mpl_connect('resize_event', self.onResize)

        self.make_picks_station_choice.clear()

        self.prev_stat.clicked.connect(self.decrementStation)
        self.next_stat.clicked.connect(self.incrementStation)

        self.filter_combo_box.addItem('Raw Data')
        self.filter_combo_box.addItem('Bandpass Filter')
        self.filter_combo_box.addItem('Spectrogram of Raw Data')
        self.filter_combo_box.addItem('Difference Filter')

        self.filter_combo_box.currentTextChanged.connect(partial(self.chooseFilter, self.filter_combo_box))

        self.export_to_csv.clicked.connect(self.exportCSV)

        for stn in self.stn_list:
            self.make_picks_station_choice.addItem("{:}-{:}".format(stn.network, stn.code))

        self.make_picks_station_choice.currentTextChanged.connect(self.navStats)

        plt.style.use('dark_background')
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(8, 5)
        self.station_ax = fig.add_subplot(1, 1, 1)
        
        self.make_picks_waveform_canvas.scene().sigMouseClicked.connect(self.mouseClicked)
        pg.QtGui.QApplication.processEvents()
        # Plot all waveforms
        plotAllWaveforms(self.dir_path, list(self.stn_list), setup, sounding, ax=self.station_ax)#, \
            #waveform_window=self.waveform_window, difference_filter_all=setup.difference_filter_all)

        self.make_picks_top_graphs.removeWidget(self.make_picks_station_graph_canvas)
        self.make_picks_station_graph_canvas = FigureCanvas(Figure(figsize=(1, 1)))
        self.make_picks_station_graph_canvas = FigureCanvas(fig)
        self.make_picks_station_graph_canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        self.make_picks_top_graphs.addWidget(self.make_picks_station_graph_canvas)    
        self.make_picks_station_graph_canvas.draw()
        SolutionGUI.update(self)

    def keyPressEvent(self, event):
        
        if event.key() == QtCore.Qt.Key_Control:
            self.ctrl_pressed = True

        elif event.key() == QtCore.Qt.Key_D:
            self.incrementStation()

        elif event.key() == QtCore.Qt.Key_A:
            self.decrementStation()

        elif event.key() == QtCore.Qt.Key_W:
            self.make_picks_waveform_canvas.clear()
            self.filterBandpass(event=event)

        elif event.key() == QtCore.Qt.Key_S:
            self.showSpectrogram(event=event)

        elif event.key() == QtCore.Qt.Key_C:
            self.make_picks_waveform_canvas.clear()
            self.filterConvolution(event=event)

        elif event.key() == QtCore.Qt.Key_Plus:
            
            # Increment the pick group
            self.pick_group += 1

            self.updatePlot()


        elif event.key() == QtCore.Qt.Key_Minus:
            # Decrement the pick group

            if self.pick_group > 0:
                self.pick_group -= 1

            self.updatePlot()


    def keyReleaseEvent(self, event):

        if event.key() == QtCore.Qt.Key_Control:
            self.ctrl_pressed = False



    def onResize(self, event):

        # Perform tight layout when window is resized
        plt.tight_layout()




    def addPick(self, pick_group, station_no, pick_time):
        """ Adds the pick to the list of picks. """

        self.pick_list.append([pick_group, station_no, pick_time])

        self.updatePickList()


    def removePick(self, station_no, pick_time_remove):
        """ Removes the pick from the list of picks with the closest time. """


        if len(self.pick_list):

            closest_pick_indx = None
            min_time_diff = np.inf

            # Go though all stations and find the pick with closest to the given time
            for i, entry in enumerate(self.pick_list):

                pick_grp, stat_no, pick_time = entry

                # Check if current station
                if stat_no == self.current_station:

                    time_diff = abs(pick_time - pick_time_remove)

                    # Store the minimum time difference
                    if time_diff < min_time_diff:

                        min_time_diff = time_diff
                        closest_pick_indx = i



            if closest_pick_indx is not None:
                
                # Remove pick on the given station closest the given time
                self.pick_list.pop(closest_pick_indx)


                self.updatePickList()




    def updatePickList(self):
        """ Updates the list of picks on the screen for the given station. """

        stations_with_picks = []
        all_pick_times = []

        for entry in self.pick_list:

            pick_grp, station_no, pick_time = entry

            stations_with_picks.append(station_no)
            all_pick_times.append(pick_time)


        # Remove old picks on all wavefrom plot
        if self.all_waves_picks_handle is not None:
            self.all_waves_picks_handle.remove()


        # Get distances of of pick stations
        dists_with_picks = [self.source_dists[stat_no] for stat_no in stations_with_picks]

        # Mark picks on the all waveform plot
        self.all_waves_picks_handle = self.ax_all_waves.scatter(dists_with_picks, all_pick_times, \
            marker='*', s=50, c='r')

        self.updatePlot(draw_waveform=False)


    def updatePickTextAndWaveMarker(self):
        """ Updates the list of picks on the screen. """


        current_station_groups = []
        current_station_picks = []

        for entry in self.pick_list:

            pick_grp, station_no, pick_time = entry

            # Take picks taken on the current station
            if station_no == self.current_station:
                current_station_groups.append(pick_grp)
                current_station_picks.append(pick_time)

        # Remove old pick text
        if self.pick_text_handle is not None:
            self.pick_text_handle.remove()

        # Generate the pick string
        pick_txt_str  = 'Change group: +/-\n'
        pick_txt_str += 'Add/remove pick: CTRL + left/right click \n'
        pick_txt_str += '\n'
        pick_txt_str += 'Current group: {:5d}\n\n'.format(self.pick_group)
        pick_txt_str += 'Picks: Group, Time\n'
        pick_txt_str += '------\n'
        pick_txt_str += "\n".join(["{:5d},   {:.2f}".format(gr, pt) for gr, pt in zip(current_station_groups, \
            current_station_picks)])

        # Print picks on screen
        #self.pick_text_handle = self.picks_ax.text(0, 1, pick_txt_str, va='top', fontsize=7)


        # Remove old pick markers
        for handle in self.pick_markers_handles:
            try:
                handle.remove()
            except:
                pass


        self.pick_markers_handles = []


        if len(current_station_picks) > 0:

            # Get a list of colors per groups
            color_list = [self.pick_group_colors[grp%len(self.pick_group_colors)] \
                for grp in current_station_groups]


            # Get the Y coordinate of the pick in the waveform plot
            if self.current_waveform_processed is not None:
                pick_y_list = []
                for pick_time in current_station_picks:
                    
                    pick_indx = np.abs(self.current_waveform_time - pick_time).argmin()
                    pick_y = self.current_waveform_processed[pick_indx]
                    pick_y_list.append(pick_y)

            else:
                pick_y_list = [0]*len(current_station_picks)

            # Set pick marker on the current wavefrom
            scat_handle = self.ax_wave.scatter(current_station_picks, pick_y_list, marker='*', \
                c=color_list, s=50)

            self.pick_markers_handles.append(scat_handle)

            # Plot group numbers above picks
            #self.pick_wavefrom_text_handles = []
            for c, grp, pt in zip(color_list, current_station_groups, current_station_picks):
                txt_handle = self.ax_wave.text(pt, 0, str(grp), color=c, ha='center', va='bottom')

                self.pick_markers_handles.append(txt_handle)


    def onWaveMousePress(self, event):

        # Check if the mouse was pressed within the waveform axis
        if event.inaxes == self.ax_wave:

            # Check if CTRL is pressed
            if self.ctrl_pressed:

                pick_time = event.xdata

                # Check if left button was pressed
                if event.button == 1:

                    # Extract network and station code
                    net, station_code = self.stn_list[self.current_station].network, self.stn_list[self.current_station].code

                    print('Adding pick on station {:s} at {:.2f}'.format(net + ": " + station_code, \
                        pick_time))

                    self.addPick(self.pick_group, self.current_station, pick_time)


                # Check if right button was pressed
                elif event.button == 3:
                    print('Removing pick...')

                    self.removePick(self.current_station, pick_time)





    def incrementStation(self, event=None):
        """ Increments the current station index. """

        self.make_picks_waveform_canvas.clear()

        self.current_station += 1

        if self.current_station >= len(self.stn_list):
            self.current_station = 0

        while self.checkExists() == False:
            self.current_station += 1
            if self.current_station >= len(self.stn_list):
                self.current_station = 0

        self.updatePlot()


    def decrementStation(self, event=None):
        """ Decrements the current station index. """

        self.make_picks_waveform_canvas.clear()

        self.current_station -= 1

        if self.current_station < 0:
            self.current_station = len(self.stn_list) - 1

        while self.checkExists() == False:
            self.current_station -= 1
            if self.current_station < 0:
                self.current_station = len(self.stn_list) - 1

        self.updatePlot()



    def markCurrentStation(self):
        """ Mark the position of the current station on the map. """


        # Extract current station
        stn = self.stn_list[self.current_station]

        if self.current_station_scat is None:

            # Mark the current station on the map
            self.current_station_scat = self.m.scatter([stn.position.lat_r], [stn.position.lon_r], s=20, \
                edgecolors='r', facecolors='none')

        else:

            # Calculate map coordinates
            stat_x, stat_y = self.m.m(stn.position.lon, stn.position.lat)

            # Set the new position
            self.current_station_scat.set_offsets([stat_x, stat_y])


    def checkExists(self):
        """
        Checks if the current waveform is readable
        """

        # Extract current station
        stn = self.stn_list[self.current_station]

        # Get the miniSEED file path
        mseed_file_path = os.path.join(self.dir_path, stn.file_name)

        try:
            
            if os.path.isfile(mseed_file_path):
                pass
            else:
                print('File {:s} does not exist!'.format(mseed_file_path))
                return False

        except TypeError as e:
            
            print('Opening file {:s} failed with error: {:s}'.format(mseed_file_path, str(e)))
            return False

        try:
            obspy.read(mseed_file_path)

        except TypeError:
            print('mseed file could not be read:', mseed_file_path)
            return False

        return True
    
    def mouseClicked(self, evt):

        mousePoint = self.make_picks_waveform_canvas.vb.mapSceneToView(evt.pos())

        self.make_picks_waveform_canvas.scatterPlot(x=[mousePoint.x()], y=[mousePoint.y()], pen='r', update=True)

        
        print((mousePoint.x(), mousePoint.y()))


    def drawWaveform(self, waveform_data=None):
        """ Draws the current waveform from the current station in the wavefrom window. Custom wavefrom 
            can be given an drawn, which is used when bandpass filtering is performed. 

        """
        setup = self.setup

        # Clear waveform axis
        self.make_picks_waveform_canvas.clear()

        # Extract current station
        stn = self.stn_list[self.current_station]

        # Get the miniSEED file path
        mseed_file_path = os.path.join(self.dir_path, stn.file_name)

        # Try reading the mseed file, if it doesn't work, skip to the next frame
        try:
            mseed = obspy.read(mseed_file_path)

        except TypeError:
            print('mseed file could not be read:', mseed_file_path)
            #self.incrementStation()
            return None

        # Unpact miniSEED data
        delta = mseed[0].stats.delta
        start_time = mseed[0].stats.starttime
        end_time = mseed[0].stats.endtime
        

        # Check if the waveform data is already given or not
        if waveform_data is None:
            waveform_data = mseed[0].data

            # Store raw data for bookkeeping on first open
            self.current_wavefrom_raw = waveform_data


        # Convert the beginning and the end time to datetime objects
        start_datetime = start_time.datetime
        end_datetime = end_time.datetime

        self.current_wavefrom_delta = delta
        self.current_waveform_time = np.arange(0, (end_datetime - start_datetime).total_seconds() + delta, \
            delta)


        # ### BANDPASS FILTERING ###

        # # Init the butterworth bandpass filter
        # butter_b, butter_a = butterworthBandpassFilter(self.bandpass_low_default, \
        #     self.bandpass_high_default, 1.0/delta, order=6)

        # # Filter the data
        # waveform_data = scipy.signal.filtfilt(butter_b, butter_a, waveform_data)

        # ##########################


        # Construct time array, 0 is at start_datetime
        time_data = np.copy(self.current_waveform_time)

        # Cut the waveform data length to match the time data
        waveform_data = waveform_data[:len(time_data)]
        time_data = time_data[:len(waveform_data)]


        # Trim the ends to avoid issues with differenced data
        # time_data = time_data[1:-1]
        # waveform_data = waveform_data[1:-1]

        # Store currently plotted waveform
        self.current_waveform_processed = waveform_data

        # Calculate the time of arrival assuming constant propagation with the given speed of sound
        t_arrival = self.source_dists[self.current_station]/(self.v_sound/1000) + self.t0
        #t_arrival = 0

        # Calculate the limits of the plot to be within the given window limit
        time_win_min = t_arrival - self.waveform_window/2
        time_win_max = t_arrival + self.waveform_window/2

        # Plot the wavefrom
        #self.wave_ax.plot(time_data, waveform_data, color='k', linewidth=0.2, zorder=3)
        self.make_picks_waveform_canvas.plot(x=time_data, y=waveform_data, pen='w')
        #self.make_picks_waveform_canvas.setXRange(t_arrival, t_arrival, padding=150)
        # proxy = pg.SignalProxy(self.make_picks_waveform_canvas.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)
        #self.atm_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)


        SolutionGUI.update(self)
        # Initialize variables
        b_time = 0

        # print('####################')
        # print("Current Station: {:}".format(stn.name))
        # print("Channel: {:}".format(stn.channel))
        # # If manual ballistic search is on
        # if setup.show_ballistic_waveform:

        #     # Plot Ballistic Prediction
        #     b_time = self.arrTimes[0, self.current_station, 0, 0]
            
        #     # check if nan
        #     if b_time == b_time:
        #         self.wave_ax.plot([b_time]*2, [np.min(waveform_data), np.max(waveform_data)], c='b', label='Ballistic', zorder=3)
                
        #         print("Ballistic Arrival: {:.3f} s".format(b_time))
        #     else:
        #         print("No Ballistic Arrival")

        #     for i in range(setup.perturb_times):
        #         if i >= 1:
        #             try:
        #                 self.wave_ax.plot([self.arrTimes[i, self.current_station, 0, 0]]*2, \
        #                  [np.min(waveform_data), np.max(waveform_data)], alpha=0.3, c='b', zorder=3)
        #             except:
        #                 pass
        # # Fragmentation Prediction

        # # If manual fragmentation search is on
        # if setup.show_fragmentation_waveform:

        #     for i, line in enumerate(setup.fragmentation_point):

        #         f_time = self.arrTimes[0, self.current_station, 1, i]
        #     #     # check if nan
        #         if f_time == f_time:
        #             # Plot Fragmentation Prediction
        #             self.wave_ax.plot([f_time]*2, [np.min(waveform_data), np.max(waveform_data)], c=self.pick_group_colors[(i+1)%4], label='Fragmentation', zorder=3)
                    
        #             if len(setup.fragmentation_point) > 1:
        #                 #self.ax_wave.text(f_time, np.min(waveform_data), 'Frag{:}'.format(i+1))
        #                 self.wave_ax.text(f_time, np.min(waveform_data) + int(i)/(len(setup.fragmentation_point))*(np.max(waveform_data) - np.min(waveform_data)), '{:.1f} km'.format(line[2]/1000))
                
        #             print('Fragmentation {:} Arrival: {:.3f} s'.format(i+1, f_time))

        #         else:
        #             print('No Fragmentation {:} Arrival'.format(i+1))

        #         for j in range(setup.perturb_times):
        #             if j >= 1:
        #                 try:
        #                     self.wave_ax.plot([self.arrTimes[j, self.current_station, 1, i]]*2, [np.min(waveform_data),\
        #                          np.max(waveform_data)], alpha=0.3,\
        #                          c=self.pick_group_colors[(i+1)%4], zorder=3)
        #                 except:
        #                     pass

        # # Set the time limits to be within the given window
        # self.wave_ax.set_xlim(time_win_min, time_win_max)

        # self.wave_ax.grid(color='#ADD8E6', linestyle='dashed', linewidth=0.5, alpha=0.5)

        # #self.ax_wave.legend()

        # # Add text with station label
        # if stn.code in setup.high_f:
        #     self.wave_ax.text(time_win_min, np.max(waveform_data), stn.network + ": " + stn.code \
        #         + "(" + stn.channel + ")" +", {:d} km".format(int(self.source_dists[self.current_station])) , va='top', ha='left', color='g')

        # elif stn.code in setup.high_b:
        #     self.wave_ax.text(time_win_min, np.max(waveform_data), stn.network + ": " + stn.code \
        #         + "(" + stn.channel + ")" + ", {:d} km".format(int(self.source_dists[self.current_station])) , va='top', ha='left', color='b')

        # else:
        #     self.wave_ax.text(time_win_min, np.max(waveform_data), stn.network + ": " + stn.code \
        #         + "(" + stn.channel + ")" + ", {:d} km".format(int(self.source_dists[self.current_station])) , va='top', ha='left', color='k')

        # # self.make_picks_top_graphs.removeWidget(self.make_picks_station_graph_canvas)
        # # self.make_picks_station_graph_canvas = FigureCanvas(Figure(figsize=(3, 3)))
        # # self.make_picks_station_graph_canvas = FigureCanvas(fig)
        # # self.make_picks_station_graph_canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        # # self.make_picks_top_graphs.addWidget(self.make_picks_station_graph_canvas)    
        # # self.make_picks_station_graph_canvas.draw()
        # # SolutionGUI.update(self)

    def markStationWaveform(self):
        """ Mark the currently shown waveform in the plot of all waveform. """
        
        if self.current_station_all_markers is not None:
            for marker in self.current_station_all_markers:
                marker.remove()


        # Calculate the position
        dist = self.source_dists[self.current_station]

        # Calcualte the time of arrival
        t_arrival = self.source_dists[self.current_station]/(self.v_sound/1000) + self.t0

        # Plot the marker
        marker1 = self.station_ax.scatter(dist, t_arrival - 500, marker='^', s=100, linewidths=3, c='w', 
            alpha=1, zorder=3)

        marker2 = self.station_ax.scatter(dist, t_arrival + 500, marker='v', s=100, linewidths=3, c='w', 
            alpha=1, zorder=3)


        self.current_station_all_markers = [marker1, marker2]


    def showSpectrogram(self, event=None):
        """ Show the spectrogram of the waveform in the current window. """


        # Get time limits of the shown waveform
        #x_min, x_max = self.wave_ax.get_xlim()

        # Extract the time and waveform
        #crop_window = (self.current_waveform_time >= x_min) & (self.current_waveform_time <= x_max)
        wave_arr = self.current_wavefrom_raw#[crop_window]


        ### Show the spectrogram ###
        
        fig = plt.figure()
        ax_spec = fig.add_subplot(111)

        ax_spec.specgram(wave_arr, Fs=1.0/self.current_wavefrom_delta, cmap=plt.cm.inferno)

        ax_spec.set_xlabel('Time (s)')
        ax_spec.set_ylabel('Frequency (Hz)')

        fig.show()

        ###


    def filterBandpass(self, event=None):
        """ Run bandpass filtering using values set on sliders. """

        # Get bandpass filter values
        bandpass_low = self.low_bandpass_slider.value()*self.bandpass_scale
        bandpass_high = self.high_bandpass_slider.value()*self.bandpass_scale


        # Limit the high frequency to be lower than the Nyquist frequency
        max_freq = (1.0/self.current_wavefrom_delta)/2

        if bandpass_high > max_freq:
            bandpass_high = max_freq - 0.1

            self.high_bandpass_slider.setValue(bandpass_high*self.bandpass_scale)
        

        # Init the butterworth bandpass filter
        butter_b, butter_a = butterworthBandpassFilter(bandpass_low, bandpass_high, \
            1.0/self.current_wavefrom_delta, order=6)

        # Filter the data
        waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(self.current_wavefrom_raw))


        # Plot the updated waveform
        self.drawWaveform(waveform_data)


    def filterConvolution(self, event=None):
        """ Apply the convolution filter on data as suggested in Kalenda et al. (2014). """

        waveform_data = convolutionDifferenceFilter(self.current_wavefrom_raw)

        self.drawWaveform(waveform_data)


    def updatePlot(self, draw_waveform=True):
        """ Update the plot after changes. """

        self.make_picks_waveform_canvas.clear()

        # Mark the position of the current station on the map

        self.make_picks_station_choice.setCurrentIndex(self.current_station)

        # Plot the wavefrom from the current station
        if draw_waveform:
            self.drawWaveform()

        # Set an arrow pointing to the current station on the waveform
        self.markStationWaveform()
        self.markCurrentStation()

        # Update the pick list text and plot marker on the waveform
        self.updatePickTextAndWaveMarker()

        # Reset bandpass filter values to default
        # self.make_picks.set_val(self.bandpass_low_default)
        # self.bandpass_high_slider.set_val(self.bandpass_high_default)

        SolutionGUI.update(self)

    def exportCSV(self, event):
        """ Save picks to a CSV file. """

        # Open the output CSV
        with open(os.path.join(self.dir_path, OUTPUT_CSV), 'w') as f:

            # Write the header
            f.write('Pick group, Network, Code, Lat, Lon, Elev, Pick JD, Pick time, station_number \n')

            # Go through all picks
            for entry in self.pick_list:

                # Unpack pick data
                pick_group, station_no, pick_time = entry

                # # Extract current station
                stn = self.stn_list[station_no]

                # # Unpack the station entry
                # net, station_code, stat_lat, stat_lon, stat_elev, station_name, channel, mseed_file = stat_entry

                # Get the miniSEED file path
                mseed_file_path = os.path.join(self.dir_path, stn.file_name)

                try:
                    
                    if os.path.isfile(mseed_file_path):
                        
                        # Read the miniSEED file
                        mseed = obspy.read(mseed_file_path)

                    else:
                        print('File {:s} does not exist!'.format(mseed_file_path))
                        continue

                except TypeError as e:

                    print('Opening file {:s} failed with error: {:s}'.format(mseed_file_path, e))
                    continue

                # Find datetime of the beginning of the file
                start_datetime = mseed[0].stats.starttime.datetime

                # Calculate Julian date of the pick time
                pick_jd = datetime2JD(start_datetime + datetime.timedelta(seconds=pick_time))


                # Write the CSV entry
                f.write("{:d}, {:s}, {:s}, {:.6f}, {:.6f}, {:.2f}, {:.8f}, {:}, {:}\n".format(pick_group, stn.network, \
                    stn.code, stn.position.lat, stn.position.lon, stn.position.elev, pick_jd, pick_time, station_no))


        print('CSV written to:', OUTPUT_CSV)

    def addMakePicksWidgets(self):
        make_picks_master_tab = QWidget()
        make_picks_master = QVBoxLayout()
        make_picks_master_tab.setLayout(make_picks_master)

        self.make_picks_top_graphs = QHBoxLayout()
        self.make_picks_bottom_graphs = QVBoxLayout()
        make_picks_master.addLayout(self.make_picks_top_graphs)
        make_picks_master.addLayout(self.make_picks_bottom_graphs)

        self.make_picks_station_graph_canvas = FigureCanvas(Figure(figsize=(1, 1)))
        self.make_picks_map_graph_canvas = FigureCanvas(Figure(figsize=(1, 1)))
        self.make_picks_top_graphs.addWidget(self.make_picks_station_graph_canvas)
        self.make_picks_top_graphs.addWidget(self.make_picks_map_graph_canvas)


        self.make_picks_waveform_view = pg.GraphicsLayoutWidget()
        self.make_picks_waveform_canvas = self.make_picks_waveform_view.addPlot()
        self.make_picks_bottom_graphs.addWidget(self.make_picks_waveform_view)
        self.make_picks_waveform_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        make_picks_control_panel = QHBoxLayout()
        self.make_picks_bottom_graphs.addLayout(make_picks_control_panel)

        make_picks_station_group = QGroupBox("Station Navigation")
        make_picks_control_panel.addWidget(make_picks_station_group)

        station_group_layout = QGridLayout()
        make_picks_station_group.setLayout(station_group_layout)

        self.make_picks_station_choice = QComboBox()
        station_group_layout.addWidget(self.make_picks_station_choice, 0, 0, 1, 2)

        self.prev_stat = QPushButton('Prev')
        station_group_layout.addWidget(self.prev_stat, 1, 0, 2, 1)

        self.next_stat = QPushButton('Next')
        station_group_layout.addWidget(self.next_stat, 1, 1, 2, 1)

        launch = QPushButton('Load Station Data')
        station_group_layout.addWidget(launch, 3, 0, 1, 2)
        launch.clicked.connect(self.makePicks)

        make_picks_filter_group = QGroupBox("Waveform Filtering")
        make_picks_control_panel.addWidget(make_picks_filter_group)

        filter_group_layout = QGridLayout()
        make_picks_filter_group.setLayout(filter_group_layout)

        self.low_bandpass_label = QLabel('Low: 2 Hz')
        filter_group_layout.addWidget(self.low_bandpass_label, 0, 0, 1, 1)

        self.low_bandpass_slider = QSlider(Qt.Horizontal)
        filter_group_layout.addWidget(self.low_bandpass_slider, 0, 1, 1, 1)

        self.high_bandpass_label = QLabel("High: 8 Hz")
        filter_group_layout.addWidget(self.high_bandpass_label, 1, 0, 1, 1)

        self.high_bandpass_slider = QSlider(Qt.Horizontal)
        filter_group_layout.addWidget(self.high_bandpass_slider, 1, 1, 1, 1)

        self.low_bandpass_slider.setMinimum(0/self.bandpass_scale)
        self.low_bandpass_slider.setMaximum(5/self.bandpass_scale)
        self.low_bandpass_slider.setValue(2/self.bandpass_scale)
        self.low_bandpass_slider.setTickInterval(0.5)
        self.low_bandpass_slider.valueChanged.connect(partial(self.makeValueChange, self.low_bandpass_label, self.low_bandpass_slider))

        self.high_bandpass_slider.setMinimum(3/self.bandpass_scale)
        self.high_bandpass_slider.setMaximum(40/self.bandpass_scale)
        self.high_bandpass_slider.setValue(8/self.bandpass_scale)
        self.high_bandpass_slider.setTickInterval(0.5)
        self.high_bandpass_slider.valueChanged.connect(partial(self.makeValueChange, self.high_bandpass_label, self.high_bandpass_slider))


        self.filter_combo_box = QComboBox()
        filter_group_layout.addWidget(self.filter_combo_box, 2, 0, 1, 2)

        make_picks_picks_group = QGroupBox("Arrival Picks")
        make_picks_control_panel.addWidget(make_picks_picks_group)

        pick_group_layout = QGridLayout()
        make_picks_picks_group.setLayout(pick_group_layout)

        self.export_to_csv = QPushButton('Export to CSV')
        pick_group_layout.addWidget(self.export_to_csv)

        self.tab_widget.addTab(make_picks_master_tab, 'Make Picks')         

    def changeRows(self, obj, change):
        n_rows = obj.rowCount()
        if change == 1:
            obj.insertRow(n_rows)
        else:
            obj.removeRow(n_rows-1)

    def fileSearch(self, filters, obj):
        # obj - where the text goes
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setNameFilters(filters)
        dlg.selectNameFilter(filters[0])
        #filenames = QStringList()

        dlg.exec_()

        filename = dlg.selectedFiles()

        obj.setText(filename[0])

    def folderSearch(self, obj):
        # dlg = QFileDialog()
        # dlg.getExistingDirectory(self, "Select Directory")

        # obj - where the text goes
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.Directory)
        #dlg.setNameFilters()
        # dlg.selectNameFilter(filters[0])
        #filenames = QStringList()

        dlg.exec_()

        filename = dlg.selectedFiles()

        obj.setText(filename[0])

    def supSearch(self):
        setup = self.saveINI(write=False)

        setup.search = [float(self.sup_south_edits.text()),
                        float(self.sup_north_edits.text()),
                        float(self.sup_west_edits.text()),
                        float(self.sup_east_edits.text()),
                        float(self.sup_height_min_edits.text()),
                        float(self.sup_height_max_edits.text())]
        setup.min_time = float(self.sup_min_time_edits.text())
        setup.max_time = float(self.sup_max_time_edits.text())
        setup.start_datetime = self.sup_ref_time_edits.dateTime().toPyDateTime()
        setup.station_picks_file = self.sup_picks_file_edit.text()

        s_info, s_name, weights, ref_pos = convStationDat(setup.station_picks_file)
        ref_pos = position(ref_pos[0], ref_pos[1], ref_pos[2])

        consts = Constants()      
        dataset = parseWeather(setup, consts)
        #dataset = findECMWFSound(lat, lon, sounding)

        setup.ref_pos = ref_pos

        results = psoSearch(s_info, weights, s_name, setup, dataset, consts)

        n_stations = len(s_info)
        xstn = s_info[0:n_stations, 0:3]
            
        self.scatterPlot(setup, results, n_stations, xstn, s_name, dataset, manual=False)

        self.residPlot(results, s_name, xstn, setup.working_directory, n_stations, manual=False)

        # set row count
        self.sup_results_table.setRowCount(n_stations + 1)

        # set column count
        self.sup_results_table.setColumnCount(5)

        self.sup_results_table.setItem(0, 0, QTableWidgetItem("Station Name"))
        self.sup_results_table.setItem(0, 1, QTableWidgetItem("Latitude"))
        self.sup_results_table.setItem(0, 2, QTableWidgetItem("Longitude"))
        self.sup_results_table.setItem(0, 3, QTableWidgetItem("Elevation"))
        self.sup_results_table.setItem(0, 4, QTableWidgetItem("Residuals"))

        for i in range(n_stations):
            self.sup_results_table.setItem(i+1, 0, QTableWidgetItem(s_name[i]))
            self.sup_results_table.setItem(i+1, 1, QTableWidgetItem(str(xstn[i][0])))
            self.sup_results_table.setItem(i+1, 2, QTableWidgetItem(str(xstn[i][1])))
            self.sup_results_table.setItem(i+1, 3, QTableWidgetItem(str(xstn[i][2])))
            self.sup_results_table.setItem(i+1, 4, QTableWidgetItem(str(results.r[i])))

    def supraSearch(self):

        setup = self.saveINI(write=False)

        supra_pos = [float(self.lat_edit.text()),
                     float(self.lon_edit.text()),
                     float(self.elev_edit.text()),
                     float(self.time_edit.text())]
        picks_file = self.picks_file_edit.text()
        atm_file = self.atmospheric_file_edit.text()

        setup.sounding_file = atm_file

        s_info, s_name, weights, ref_pos = convStationDat(picks_file)
        ref_pos = position(ref_pos[0], ref_pos[1], ref_pos[2])

        consts = Constants()      
        dataset = parseWeather(setup, consts)
        #dataset = findECMWFSound(lat, lon, sounding)


        setup.manual_fragmentation_search = supra_pos
        setup.ref_pos = ref_pos
        setup.start_datetime = self.ref_edit.dateTime().toPyDateTime()

        results = psoSearch(s_info, weights, s_name, setup, dataset, consts)

        n_stations = len(s_info)
        xstn = s_info[0:n_stations, 0:3]
            
        self.scatterPlot(setup, results, n_stations, xstn, s_name, dataset)

        self.residPlot(results, s_name, xstn, setup.working_directory, n_stations)

        # set row count
        self.tableWidget.setRowCount(n_stations + 1)

        # set column count
        self.tableWidget.setColumnCount(5)

        self.tableWidget.setItem(0, 0, QTableWidgetItem("Station Name"))
        self.tableWidget.setItem(0, 1, QTableWidgetItem("Latitude"))
        self.tableWidget.setItem(0, 2, QTableWidgetItem("Longitude"))
        self.tableWidget.setItem(0, 3, QTableWidgetItem("Elevation"))
        self.tableWidget.setItem(0, 4, QTableWidgetItem("Residuals"))

        for i in range(n_stations):
            self.tableWidget.setItem(i+1, 0, QTableWidgetItem(s_name[i]))
            self.tableWidget.setItem(i+1, 1, QTableWidgetItem(str(xstn[i][0])))
            self.tableWidget.setItem(i+1, 2, QTableWidgetItem(str(xstn[i][1])))
            self.tableWidget.setItem(i+1, 3, QTableWidgetItem(str(xstn[i][2])))
            self.tableWidget.setItem(i+1, 4, QTableWidgetItem(str(results.r[i])))

    def comboSet(self, obj, text):
        
        index = obj.findText(str(text).lower(), Qt.MatchFixedString)
        if index >= 0:
            obj.setCurrentIndex(index)

    def toTable(self, obj, table):
        
        if len(table) > 0:
            X = len(table)
            Y = len(table[0])

            obj.setRowCount(X)
            obj.setColumnCount(Y)

            for x in range(X):
                for y in range(Y):
                    obj.setItem(x, y, QTableWidgetItem(str(table[x][y])))
        else:
            print("Warning: Table has no length")

    def fromTable(self, obj):

        X = obj.rowCount()
        Y = obj.columnCount()

        table = []
    
        for x in range(X):
            line = [None]*Y
            for y in range(Y):
                try:
                    line[y] = float(obj.item(x, y).text())
                except:
                    line[y] = obj.item(x, y).text()
            table.append(line)
        return table

    def toTableFromStn(self, obj, table):
        
        X = len(table)

        obj.setRowCount(X)

        for x in range(X):
            stn = table[x]
            obj.setItem(x, 0, QTableWidgetItem(str(stn.network)))
            obj.setItem(x, 1, QTableWidgetItem(str(stn.code)))
            obj.setItem(x, 2, QTableWidgetItem(str(stn.position.lat)))
            obj.setItem(x, 3, QTableWidgetItem(str(stn.position.lon)))
            obj.setItem(x, 4, QTableWidgetItem(str(stn.position.elev)))
            obj.setItem(x, 5, QTableWidgetItem(str(stn.channel)))
            obj.setItem(x, 6, QTableWidgetItem(str(stn.name)))
            obj.setItem(x, 7, QTableWidgetItem(str(stn.file_name)))

    def tryFloat(self, fr):

        try:
            return float(fr)
        except:
            return None

    def tryInt(self, fr):

        try:
            return int(fr)
        except:
            return None

    def saveINI(self, write=True):

        """if write == True - used for saving data to setup obj, but not writing
        """
        if write:
            dlg = QFileDialog.getSaveFileName(self, 'Save File')

        class Config:

            def __init__(self):
                pass

        setup = Config()

        setup.fireball_name = self.fireball_name_edits.text()
        setup.difference_filter_all = self.difference_filter_edits.currentText()
        setup.get_data = self.get_data_edits.currentText()
        setup.run_mode = self.run_mode_edits.currentText()
        setup.debug = self.debug_edits.currentText()

        setup.working_directory = self.working_directory_edits.text()
        setup.arrival_times_file = self.arrival_times_edits.text()
        setup.sounding_file = self.sounding_file_edits.text()
        setup.perturbation_spread_file = self.perturbation_file_edits.text()
        setup.station_picks_file = self.station_picks_file_edits.text()
        setup.replot_points_file = self.points_name_edits.text()

        setup.lat_centre = self.tryFloat(self.lat_centre_edits.text())
        setup.lon_centre = self.tryFloat(self.lon_centre_edits.text())
        setup.deg_radius = self.tryFloat(self.deg_radius_edits.text())

        setup.start_datetime = self.start_datetime_edits.dateTime().toPyDateTime()
        setup.end_datetime = self.end_datetime_edits.dateTime().toPyDateTime()
        
        setup.v_sound = self.tryFloat(self.v_sound_edits.text())

        setup.t0 = self.tryFloat(self.t0_edits.text())
        setup.v = self.tryFloat(self.v_edits.text())
        setup.azim = self.tryFloat(self.azim_edits.text())
        setup.zangle = self.tryFloat(self.zangle_edits.text())
        setup.lat_i = self.tryFloat(self.lat_i_edits.text())
        setup.lon_i = self.tryFloat(self.lon_i_edits.text())
        setup.elev_i = self.tryFloat(self.elev_i_edits.text())
        setup.lat_f = self.tryFloat(self.lat_f_edits.text())
        setup.lon_f = self.tryFloat(self.lon_f_edits.text())
        setup.elev_f = self.tryFloat(self.elev_f_edits.text())

        setup.show_ballistic_waveform = self.show_ballistic_waveform_edits.currentText()

        setup.manual_fragmentation_search = []
        setup.fragmentation_point = self.fromTable(self.fragmentation_point)
        setup.show_fragmentation_waveform = self.show_fragmentation_waveform_edits.currentText()

        setup.manual_fragmentation_search.append(self.tryFloat(self.lat_frag_edits.text()))
        setup.manual_fragmentation_search.append(self.tryFloat(self.lon_frag_edits.text()))
        setup.manual_fragmentation_search.append(self.tryFloat(self.elev_frag_edits.text()))
        setup.manual_fragmentation_search.append(self.tryFloat(self.time_frag_edits.text()))

        setup.v_fixed = self.tryFloat(self.v_fixed_edits.text())
        setup.max_error = self.tryFloat(self.max_error_edits.text())
        setup.enable_restricted_time = self.restricted_time_check.isChecked()
        setup.restricted_time = self.restricted_time_edits.dateTime().toPyDateTime()
        setup.traj_tol = self.tryFloat(self.traj_tol_edits.text())
        setup.restrict_to_trajectory = self.restrict_to_trajectory_edits.currentText()

        setup.search_area = []
        setup.search_height = []
        setup.azimuth_min = self.tryFloat(self.azimuth_min_edits.text())
        setup.azimuth_max = self.tryFloat(self.azimuth_max_edits.text())
        setup.zangle_min = self.tryFloat(self.zangle_min_edits.text())
        setup.zangle_max = self.tryFloat(self.zangle_max_edits.text())
        setup.x_min = self.tryFloat(self.x_min_edits.text())
        setup.x_max = self.tryFloat(self.x_max_edits.text())
        setup.y_min = self.tryFloat(self.y_min_edits.text())
        setup.y_max = self.tryFloat(self.y_max_edits.text())
        setup.t_min = self.tryFloat(self.t_min_edits.text())
        setup.t_max = self.tryFloat(self.t_max_edits.text())
        setup.v_min = self.tryFloat(self.v_min_edits.text())
        setup.v_max = self.tryFloat(self.v_max_edits.text())
        setup.weight_distance_min = self.tryFloat(self.weight_distance_min_edits.text())
        setup.weight_distance_max = self.tryFloat(self.weight_distance_max_edits.text())
        setup.min_time = self.tryFloat(self.search_time_min_edits.text())
        setup.max_time = self.tryFloat(self.search_time_max_edits.text())
        setup.search_area.append(self.tryFloat(self.search_lat_min_edits.text()))
        setup.search_area.append(self.tryFloat(self.search_lat_max_edits.text()))
        setup.search_area.append(self.tryFloat(self.search_lon_min_edits.text()))
        setup.search_area.append(self.tryFloat(self.search_lon_max_edits.text()))
        setup.search_height.append(self.tryFloat(self.search_elev_min_edits.text()))
        setup.search_height.append(self.tryFloat(self.search_elev_max_edits.text()))

        setup.enable_winds = self.enable_winds_edits.currentText()
        setup.weather_type = self.weather_type_edits.currentText()
        setup.grid_size = self.tryFloat(self.grid_size_edits.text())

        setup.perturb_times = self.tryInt(self.perturb_times_edits.text())

        setup.observe_frag_no = self.tryInt(self.frag_no_edits.text())
        setup.perturb = self.perturb_edits.currentText()
        setup.perturb_method = self.perturb_method_edits.currentText()

        setup.fast_ballistic = self.fast_ballistic_edits.currentText()
        setup.fit_type = self.tryInt(self.fit_type_edits.text())
        setup.n_theta = self.tryInt(self.n_theta_edits.text())
        setup.n_phi = self.tryInt(self.n_phi_edits.text())
        setup.angle_precision = self.tryFloat(self.angle_precision_edits.text())
        setup.angle_error_tol = self.tryFloat(self.angle_error_tol_edits.text())

        setup.maxiter = self.tryInt(self.maxiter_edits.text())
        setup.swarmsize = self.tryInt(self.swarmsize_edits.text())
        setup.run_times = self.tryInt(self.run_times_edits.text())
        setup.minfunc = self.tryFloat(self.minfunc_edits.text())
        setup.minstep = self.tryFloat(self.minstep_edits.text())
        setup.phip = self.tryFloat(self.phip_edits.text())
        setup.phig = self.tryFloat(self.phig_edits.text())
        setup.omega = self.tryFloat(self.omega_edits.text())
        setup.pso_debug = self.pso_debug_edits.currentText()

        setup.plot_all_stations = self.plot_all_stations_edits.currentText()
        setup.colortoggle = self.color_toggle_edits.currentText()
        setup.dot_tol = self.tryInt(self.dot_tol_edits.text())
        setup.contour_res = self.tryInt(self.contour_res_edits.text())
        setup.high_f = self.high_f_edits.text()
        setup.high_b = self.high_b_edits.text()
        setup.rm_stat = self.rm_stat_edits.text()
        setup.img_dim = self.tryInt(self.img_dim_edits.text())
        setup.reported_points = self.fromTable(self.reported_points)

        setup.stations = self.fromTable(self.extra_point)

        # setup.arrival_times_file = os.path.join(setup.working_directory, setup.fireball_name, setup.arrival_times_file)
        # setup.sounding_file = os.path.join(setup.working_directory, setup.fireball_name, setup.sounding_file)
        # setup.perturbation_spread_file = os.path.join(setup.working_directory, setup.fireball_name, setup.perturbation_spread_file)
        # setup.station_picks_file = os.path.join(setup.working_directory, setup.fireball_name, setup.station_picks_file)
        # setup.replot_points_file = os.path.join(setup.working_directory, setup.fireball_name, setup.replot_points_file)

        if write:
            if '.ini' not in dlg[0]:
                configWrite(dlg[0] + '.ini', setup)
            else:
                configWrite(dlg[0], setup)

        setup.restrict_to_trajectory = (setup.restrict_to_trajectory.lower() == 'true')
        setup.show_ballistic_waveform = (setup.show_ballistic_waveform.lower() == 'true')
        setup.show_fragmentation_waveform = (setup.show_fragmentation_waveform.lower() == 'true')    

        return setup



    def loadINI(self):
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setNameFilters(['INI File (*.ini)'])
        #filenames = QStringList()

        dlg.exec_()

        filename = dlg.selectedFiles()
        setup = configRead(filename[0])
        
        self.fireball_name_edits.setText(setup.fireball_name)
        self.comboSet(self.difference_filter_edits, setup.difference_filter_all)
        self.comboSet(self.get_data_edits, setup.get_data)
        self.comboSet(self.run_mode_edits, setup.run_mode)
        self.comboSet(self.debug_edits, setup.debug)

        self.working_directory_edits.setText(setup.working_directory)
        self.arrival_times_edits.setText(setup.arrival_times_file)
        self.sounding_file_edits.setText(setup.sounding_file)
        self.perturbation_file_edits.setText(setup.perturbation_spread_file)
        self.station_picks_file_edits.setText(setup.station_picks_file)
        self.points_name_edits.setText(setup.replot_points_file)

        self.lat_centre_edits.setText(str(setup.lat_centre))
        self.lon_centre_edits.setText(str(setup.lon_centre))
        self.deg_radius_edits.setText(str(setup.deg_radius))
        self.start_datetime_edits.setDateTime(setup.start_datetime)
        self.end_datetime_edits.setDateTime(setup.end_datetime)
        self.v_sound_edits.setText(str(setup.v_sound))

        self.t0_edits.setText(str(setup.t0))
        self.v_edits.setText(str(setup.v))
        self.azim_edits.setText(str(setup.azim))
        self.zangle_edits.setText(str(setup.zangle))
        self.lat_i_edits.setText(str(setup.lat_i))
        self.lon_i_edits.setText(str(setup.lon_i))
        self.elev_i_edits.setText(str(setup.elev_i))
        self.lat_f_edits.setText(str(setup.lat_f))
        self.lon_f_edits.setText(str(setup.lon_f))
        self.elev_f_edits.setText(str(setup.elev_f))
        self.comboSet(self.show_ballistic_waveform_edits, setup.show_ballistic_waveform)

        self.toTable(self.fragmentation_point, setup.fragmentation_point)
        self.comboSet(self.show_fragmentation_waveform_edits, setup.show_fragmentation_waveform)
        self.lat_frag_edits.setText(str(setup.manual_fragmentation_search[0]))
        self.lon_frag_edits.setText(str(setup.manual_fragmentation_search[1]))
        self.elev_frag_edits.setText(str(setup.manual_fragmentation_search[2]))
        self.time_frag_edits.setText(str(setup.manual_fragmentation_search[3]))

        self.v_fixed_edits.setText(str(setup.v_fixed))
        self.max_error_edits.setText(str(setup.max_error))
        self.restricted_time_check.setChecked(setup.enable_restricted_time)
        self.restricted_time_edits.setDateTime(setup.restricted_time)
        self.traj_tol_edits.setText(str(setup.traj_tol))
        self.comboSet(self.restrict_to_trajectory_edits, setup.restrict_to_trajectory)

        self.azimuth_min_edits.setText(str(setup.azimuth_min))
        self.azimuth_max_edits.setText(str(setup.azimuth_max))
        self.zangle_min_edits.setText(str(setup.zangle_min))
        self.zangle_max_edits.setText(str(setup.zangle_max))
        self.x_min_edits.setText(str(setup.x_min))
        self.x_max_edits.setText(str(setup.x_max))
        self.y_min_edits.setText(str(setup.y_min))
        self.y_max_edits.setText(str(setup.y_max))
        self.t_min_edits.setText(str(setup.t_min))
        self.t_max_edits.setText(str(setup.t_max))
        self.v_min_edits.setText(str(setup.v_min))
        self.v_max_edits.setText(str(setup.v_max))
        self.weight_distance_min_edits.setText(str(setup.weight_distance_min))
        self.weight_distance_max_edits.setText(str(setup.weight_distance_max))
        self.search_time_min_edits.setText(str(setup.min_time))
        self.search_time_max_edits.setText(str(setup.max_time))

        if len(setup.search_area) > 0:
            self.search_lat_min_edits.setText(str(setup.search_area[0]))
            self.search_lat_max_edits.setText(str(setup.search_area[1]))
            self.search_lon_min_edits.setText(str(setup.search_area[2]))
            self.search_lon_max_edits.setText(str(setup.search_area[3]))
        else:
            print("Warning: Unable to detect search_area")

        if len(setup.search_area) > 0:
            self.search_elev_min_edits.setText(str(setup.search_height[0]))
            self.search_elev_max_edits.setText(str(setup.search_height[1]))
        else:
            print("Warning: Unable to detect search_height")

        self.comboSet(self.enable_winds_edits, setup.enable_winds)
        self.comboSet(self.weather_type_edits, setup.weather_type)
        self.grid_size_edits.setText(str(setup.grid_size))

        self.perturb_times_edits.setText(str(setup.perturb_times))
        self.frag_no_edits.setText(str(setup.observe_frag_no))
        self.comboSet(self.perturb_edits, setup.perturb)
        self.comboSet(self.perturb_method_edits, setup.perturb_method)

        self.comboSet(self.fast_ballistic_edits, setup.fast_ballistic)
        self.fit_type_edits.setText(str(setup.fit_type))
        self.n_theta_edits.setText(str(setup.n_theta))
        self.n_phi_edits.setText(str(setup.n_phi))
        self.angle_precision_edits.setText(str(setup.angle_precision))
        self.angle_error_tol_edits.setText(str(setup.angle_error_tol))

        self.maxiter_edits.setText(str(setup.maxiter))
        self.swarmsize_edits.setText(str(setup.swarmsize))
        self.run_times_edits.setText(str(setup.run_times))
        self.minfunc_edits.setText(str(setup.minfunc))
        self.minstep_edits.setText(str(setup.minstep))
        self.phip_edits.setText(str(setup.phip))
        self.phig_edits.setText(str(setup.phig))
        self.omega_edits.setText(str(setup.omega))
        self.comboSet(self.pso_debug_edits, setup.pso_debug)

        self.comboSet(self.plot_all_stations_edits, setup.plot_all_stations)
        self.comboSet(self.color_toggle_edits, setup.colortoggle)
        self.dot_tol_edits.setText(str(setup.dot_tol))
        self.contour_res_edits.setText(str(setup.contour_res))
        self.high_f_edits.setText(str(setup.high_f))
        self.high_b_edits.setText(str(setup.high_b))
        self.rm_stat_edits.setText(str(setup.rm_stat))
        self.img_dim_edits.setText(str(setup.img_dim))
        self.toTable(self.reported_points, setup.reported_points)

        self.toTableFromStn(self.extra_point, setup.stations)

    def createGrad(self, resid):

        max_error = max(resid)
        min_error = min(resid)

        error_range = max_error - min_error 

        resid = (resid - min_error)/error_range

        return resid

    def scatterPlot(self, setup, results, n_stations, xstn, s_name, dataset, manual=True):

        """ Outputs a scatter plot of the search data and the optimal supracenter

        Arguments:
            single_point: [list] lat, lon, and height of a manual search point
            n_stations: [int] number of stations
            xstn: [list] lat, lon, and height of all the stations
            s_name: [list] name of all the stations
            r: [list] residuals to each station
            x_opt: [list] lat, lon, and height of the optimal supracenter
            reported_points: [list] name, lat, lon, and height of extra reference points to be plotted
            search: [list] min/max lat, lon, and height of the search boundary
            output_name: [string] folder name to save output files in
            ref_pos: [list] mean position of the stations to act as the origin for the local coordinate system
            sup: [ndarray] array of points searched with the search algorithm used for plotting
            errors: [ndarray] array of the errors in each of the points in sup for plotting

        Returns:
            min_search: [list] min lat, lon and height of the search area
            max_search: [list] max lat, lon and height of the search area
        """


        r = results.r
        x_opt = results.x_opt
        ref_pos = [setup.ref_pos.lat, setup.ref_pos.lon, setup.ref_pos.elev]
        sup = results.sup
        errors = results.errors

        plt.style.use('dark_background')
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(20.9, 11.7)
        ax = fig.add_subplot(1, 1, 1, projection='3d')


        ### Labels
        ax.set_title("Supracenter Locations")
        ax.set_xlabel("Latitude (deg N)", linespacing=3.1)
        ax.set_ylabel("Longitude (deg E)", linespacing=3.1)
        ax.set_zlabel('Elevation (m)', linespacing=3.1)

        # plot station names and residuals
        for h in range(n_stations):

            # Convert station locations to geographic
            xstn[h, 0], xstn[h, 1], xstn[h, 2] = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], xstn[h, :])

            # Add station names
            ax.text(xstn[h, 0], xstn[h, 1], xstn[h, 2],  '%s' % (s_name[h]), size=10, zorder=1, color='w')

        # Add stations with color based off of residual
        res = ax.scatter(xstn[:, 0], xstn[:, 1], xstn[:, 2], c=abs(r), marker='^', cmap='viridis_r', depthshade=False)

        # Add point and label
        ax.scatter(x_opt[0], x_opt[1], x_opt[2], c = 'r', marker='*')
        ax.text(x_opt[0], x_opt[1], x_opt[2], '%s' % ('Supracenter'), zorder=1, color='w')

        if not manual:
            #for i in range(len(sup)):
                #sup[i, 0], sup[i, 1], sup[i, 2] = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], sup[i, :])
            sc = ax.scatter(sup[:, 0], sup[:, 1], sup[:, 2], c=errors, cmap='inferno_r', depthshade=False)
            a = plt.colorbar(sc, ax=ax)
            a.set_label("Error in Supracenter (s)")

        # colorbars
        b = plt.colorbar(res, ax=ax)
        b.set_label("Station Residuals (s)")

        if manual:
            try:
                self.plots.removeWidget(self.threelbar)
            except:
                pass
            self.plots.removeWidget(self.three_canvas)
            self.three_canvas = FigureCanvas(Figure(figsize=(3, 3)))
            self.three_canvas = FigureCanvas(fig)
            self.three_canvas.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
            self.threelbar = NavigationToolbar(self.three_canvas, self)
            self.plots.addWidget(self.three_canvas)
            self.plots.addWidget(self.threelbar)
            self.three_canvas.draw()
        else:
            try:
                self.sup_plots.removeWidget(self.sup_threelbar)
            except:
                pass
            self.sup_plots.removeWidget(self.sup_three_canvas)
            self.sup_three_canvas = FigureCanvas(Figure(figsize=(3, 3)))
            self.sup_three_canvas = FigureCanvas(fig)
            self.sup_three_canvas.setSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred)
            self.sup_threelbar = NavigationToolbar(self.sup_three_canvas, self)
            self.sup_plots.addWidget(self.sup_three_canvas)
            self.sup_plots.addWidget(self.sup_threelbar)
            self.sup_three_canvas.draw()
        ax.mouse_init()
        SolutionGUI.update(self)

    def residPlot(self, results_arr, s_name, xstn, output_name, n_stations, manual=True):
        """ outputs a 2D residual plot of the stations with the optimal supracenter

        Arguments:
            x_opt: [list] optimal supracenter position
            s_name: [list] list of station names
            xstn: [list] list of lat, lon, height positions of each station
            resid: [list] list of residuals to each station
            output_name: [string] folder to store the data in
            n_stations: [int] number of stations
        """
        
        x_opt = results_arr.x_opt
        resid = results_arr.r

        plt.style.use('dark_background')
        fig = plt.figure(figsize=plt.figaspect(0.5))
        fig.set_size_inches(20.9, 11.7)
        ax = fig.add_subplot(1, 1, 1)
        res = ax.scatter(xstn[:, 0], xstn[:, 1], c=abs(resid), marker='^', cmap='viridis_r', s=21)
        
        for h in range(n_stations):

            # Add station names
            ax.text(xstn[h, 0], xstn[h, 1],  '%s' % (s_name[h]), size=10, zorder=1, color='w')

        # Add point and label
        ax.scatter(x_opt[0], x_opt[1], c = 'r', marker='*', s=21)
        ax.text(x_opt[0], x_opt[1],  '%s' % ('Supracenter'), size=10, zorder=1, color='w')

        ax.set_xlabel("Latitude (deg N)")
        ax.set_ylabel("Longitude (deg E)")
        ax.set_title("Station Residuals")

        c = plt.colorbar(res, ax=ax)
        c.set_label("Station Residuals (s)")

        if manual:
            self.plots.removeWidget(self.two_canvas)
            try:
                self.plots.removeWidget(self.twolbar)
            except:
                pass
            self.two_canvas = FigureCanvas(Figure(figsize=(3, 3)))
            self.two_canvas = FigureCanvas(fig)
            self.two_canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
            self.twolbar = NavigationToolbar(self.two_canvas, self)
            self.plots.addWidget(self.two_canvas)    
            self.plots.addWidget(self.twolbar)
            self.two_canvas.draw()
        else:

            self.sup_plots.removeWidget(self.sup_two_canvas)
            try:
                self.sup_plots.removeWidget(self.sup_twolbar)
            except:
                pass
            self.sup_two_canvas = FigureCanvas(Figure(figsize=(3, 3)))
            self.sup_two_canvas = FigureCanvas(fig)
            self.sup_two_canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
            self.sup_twolbar = NavigationToolbar(self.sup_two_canvas, self)
            self.sup_plots.addWidget(self.sup_two_canvas)    
            self.sup_plots.addWidget(self.sup_twolbar)
            self.sup_two_canvas.draw()
        SolutionGUI.update(self)

    def tabChange(self):

        setup = self.saveINI(False)

        # manual Supracenter
        self.lat_edit.setText(str(setup.manual_fragmentation_search[0]))
        self.lon_edit.setText(str(setup.manual_fragmentation_search[1]))
        self.elev_edit.setText(str(setup.manual_fragmentation_search[2]))
        self.time_edit.setText(str(setup.manual_fragmentation_search[3]))
        self.ref_edit.setDateTime(setup.start_datetime)
        self.picks_file_edit.setText(os.path.join(setup.working_directory, setup.fireball_name, setup.station_picks_file))
        self.atmospheric_file_edit.setText(os.path.join(setup.working_directory, setup.fireball_name, setup.sounding_file))

        # atmospheric profile
        self.atm_atm_file_edits.setText(os.path.join(setup.working_directory, setup.fireball_name, setup.sounding_file))
        self.comboSet(self.atm_weather_type_edits, setup.weather_type)
        self.atm_perturbation_file_edits.setText(os.path.join(setup.working_directory, setup.fireball_name, setup.perturbation_spread_file))
        self.atm_perturb_times_edits.setText(str(setup.perturb_times))
        self.comboSet(self.atm_perturb_method_edits, setup.perturb_method)

        #PSO Supracenter
        self.sup_south_edits.setText(str(setup.search_area[0]))
        self.sup_north_edits.setText(str(setup.search_area[1]))
        self.sup_west_edits.setText(str(setup.search_area[2]))
        self.sup_east_edits.setText(str(setup.search_area[3]))
        self.sup_height_min_edits.setText(str(setup.search_height[0]))
        self.sup_height_max_edits.setText(str(setup.search_height[1]))
        self.sup_min_time_edits.setText(str(setup.min_time))
        self.sup_max_time_edits.setText(str(setup.max_time))
        self.sup_picks_file_edit.setText(os.path.join(setup.working_directory, setup.fireball_name, setup.station_picks_file))
        self.sup_ref_time_edits.setDateTime(setup.start_datetime)

if __name__ == '__main__':

    app = QApplication(sys.argv)

    gui = SolutionGUI()

    w = 1280; h = 1024
    gui.resize(w, h)
    gui.show()

    app.exec_()