
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time
import datetime

from netCDF4 import Dataset

from PyQt5.QtWidgets import *
import sys
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph as pg
import pyqtgraph.opengl as gl

from functools import partial

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from mpl_toolkits.mplot3d import Axes3D

from supra.Supracenter.stationDat import convStationDat
from supra.Supracenter.psoSearch import psoSearch
from supra.Fireballs.Program import position, configRead, configWrite
from supra.Fireballs.SeismicTrajectory import Constants, parseWeather
from supra.Supracenter.angleConv import loc2Geo, geo2Loc, angle2NDE
from supra.Supracenter.netCDFconv import findECMWFSound
from supra.Supracenter.SPPT import perturb as perturbation_method
from supra.Supracenter.fetchCopernicus import copernicusAPI

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

        self.tab_widget = QTabWidget()
        self.tab_widget.blockSignals(True)
        self.tab_widget.currentChanged.connect(self.tabChange)
        # self.tab_widget.blockSignals(False)

        self.addIniWidgets()
        self.addSupraWidgets()
        self.addSupWidgets()
        self.addFetchATMWidgets()
        self.addProfileWidgets()
        self.addDocsWidgets()

        self.var_typ = 't'

        self.tab_widget.blockSignals(False)
        layout.addWidget(self.tab_widget, 1, 1)

        # extractAction = QtGui.QAction("&GET TO THE CHOPPAH!!!", self)
        # extractAction.setShortcut("Ctrl+Q")
        # extractAction.setStatusTip('Leave The App')
        # extractAction.triggered.connect(self.close_application)

        # self.statusBar()

        # mainMenu = self.menuBar()
        # fileMenu = mainMenu.addMenu('&File')
        # fileMenu.addAction(extractAction)

        menu_bar = self.menuBar() 
        layout.addWidget(menu_bar, 0, 1)
        fileMenu = menu_bar.addMenu('&File')
        aboutMenu = menu_bar.addMenu('&About')
        #fileMenu.addAction(exitAction)

        stylesheet = """ 
        QTabWidget>QWidget>QWidget{background: gray;}
        QLabel{color: white;}
        QCheckBox(color: white;)
        """

        self.setStyleSheet(stylesheet)

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

        setup.lat_centre = self.fatm_lat_slide.value()
        setup.lon_centre = self.fatm_lon_slide.value()
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
            print('error, atmPlotProfile')
            print(self.fatm_variable_combo.currentText())

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
            obj.setText('Latitude: {:.1f}'.format(slider.value()))
        elif obj == self.fatm_lon_label:
            obj.setText('Longitude: {:.1f}'.format(slider.value()))
        else:
            self.errorMessage(self, 'Bad atm slider pass in fatmValueChange', 2)

        #self.atmPlotProfile(self.atm_lat_slide.value(), self.atm_lon_slide.value(), self.var_typ)
    
    def fatmPrint(self):
        
        filename = QFileDialog.getSaveFileName(self, 'Save File', '', 'Text File (*.txt)')

        consts = Constants()
        setup = self.saveINI(False)

        setup.lat_centre = self.fatm_lat_slide.value()
        setup.lon_centre = self.fatm_lon_slide.value()
        setup.sounding_file = self.fatm_name_edits.text()
        setup.weather_type = 'ecmwf'

        try:
            dataset = parseWeather(setup, consts)
            sounding = findECMWFSound(setup.lat_centre, setup.lon_centre, dataset)
        except:
            self.errorMessage('Error reading weather profile in fatmPrint', 2)
            return None

        if '.txt' not in filename[0]:
            filename[0] = filename[0] + '.txt'

        with open(str(filename[0]), 'w') as f:
            
            f.write('Height (m) , Speed of Sound (m/s), Wind Magnitude (m/s), Wind Direction (deg f N)')

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
        self.fatm_lat_slide.setMinimum(-90)
        self.fatm_lat_slide.setMaximum(90)
        self.fatm_lat_slide.setValue(0)
        self.fatm_lat_slide.setTickInterval(0.5)
        self.fatm_lat_slide.valueChanged.connect(partial(self.fatmValueChange, self.fatm_lat_label, self.fatm_lat_slide))

        self.fatm_lon_label = QLabel("Longitude: 0")
        self.fatm_lon_slide = QSlider(Qt.Horizontal)
        fetch_content.addWidget(self.fatm_lon_label, 8, 1)
        fetch_content.addWidget(self.fatm_lon_slide, 8, 2,)
        self.fatm_lon_slide.setMinimum(-180)
        self.fatm_lon_slide.setMaximum(180)
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


        # folder1_2 = QTreeWidgetItem(folder1, ["Living room", "Approved by client"])
        # folder1_2.setData(2, Qt.EditRole, '') 

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
            print('error, atmPlotProfile')

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
            obj.setText('Latitude: {:.1f}'.format(slider.value()))
        elif obj == self.atm_lon_label:
            obj.setText('Longitude: {:.1f}'.format(slider.value()))
        else:
            self.errorMessage(self, 'Bad atm slider pass in atmValueChange', 2)

        self.atmPlotProfile(self.atm_lat_slide.value(), self.atm_lon_slide.value(), self.var_typ)
    
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
        self.atm_lat_slide.setMinimum(-90)
        self.atm_lat_slide.setMaximum(90)
        self.atm_lat_slide.setValue(0)
        self.atm_lat_slide.setTickInterval(0.5)
        self.atm_lat_slide.valueChanged.connect(partial(self.atmValueChange, self.atm_lat_label, self.atm_lat_slide))

        self.atm_lon_label = QLabel("Longitude: 0")
        self.atm_lon_slide = QSlider(Qt.Horizontal)
        profile_tab_content.addWidget(self.atm_lon_label, 7, 1)
        profile_tab_content.addWidget(self.atm_lon_slide, 7, 2, 1, 3)
        self.atm_lon_slide.setMinimum(-180)
        self.atm_lon_slide.setMaximum(180)
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
        self.atm_T_button.clicked.connect(partial(self.atmPlotProfile, self.atm_lat_slide.value(), self.atm_lon_slide.value(), var_typ='t'))

        self.atm_mag_button = QPushButton('Wind Magnitude')
        profile_tab_content_graph.addWidget(self.atm_mag_button)
        self.atm_mag_button.clicked.connect(partial(self.atmPlotProfile, self.atm_lat_slide.value(), self.atm_lon_slide.value(), var_typ='m'))

        self.atm_dir_button = QPushButton('Wind Direction')
        profile_tab_content_graph.addWidget(self.atm_dir_button)
        self.atm_dir_button.clicked.connect(partial(self.atmPlotProfile, self.atm_lat_slide.value(), self.atm_lon_slide.value(), var_typ='d'))

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
        
        X = len(table)
        Y = len(table[0])

        obj.setRowCount(X)
        obj.setColumnCount(Y)

        for x in range(X):
            for y in range(Y):
                obj.setItem(x, y, QTableWidgetItem(str(table[x][y])))

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
        setup.dot_tol = self.tryFloat(self.dot_tol_edits.text())
        setup.contour_res = self.tryFloat(self.contour_res_edits.text())
        setup.high_f = self.high_f_edits.text()
        setup.high_b = self.high_b_edits.text()
        setup.rm_stat = self.rm_stat_edits.text()
        setup.img_dim = self.tryFloat(self.img_dim_edits.text())
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

        self.search_lat_min_edits.setText(str(setup.search_area[0]))
        self.search_lat_max_edits.setText(str(setup.search_area[1]))
        self.search_lon_min_edits.setText(str(setup.search_area[2]))
        self.search_lon_max_edits.setText(str(setup.search_area[3]))
        self.search_elev_min_edits.setText(str(setup.search_height[0]))
        self.search_elev_max_edits.setText(str(setup.search_height[1]))

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
            for i in range(len(sup)):
                print(sup)
                sup[i, 0], sup[i, 1], sup[i, 2] = loc2Geo(ref_pos[0], ref_pos[1], ref_pos[2], sup[i, :])
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
            self.three_canvas.setSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred)
            self.threelbar = NavigationToolbar(self.three_canvas, self)
            self.plots.addWidget(self.three_canvas)
            self.plots.addWidget(self.threelbar)
            self.three_canvas.draw()
        else:
            try:
                self.sup_plots.removeWidget(self.sup_threelbar)
            except:
                pass
            self.sup_plots.removeWidget(self.three_canvas)
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