
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
from wmpl.Utils.TrajConversions import datetime2JD

global arrTimes 
global sounding

DATA_FILE = 'data.txt'
OUTPUT_CSV = 'data_picks.csv'

class Pick:
    def __init__(self, time, stn, stn_no):
        self.time = time
        self.stn = stn
        self.stn_no = stn_no

             
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
        self.addPicksReadWidgets()
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

    def addPicksReadWidgets(self):
        picks_read_tab = QWidget()
        picks_read_tab_content = QGridLayout()
        picks_read_tab.setLayout(picks_read_tab_content)

        self.tab_widget.addTab(picks_read_tab, "Picks Read")

        self.csv_table = QTableWidget(0, 9)
        picks_read_tab_content.addWidget(self.csv_table, 1, 1, 1, 4)
        self.csv_table.setHorizontalHeaderLabels(['Pick Group', 'Network', 'Code', 'Latitude', 'Longitude', 'Elevation', 'Pick JD', 'Pick Time', 'station_number'])
        header = self.csv_table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.Stretch)

        self.csv_table_add = QPushButton("+")
        picks_read_tab_content.addWidget(self.csv_table_add, 2, 2, 1, 1)
        self.csv_table_add.clicked.connect(partial(self.changeRows, self.csv_table, 1))
        self.csv_table_add.setToolTip("Add row")

        self.csv_table_min = QPushButton("-")
        picks_read_tab_content.addWidget(self.csv_table_min, 2, 1, 1, 1)
        self.csv_table_min.clicked.connect(partial(self.changeRows, self.csv_table, -1))
        self.csv_table_min.setToolTip("Remove row")

        self.csv_table_load = QPushButton("Load")
        picks_read_tab_content.addWidget(self.csv_table_load, 2, 3, 1, 1)
        self.csv_table_load.clicked.connect(self.csvLoad)

        self.csv_table_save = QPushButton("Save")
        picks_read_tab_content.addWidget(self.csv_table_save, 2, 4, 1, 1)
        self.csv_table_save.clicked.connect(self.csvSave)

    def csvLoad(self):
        
        dlg = QFileDialog()
        dlg.setFileMode(QFileDialog.AnyFile)
        dlg.setNameFilters(['CSV File (*.csv)'])
        dlg.exec_()

        filename = dlg.selectedFiles()

        data_table = []

        with open(filename[0]) as f:
            
            next(f)
            
            for line in f:
                data_table.append(line.split(','))

        self.toTable(self.csv_table, data_table)


    def csvSave(self):

        dlg = QFileDialog.getSaveFileName(self, 'Save File')


        if '.csv' not in dlg[0]:
            file_name = dlg[0] + '.csv'
        else:
            file_name = dlg[0]

        data_set = self.fromTable(self.csv_table)
        # Open the output CSV
        with open(os.path.join(file_name), 'w') as f:

            # Write the header
            f.write('Pick group, Network, Code, Lat, Lon, Elev, Pick JD, Pick time, station_number \n')

            # Go through all picks
            for line in data_set:

                # Write the CSV entry
                f.write("{:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}, {:}\n".format(*line))

        self.errorMessage('Output to CSV!', 0, title='Exported!')


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
        self.station_picks_edits.setText(self.picks_file_edit.text())
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

        self.three_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.three_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.plots.addWidget(self.three_canvas)

        self.search_button = QPushButton('Search')
        self.supra_tab_content.addWidget(self.search_button, 11, 3, 1, 3)
        self.search_button.clicked.connect(self.supraSearch)

        self.lat_label, self.lat_edit = self.createLabelEditObj("Latitude: ", self.supra_tab_content, 1, h_shift=2, width=2)
        self.lon_label, self.lon_edit = self.createLabelEditObj("Longitude: ", self.supra_tab_content, 2, h_shift=2, width=2)
        self.elev_label, self.elev_edit = self.createLabelEditObj("Elevation: ", self.supra_tab_content, 3, h_shift=2, width=2)
        self.time_label, self.time_edit = self.createLabelEditObj("Time: ", self.supra_tab_content, 4, h_shift=2, width=2)

        self.supra_save_changes_button = QPushButton('Copy to INI Builder')
        self.supra_tab_content.addWidget(self.supra_save_changes_button, 8, 3, 1, 1)
        self.supra_save_changes_button.clicked.connect(self.supraSaveChanges)

        self.results_label = QLabel("Results: ")
        self.supra_tab_content.addWidget(self.results_label, 9, 1, 1, 1)

        self.tableWidget = QTableWidget(0, 0)
        self.supra_tab_content.addWidget(self.tableWidget, 10, 1, 1, 10)

        self.ref_label = QLabel("Reference Datetime: ")
        self.supra_tab_content.addWidget(self.ref_label, 5, 3, 1, 1)
        self.ref_edit = QDateTimeEdit()
        self.supra_tab_content.addWidget(self.ref_edit, 5, 4, 1, 2)
        self.ref_edit.setCalendarPopup(True)
        
        self.picks_file_label, self.picks_file_edit, self.picks_file_buton = self.createFileSearchObj("Picks File: ", self.supra_tab_content, 6, h_shift=2)
        self.picks_file_buton.clicked.connect(partial(self.fileSearch, ["CSV Picks File (*.csv)"], self.picks_file_edit))

        self.atmospheric_file_label, self.atmospheric_file_edit, self.atmospheric_file_buton = self.createFileSearchObj("Atmospheric File: ", self.supra_tab_content, 7, h_shift=2)
        self.atmospheric_file_buton.clicked.connect(partial(self.fileSearch, ["NetCDF (*.nc)"], self.atmospheric_file_edit))
        
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
            
        sounding = self.parseGeneralECMWF(setup.sounding_file, setup.lat_centre, setup.lon_centre, atm_time, variables)

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

        self.fatm_name_label, self.fatm_name_edits = self.createLabelEditObj('Name:', fetch_content, 1)

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


    def createLabelEditObj(self, label_name, parent, row, width=1, h_shift=0, tool_tip=''):
        """ Creates a label and line edit object beside each other

        Arguments:
        label_name [String]: string to be displayed in the label
        parent [obj]: the object to hold this item. Must be in QGridLayout
        row [int]: the row the object will sit on
        
        Keyword Arguments:
        width [int]: width of the lineedit box
        h_shift [int]: horizontal position to place the combination
        tool_tip [String]: identifier name in tool_tips.csv that the tool tip is under

        Returns:
        label_obj, edits_obj [obj]: the label and lineedit objects
        """

        label_obj = QLabel(label_name)
        edits_obj = QLineEdit('')
        parent.addWidget(label_obj, row, 1 + h_shift)
        parent.addWidget(edits_obj, row, 2 + h_shift, 1, width)

        if tool_tip != '':
            label_obj.setToolTip(self.toolTime(tool_tip))

        return label_obj, edits_obj

    def createFileSearchObj(self, label_name, parent, row, width=1, h_shift=0, tool_tip=''):

        label_obj, edits_obj = self.createLabelEditObj(label_name, parent, row, width=width, h_shift=h_shift, tool_tip=tool_tip)

        buton_obj = QPushButton('Browse')
        parent.addWidget(buton_obj, row, 3 + h_shift)

        return label_obj, edits_obj, buton_obj

    def createComboBoxObj(self, label_name, parent, row, items=[], width=1, h_shift=0, tool_tip=''):

        label_obj = QLabel(label_name)
        combo_obj = QComboBox()
        parent.addWidget(label_obj, row, 1 + h_shift)
        parent.addWidget(combo_obj, row, 2 + h_shift, 1, width)

        for item in items: combo_obj.addItem(item)

        if tool_tip != '':
            label_obj.setToolTip(self.toolTime(tool_tip))

        return label_obj, combo_obj

    def createLabelDateEditObj(self, label_name, parent, row, width=1, h_shift=0, tool_tip='', popup=True):
        
        label_obj = QLabel(label_name)
        dedit_obj = QDateTimeEdit()
        parent.addWidget(label_obj, row, 1 + h_shift)
        parent.addWidget(dedit_obj, row, 2 + h_shift, 1, width)
        dedit_obj.setCalendarPopup(popup)

        if tool_tip != '':
            label_obj.setToolTip(self.toolTime(tool_tip))

        return label_obj, dedit_obj


    def initGeneralTab(self):

        general = QWidget()
        general_content = QGridLayout()
        general.setLayout(general_content)
        self.section_widget.addTab(general, "General")

        self.fireball_name_label, self.fireball_name_edits = self.createLabelEditObj('Fireball Name:', general_content, 1, tool_tip='fireball_name')

        self.difference_filter_label, self.difference_filter_edits = self.createComboBoxObj('Difference Filter:', general_content, 2, items=['True', 'False'], tool_tip='difference_filter')
        self.get_data_label, self.get_data_edits = self.createComboBoxObj('Get Data: ', general_content, 3, items=['True', 'False'], tool_tip='get_data')
        self.run_mode_label, self.run_mode_edits = self.createComboBoxObj('Run Mode: ', general_content, 4, items=['Search', 'Replot', 'Manual'], tool_tip='run_mode')
        self.debug_label, self.debug_edits = self.createComboBoxObj('Debug: ', general_content, 5, items=['True', 'False'], tool_tip='debug')

    def initFilesTab(self):
        files = QWidget()
        files_content = QGridLayout()
        files.setLayout(files_content)
        self.section_widget.addTab(files, "Files")

        self.working_directory_label, self.working_directory_edits, self.working_directory_buton = self.createFileSearchObj('Working Directory: ', files_content, 1, width=1, h_shift=0, tool_tip='working_directory')
        self.working_directory_buton.clicked.connect(partial(self.folderSearch, self.working_directory_edits))

        self.arrival_times_label, self.arrival_times_edits, self.arrival_times_buton = self.createFileSearchObj('Arrival Times:', files_content, 2, width=1, h_shift=0, tool_tip='arrival_times_file')
        self.arrival_times_buton.clicked.connect(partial(self.fileSearch, ['Numpy Array (*.npy)'], self.arrival_times_edits))

        self.sounding_file_label, self.sounding_file_edits, self.sounding_file_buton = self.createFileSearchObj('Sounding File:', files_content, 3, width=1, h_shift=0, tool_tip='sounding_file')
        self.sounding_file_buton.clicked.connect(partial(self.fileSearch, ['NetCDF (*.nc)', 'HDF (*.HDF)'], self.sounding_file_edits))

        self.perturbation_file_label, self.perturbation_file_edits, self.perturbation_file_buton = self.createFileSearchObj('Perturbation', files_content, 4, width=1, h_shift=0, tool_tip='perturbation_spread_file')
        self.perturbation_file_buton.clicked.connect(partial(self.fileSearch, ['NetCDF (*.nc)'], self.perturbation_file_edits))

        self.station_picks_label, self.station_picks_edits, self.station_picks_buton = self.createFileSearchObj('Station Picks File: ', files_content, 5, width=1, h_shift=0, tool_tip='station_picks_file')
        self.station_picks_buton.clicked.connect(partial(self.fileSearch, ['CSV (*.csv)', 'Text File (*.txt)'], self.station_picks_edits))

        self.points_name_label, self.points_name_edits, self.points_name_buton = self.createFileSearchObj('Replot Points File: ', files_content, 6, width=1, h_shift=0, tool_tip='points_name')
        self.points_name_buton.clicked.connect(partial(self.fileSearch, ['CSV (*.csv)'], self.points_name_edits))

    def initParametersTab(self):
        params = QWidget()
        params_content = QGridLayout()
        params.setLayout(params_content)
        self.section_widget.addTab(params, "Parameters")

        self.lat_centre_label, self.lat_centre_edits = self.createLabelEditObj('Latitude Center:', params_content, 1, tool_tip='lat_centre')
        self.lon_centre_label, self.lon_centre_edits = self.createLabelEditObj('Longitude Center:', params_content, 2, tool_tip='lon_centre')
        self.deg_radius_label, self.deg_radius_edits = self.createLabelEditObj('Degrees in Search Radius:', params_content, 3, tool_tip='deg_radius')

        self.start_datetime_label, self.start_datetime_edits = self.createLabelDateEditObj("Start Datetime", params_content, 4, tool_tip='start_datetime')
        self.end_datetime_label, self.end_datetime_edits = self.createLabelDateEditObj("End Datetime", params_content, 5, tool_tip='end_datetime')

        self.v_sound_label, self.v_sound_edits = self.createLabelEditObj('Average Speed of Sound:', params_content, 6, tool_tip='v_sound')

    def initBallisticTab(self):

        ballistic = QWidget()
        ballistic_content = QGridLayout()
        ballistic.setLayout(ballistic_content)
        self.section_widget.addTab(ballistic, "Ballistic")

        self.t0_label, self.t0_edits = self.createLabelEditObj('t0:', ballistic_content, 1, width=3, tool_tip='t0')
        self.v_label, self.v_edits = self.createLabelEditObj('v:', ballistic_content, 2, width=3, tool_tip='v')
        self.azim_label, self.azim_edits = self.createLabelEditObj('azim:', ballistic_content, 3, width=3, tool_tip='azim')
        self.zangle_label, self.zangle_edits = self.createLabelEditObj('zangle:', ballistic_content, 4, width=3, tool_tip='zangle')
        self.lat_i_label, self.lat_i_edits = self.createLabelEditObj('lat_i:', ballistic_content, 5, tool_tip='lat_i')
        self.lon_i_label, self.lon_i_edits = self.createLabelEditObj('lon_i:', ballistic_content, 6, tool_tip='lon_i')
        self.elev_i_label, self.elev_i_edits = self.createLabelEditObj('elev_i:', ballistic_content, 7, tool_tip='elev_i')
        self.lat_f_label, self.lat_f_edits = self.createLabelEditObj('lat_f:', ballistic_content, 5, h_shift=2, tool_tip='lat_f')
        self.lon_f_label, self.lon_f_edits = self.createLabelEditObj('lon_f:', ballistic_content, 6, h_shift=2, tool_tip='lon_f')
        self.elev_f_label, self.elev_f_edits = self.createLabelEditObj('elev_f:', ballistic_content, 7, h_shift=2, tool_tip='elev_f')

        self.show_ballistic_waveform_label, self.show_ballistic_waveform_edits = self.createComboBoxObj('Show Ballistic Waveform', ballistic_content, 8, items=['True', 'False'], width=2, tool_tip='show_ballistic_waveform')

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

        self.show_fragmentation_waveform_label, self.show_fragmentation_waveform_edits = self.createComboBoxObj("Show Fragmentation Waveform: ", fragmentation_content, 3, width=2, items=['True', 'False'], tool_tip='show_fragmentation_waveform')

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

        self.v_fixed_label, self.v_fixed_edits = self.createLabelEditObj('v_fixed:', restriction_content, 1, width=3, tool_tip='v_fixed')
        self.max_error_label, self.max_error_edits = self.createLabelEditObj('max_error:', restriction_content, 2, width=3, tool_tip='max_error')

        self.restricted_time_check = QCheckBox("Enable Restricted Time: ")
        restriction_content.addWidget(self.restricted_time_check, 3, 4, 1, 1)
        self.restricted_time_label, self.restricted_time_edits = self.createLabelDateEditObj("Restricted Time: ", restriction_content, 3, width=2, tool_tip='restricted_time')

        self.traj_tol_label, self.traj_tol_edits = self.createLabelEditObj('traj_tol:', restriction_content, 4, width=3, tool_tip='traj_tol')

        self.restrict_to_trajectory_label, self.restrict_to_trajectory_edits = self.createComboBoxObj('restrict_to_trajectory', restriction_content, 5, items=['True', 'False'], width=2, tool_tip='restrict_to_trajectory')

        self.azimuth_min_label, self.azimuth_min_edits = self.createLabelEditObj('azimuth_min:', restriction_content, 6, tool_tip='azimuth_min')
        self.azimuth_max_label, self.azimuth_max_edits = self.createLabelEditObj('azimuth_max:', restriction_content, 6, h_shift=2, tool_tip='azimuth_max')
        self.zangle_min_label, self.zangle_min_edits = self.createLabelEditObj('zangle_min:', restriction_content, 7, tool_tip='zenith_min')
        self.zangle_max_label, self.zangle_max_edits = self.createLabelEditObj('zangle_max:', restriction_content, 7, h_shift=2, tool_tip='zenith_max')
        self.x_min_label, self.x_min_edits = self.createLabelEditObj('x_min:', restriction_content, 8, tool_tip='x_min')
        self.x_max_label, self.x_max_edits = self.createLabelEditObj('x_max:', restriction_content, 8, h_shift=2, tool_tip='x_max')
        self.y_min_label, self.y_min_edits = self.createLabelEditObj('y_min:', restriction_content, 9, tool_tip='y_min')
        self.y_max_label, self.y_max_edits = self.createLabelEditObj('y_max:', restriction_content, 9, h_shift=2, tool_tip='y_max')
        self.t_min_label, self.t_min_edits = self.createLabelEditObj('t_min:', restriction_content, 10, tool_tip='t_min')
        self.t_max_label, self.t_max_edits = self.createLabelEditObj('t_max:', restriction_content, 10, h_shift=2, tool_tip='t_max')
        self.v_min_label, self.v_min_edits = self.createLabelEditObj('v_min:', restriction_content, 11, tool_tip='v_min')
        self.v_max_label, self.v_max_edits = self.createLabelEditObj('v_max:', restriction_content, 11, h_shift=2, tool_tip='v_max')
        self.weight_distance_min_label, self.weight_distance_min_edits = self.createLabelEditObj('weight_distance_min', restriction_content, 12, tool_tip='weight_distance_min')
        self.weight_distance_max_label, self.weight_distance_max_edits = self.createLabelEditObj('weight_distance_max', restriction_content, 12, h_shift=2, tool_tip='weight_distance_max')
        self.search_time_min_label, self.search_time_min_edits = self.createLabelEditObj('search_time_min', restriction_content, 13, tool_tip='min_time')
        self.search_time_max_label, self.search_time_max_edits = self.createLabelEditObj('search_time_max', restriction_content, 13, h_shift=2, tool_tip='max_time')
        self.search_lat_min_label, self.search_lat_min_edits = self.createLabelEditObj('search_lat_min', restriction_content, 14, tool_tip='search_area')
        self.search_lat_max_label, self.search_lat_max_edits = self.createLabelEditObj('search_lat_max', restriction_content, 14, h_shift=2, tool_tip='search_area')
        self.search_lon_min_label, self.search_lon_min_edits = self.createLabelEditObj('search_lon_min', restriction_content, 15, tool_tip='search_area')
        self.search_lon_max_label, self.search_lon_max_edits = self.createLabelEditObj('search_lon_max', restriction_content, 15, h_shift=2, tool_tip='search_area')
        self.search_elev_min_label, self.search_elev_min_edits = self.createLabelEditObj('search_elev_min', restriction_content, 16, tool_tip='search_area')
        self.search_elev_max_label, self.search_elev_max_edits = self.createLabelEditObj('search_elev_max', restriction_content, 16, h_shift=2, tool_tip='search_area')

    def initAtmosphereTab(self):

        atmosphere = QWidget()
        atmosphere_content = QGridLayout()
        atmosphere.setLayout(atmosphere_content)
        self.section_widget.addTab(atmosphere, "Atmosphere")

        self.enable_winds_label, self.enable_winds_edits = self.createComboBoxObj('Enable Winds: ', atmosphere_content, 1, items=['True', 'False'], tool_tip='enable_winds')
        self.weather_type_label, self.weather_type_edits = self.createComboBoxObj('Weather Type: ', atmosphere_content, 2, items=['none', 'ecmwf', 'ukmo', 'merra', 'custom'], tool_tip='weather_type')

        self.grid_size_label, self.grid_size_edits = self.createLabelEditObj('Grid Size', atmosphere_content, 3, tool_tip='grid_size')

    def initPerturbationsTab(self):

        perturb = QWidget()
        perturb_content = QGridLayout()
        perturb.setLayout(perturb_content)
        self.section_widget.addTab(perturb, "Perturbations")

        self.perturb_times_label, self.perturb_times_edits = self.createLabelEditObj('Perturbation Times', perturb_content, 1, tool_tip='perturb_times')
        self.frag_no_label, self.frag_no_edits = self.createLabelEditObj('Fragmentation Number', perturb_content, 2, tool_tip='fragno')

        self.perturb_label, self.perturb_edits = self.createComboBoxObj('Perturb: ', perturb_content, 3, items=['True', 'False'], tool_tip='perturb')
        self.perturb_method_label, self.perturb_method_edits = self.createComboBoxObj('Perturb Method', perturb_content, 4, items=['none', 'bmp', 'sppt', 'temporal', 'spread', 'spread_r'], tool_tip='perturb_method')


    def initSpeedTab(self):

        speed = QWidget()
        speed_content = QGridLayout()
        speed.setLayout(speed_content)
        self.section_widget.addTab(speed, "Speed")

        self.fast_ballistic_label, self.fast_ballistic_edits = self.createComboBoxObj("Fast Ballistic: ", speed_content, 1, items=['True', 'False'], tool_tip='fast_ballistic')

        self.fit_type_label, self.fit_type_edits = self.createLabelEditObj('Fit Type:', speed_content, 2, tool_tip='fit_type')
        self.n_theta_label, self.n_theta_edits = self.createLabelEditObj('Theta Resolution', speed_content, 3, tool_tip='n_theta')
        self.n_phi_label, self.n_phi_edits = self.createLabelEditObj('Phi Resolution', speed_content, 4, tool_tip='n_phi')
        self.angle_precision_label, self.angle_precision_edits = self.createLabelEditObj('Angle Precision', speed_content, 5, tool_tip='angle_precision')
        self.angle_error_tol_label, self.angle_error_tol_edits = self.createLabelEditObj('Angle Error Tolerance', speed_content, 6, tool_tip='angle_error_tol')

    def initPSOTab(self):
        pso = QWidget()
        pso_content = QGridLayout()
        pso.setLayout(pso_content)
        self.section_widget.addTab(pso, "PSO")

        self.maxiter_label, self.maxiter_edits = self.createLabelEditObj('Max Iterations: ', pso_content, 1, tool_tip='maxiter')
        self.swarmsize_label, self.swarmsize_edits = self.createLabelEditObj('Swarm Size: ', pso_content, 2, tool_tip='swarmsize')
        self.run_times_label, self.run_times_edits = self.createLabelEditObj('Run Times:', pso_content, 3, tool_tip='run_times')
        self.minfunc_label, self.minfunc_edits = self.createLabelEditObj('minfunc:', pso_content, 4, tool_tip='minfunc')
        self.minstep_label, self.minstep_edits = self.createLabelEditObj('minstep:', pso_content, 5, tool_tip='minstep')
        self.phip_label, self.phip_edits = self.createLabelEditObj('phip:', pso_content, 6, tool_tip='phip')
        self.phig_label, self.phig_edits = self.createLabelEditObj('phig:', pso_content, 7, tool_tip='phig')
        self.omega_label, self.omega_edits = self.createLabelEditObj('omega:', pso_content, 8, tool_tip='omega')

        self.pso_debug_label, self.pso_debug_edits = self.createComboBoxObj("PSO Debug: ", pso_content, 9, tool_tip='pso_debug')

    def initGraphingTab(self):

        graphing = QWidget()
        graphing_content = QGridLayout()
        graphing.setLayout(graphing_content)
        self.section_widget.addTab(graphing, "Graphing")

        self.plot_all_stations_label, self.plot_all_stations_edits = self.createComboBoxObj("Plot All Stations: ", graphing_content, 1, items=['True', 'False'], tool_tip='plot_all_stations')
        self.color_toggle_label, self.color_toggle_edits = self.createComboBoxObj("Toggle Color: ", graphing_content, 2, items=['True', 'False'], tool_tip='colortoggle')

        self.dot_tol_label, self.dot_tol_edits = self.createLabelEditObj('Dot Product Tolerance:', graphing_content, 3, tool_tip='dot_tol')
        self.contour_res_label, self.contour_res_edits = self.createLabelEditObj('Contour Resolution:', graphing_content, 4, tool_tip='contour_res')
        self.high_f_label, self.high_f_edits = self.createLabelEditObj('Highlight Fragmentation:', graphing_content, 5, tool_tip='high_f')
        self.high_b_label, self.high_b_edits = self.createLabelEditObj('Highlight Ballistic:', graphing_content, 6, tool_tip='high_b')
        self.rm_stat_label, self.rm_stat_edits = self.createLabelEditObj('Remove Stations:', graphing_content, 7, tool_tip='rm_stat')
        self.img_dim_label, self.img_dim_edits = self.createLabelEditObj('Image Dimensions:', graphing_content, 8, tool_tip='img_dim')

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

        self.atm_atm_file_label, self.atm_atm_file_edits, self.atm_atm_file_buton = self.createFileSearchObj('Atmospheric File: ', profile_tab_content, 1, width=2, h_shift=0)
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

        self.atm_perturbation_file_label, self.atm_perturbation_file_edits, self.atm_perturbation_file_buton = self.createFileSearchObj('Perturbation File: ', profile_tab_content, 3, width=2, h_shift=0)
        self.atm_perturbation_file_buton.clicked.connect(partial(self.fileSearch, ['NetCDF (*.nc)'], self.atm_perturbation_file_edits))

        self.atm_perturb_times_label, self.atm_perturb_times_edits = self.createLabelEditObj('Perturbation Times: ', profile_tab_content, 4, width=2)

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
        self.pick_group_colors = ['r', 'g', 'm', 'w', 'y']

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

        self.pick_list = []

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


        # Init ground map
        #self.m = GroundMap(self.lat_list, self.lon_list, ax=self.map_ax, color_scheme='light')

        for stn in self.stn_list:
            self.m.scatter(stn.position.lat_r, stn.position.lon_r, c='k', s=2)

        current_stn = self.stn_list[self.current_station]
        self.m.scatter(current_stn.position.lat_r, current_stn.position.lon_r, c='r', s=2)

        self.make_picks_map_graph_canvas.draw()
        SolutionGUI.update(self)


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

        mousePoint = self.make_picks_waveform_canvas.vb.mapToView(evt.pos())

        self.make_picks_waveform_canvas.scatterPlot(x=[mousePoint.x()], y=[0], pen='r', update=True)

        pick = Pick(mousePoint.x(), self.stn_list[self.current_station], self.current_station)
        self.pick_list.append(pick)

        if self.setup.debug:
            print("New pick object made: {:} {:} {:}".format(mousePoint.x(), self.stn_list[self.current_station].code, self.current_station))


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
        self.make_picks_waveform_canvas.setXRange(t_arrival-100, t_arrival+100, padding=1)
        # proxy = pg.SignalProxy(self.make_picks_waveform_canvas.scene().sigMouseMoved, rateLimit=60, slot=self.mouseMoved)
        #self.atm_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)


        for pick in self.pick_list:
            if pick.stn_no == self.current_station:
                self.make_picks_waveform_canvas.scatterPlot(x=[pick.time], y=[0], pen='r', update=True)

        SolutionGUI.update(self)
        # Initialize variables
        b_time = 0

        print('####################')
        print("Current Station: {:}".format(stn.name))
        print("Channel: {:}".format(stn.channel))
        # If manual ballistic search is on
        if setup.show_ballistic_waveform:

            # Plot Ballistic Prediction
            b_time = self.arrTimes[0, self.current_station, 0, 0]
            
            # check if nan
            if b_time == b_time:
                self.make_picks_waveform_canvas.plot(x=[b_time]*2, y=[np.min(waveform_data), np.max(waveform_data)], pen='b', label='Ballistic')
                print("Ballistic Arrival: {:.3f} s".format(b_time))
            else:
                print("No Ballistic Arrival")

            for i in range(setup.perturb_times):
                if i >= 1:
                    try:
                        self.make_picks_waveform_canvas.plot(x=[self.arrTimes[i, self.current_station, 0, 0]]*2, \
                         y=[np.min(waveform_data), np.max(waveform_data)], alpha=0.3, pen='b')
                    except:
                        pass
            # Fragmentation Prediction

            # If manual fragmentation search is on
            if setup.show_fragmentation_waveform:

                for i, line in enumerate(setup.fragmentation_point):

                    f_time = self.arrTimes[0, self.current_station, 1, i]
                #     # check if nan
                    if f_time == f_time:
                        # Plot Fragmentation Prediction
                        self.make_picks_waveform_canvas.plot(x=[f_time]*2, y=[np.min(waveform_data), np.max(waveform_data)], pen=self.pick_group_colors[(i+1)%4], label='Fragmentation')
                        
                        #if len(setup.fragmentation_point) > 1:
                            #self.ax_wave.text(f_time, np.min(waveform_data), 'Frag{:}'.format(i+1))
                            #self.wave_ax.text(f_time, np.min(waveform_data) + int(i)/(len(setup.fragmentation_point))*(np.max(waveform_data) - np.min(waveform_data)), '{:.1f} km'.format(line[2]/1000))
                    
                        print('Fragmentation {:} Arrival: {:.3f} s'.format(i+1, f_time))

                    else:
                        print('No Fragmentation {:} Arrival'.format(i+1))

                    for j in range(setup.perturb_times):
                        if j >= 1:
                            try:
                                self.make_picks_waveform_canvas.plot(x=[self.arrTimes[j, self.current_station, 1, i]]*2, y=[np.min(waveform_data),\
                                     np.max(waveform_data)], alpha=0.3,\
                                     pen=self.pick_group_colors[(i+1)%4], zorder=3)
                            except:
                                pass

            # Set the time limits to be within the given window
            # self.wave_ax.set_xlim(time_win_min, time_win_max)

            # self.wave_ax.grid(color='#ADD8E6', linestyle='dashed', linewidth=0.5, alpha=0.5)

            #self.ax_wave.legend()

            # Add text with station label
            # if stn.code in setup.high_f:
            #     self.wave_ax.text(time_win_min, np.max(waveform_data), stn.network + ": " + stn.code \
            #         + "(" + stn.channel + ")" +", {:d} km".format(int(self.source_dists[self.current_station])) , va='top', ha='left', color='g')

            # elif stn.code in setup.high_b:
            #     self.wave_ax.text(time_win_min, np.max(waveform_data), stn.network + ": " + stn.code \
            #         + "(" + stn.channel + ")" + ", {:d} km".format(int(self.source_dists[self.current_station])) , va='top', ha='left', color='b')

            # else:
            #     self.wave_ax.text(time_win_min, np.max(waveform_data), stn.network + ": " + stn.code \
            #         + "(" + stn.channel + ")" + ", {:d} km".format(int(self.source_dists[self.current_station])) , va='top', ha='left', color='k')

            # self.make_picks_top_graphs.removeWidget(self.make_picks_station_graph_canvas)
            # self.make_picks_station_graph_canvas = FigureCanvas(Figure(figsize=(3, 3)))
            # self.make_picks_station_graph_canvas = FigureCanvas(fig)
            # self.make_picks_station_graph_canvas.setSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
            # self.make_picks_top_graphs.addWidget(self.make_picks_station_graph_canvas)    
            # self.make_picks_station_graph_canvas.draw()
            # SolutionGUI.update(self)

    def markStationWaveform(self):
        """ Mark the currently shown waveform in the plot of all waveform. """

        pass

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

        # Reset bandpass filter values to default
        # self.make_picks.set_val(self.bandpass_low_default)
        # self.bandpass_high_slider.set_val(self.bandpass_high_default)

        SolutionGUI.update(self)

    def exportCSV(self, event):
        """ Save picks to a CSV file. """

        dlg = QFileDialog.getSaveFileName(self, 'Save File')


        if '.csv' not in dlg[0]:
            file_name = dlg[0] + '.csv'
        else:
            file_name = dlg[0]

        # Open the output CSV
        with open(os.path.join(file_name), 'w') as f:

            # Write the header
            f.write('Pick group, Network, Code, Lat, Lon, Elev, Pick JD, Pick time, station_number \n')

            # Go through all picks
            for pick in self.pick_list:

                # Calculate Julian date of the pick time
                pick_jd = datetime2JD(self.setup.start_datetime + datetime.timedelta(seconds=pick.time))

                stn = pick.stn

                # Write the CSV entry
                f.write("{:d}, {:s}, {:s}, {:.6f}, {:.6f}, {:.2f}, {:.8f}, {:}, {:}\n".format(0, stn.network, \
                    stn.code, stn.position.lat, stn.position.lon, stn.position.elev, pick_jd, pick.time, pick.stn_no))

        self.errorMessage('Output to CSV!', 0, title='Exported!')


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
        setup.station_picks_file = self.station_picks_edits.text()
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
        self.station_picks_edits.setText(setup.station_picks_file)
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
    #3794