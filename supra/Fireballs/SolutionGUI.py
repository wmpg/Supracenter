
import os
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time
import datetime

from PyQt5.QtWidgets import *
import sys
from PyQt5.QtGui import *
from PyQt5.QtCore import *

from functools import partial

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from mpl_toolkits.mplot3d import Axes3D

from supra.Supracenter.stationDat import convStationDat
from supra.Supracenter.psoSearch import psoSearch
from supra.Fireballs.Program import position, configRead, configWrite
from supra.Fireballs.SeismicTrajectory import Constants
from supra.Supracenter.angleConv import loc2Geo

class SolutionGUI(QMainWindow):
    def __init__(self):
        super().__init__()

        self._main = QWidget()
        self.setCentralWidget(self._main)
        layout = QGridLayout(self._main)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        self.tab_widget = QTabWidget()
        self.tab_widget.blockSignals(True)
        self.tab_widget.currentChanged.connect(self.tabChange)
        # self.tab_widget.blockSignals(False)

        self.addIniWidgets()
        self.addSupraWidgets()

        self.tab_widget.blockSignals(False)
        layout.addWidget(self.tab_widget, 1, 1)

        stylesheet = """ 
        QTabWidget>QWidget>QWidget{background: gray;}
        QLabel{color: white;}
        """

        self.setStyleSheet(stylesheet)

    def toolTime(self, var):

        tool_time_dict = {}

        with open('supra/Fireballs/tool_tips.csv') as f:
            for line in f:
                line = line.split(':')
                tool_time_dict[line[0]] = line[1]

        return tool_time_dict[var]

    def errorMessage(self, message, level, info='', title='Yikes!', detail='Unknown Error'):
       
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

    def addSupraWidgets(self):

        #self.tab_widget.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

        supra_tab = QWidget()
        self.master_supra = QHBoxLayout()
        self.supra_tab_content = QGridLayout()
        self.plots = QVBoxLayout()


        self.three_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.three_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.plots.addWidget(self.three_canvas)
        #self.three_ax = self.three_canvas.figure.subplots()

        self.two_canvas = FigureCanvas(Figure(figsize=(0, 0)))
        self.two_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        self.plots.addWidget(self.two_canvas)
        #self.two_ax = self.two_canvas.figure.subplots()

        self.search_button = QPushButton('Search')
        self.supra_tab_content.addWidget(self.search_button, 10, 3, 1, 3)
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

        self.results_label = QLabel("Results: ")
        self.supra_tab_content.addWidget(self.results_label, 8, 1, 1, 1)

        self.tableWidget = QTableWidget(0, 0)
        self.supra_tab_content.addWidget(self.tableWidget, 9, 1, 1, 10)

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
        restriction_content.addWidget(self.restricted_time_label, 3, 1, 1, 2)
        restriction_content.addWidget(self.restricted_time_edits, 3, 3, 1, 2)
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

    def supraSearch(self):
        supra_pos = [float(self.lat_edit.text()),
                     float(self.lon_edit.text()),
                     float(self.elev_edit.text()),
                     float(self.time_edit.text())]
        picks_file = self.picks_file_edit.text()
        atm_file = self.atmospheric_file_edit.text()

        s_info, s_name, weights, ref_pos = convStationDat(picks_file)
        ref_pos = position(ref_pos[0], ref_pos[1], ref_pos[2])

        dataset = np.array([[    0.0, 310, 0.0, 0.0],
                             [10000.0, 310, 0.0, 0.0],
                             [99999.0, 310, 0.0, 0.0]])

        consts = Constants()

        class Config:

            def __init__(self):
                pass

        setup = Config()

        setup.manual_fragmentation_search = supra_pos
        setup.search = [0, 1, 0, 1, 0, 1]
        setup.ref_pos = ref_pos
        setup.start_datetime = datetime.datetime.strptime(self.ref_edit.text(), "%Y-%m-%d %H:%M:%S.%f")
        setup.restricted_time = 0
        setup.output_name = "/home/luke/Desktop/"
        setup.grid_size = 0.5
        setup.weather_type = 'none'
        setup.enable_winds = True
        setup.n_theta = 37
        setup.n_phi = 73
        setup.angle_precision = 1e-8
        setup.angle_error_tol = 1e-8
        setup.fit_type = 1

        results = psoSearch(s_info, weights, s_name, setup, dataset, consts)

        n_stations = len(s_info)
        xstn = s_info[0:n_stations, 0:3]
            
        self.scatterPlot(setup, results, n_stations, xstn, s_name, dataset)

        self.residPlot(results, s_name, xstn, setup.output_name, n_stations)

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

        setup.lat_centre = self.lat_centre_edits.text()
        setup.lon_centre = self.lon_centre_edits.text()
        setup.deg_radius = self.deg_radius_edits.text()
        setup.start_datetime = self.start_datetime_edits.dateTime().toPyDateTime()
        setup.end_datetime = self.end_datetime_edits.dateTime().toPyDateTime()
        setup.v_sound = self.v_sound_edits.text()

        setup.t0 = self.t0_edits.text()
        setup.v = self.v_edits.text()
        setup.azim = self.azim_edits.text()
        setup.zangle = self.zangle_edits.text()
        setup.lat_i = self.lat_i_edits.text()
        setup.lon_i = self.lon_i_edits.text()
        setup.elev_i = self.elev_i_edits.text()
        setup.lat_f = self.lat_f_edits.text()
        setup.lon_f = self.lon_f_edits.text()
        setup.elev_f = self.elev_f_edits.text()
        setup.show_ballistic_waveform = self.show_ballistic_waveform_edits.currentText()

        setup.manual_fragmentation_search = []
        setup.fragmentation_point = self.fromTable(self.fragmentation_point)
        setup.show_fragmentation_waveform = self.show_fragmentation_waveform_edits.currentText()

        setup.manual_fragmentation_search.append(self.lat_frag_edits.text())
        setup.manual_fragmentation_search.append(self.lon_frag_edits.text())
        setup.manual_fragmentation_search.append(self.elev_frag_edits.text())
        setup.manual_fragmentation_search.append(self.time_frag_edits.text())

        setup.v_fixed = self.v_fixed_edits.text()
        setup.max_error = self.max_error_edits.text()
        setup.restricted_time = self.restricted_time_edits.dateTime().toPyDateTime()
        setup.traj_tol = self.traj_tol_edits.text()
        setup.restrict_to_trajectory = self.restrict_to_trajectory_edits.currentText()

        setup.search_area = []
        setup.search_height = []
        setup.azimuth_min = self.azimuth_min_edits.text()
        setup.azimuth_max = self.azimuth_max_edits.text()
        setup.zangle_min = self.zangle_min_edits.text()
        setup.zangle_max = self.zangle_max_edits.text()
        setup.x_min = self.x_min_edits.text()
        setup.x_max = self.x_max_edits.text()
        setup.y_min = self.y_min_edits.text()
        setup.y_max = self.y_max_edits.text()
        setup.t_min = self.t_min_edits.text()
        setup.t_max = self.t_max_edits.text()
        setup.v_min = self.v_min_edits.text()
        setup.v_max = self.v_max_edits.text()
        setup.weight_distance_min = self.weight_distance_min_edits.text()
        setup.weight_distance_max = self.weight_distance_max_edits.text()
        setup.min_time = self.search_time_min_edits.text()
        setup.max_time = self.search_time_max_edits.text()
        setup.search_area.append(self.search_lat_min_edits.text())
        setup.search_area.append(self.search_lat_max_edits.text())
        setup.search_area.append(self.search_lon_min_edits.text())
        setup.search_area.append(self.search_lon_max_edits.text())
        setup.search_height.append(self.search_elev_min_edits.text())
        setup.search_height.append(self.search_elev_max_edits.text())

        setup.enable_winds = self.enable_winds_edits.currentText()
        setup.weather_type = self.weather_type_edits.currentText()
        setup.grid_size = self.grid_size_edits.text()

        setup.perturb_times = self.perturb_times_edits.text()
        setup.observe_frag_no = self.frag_no_edits.text()
        setup.perturb = self.perturb_edits.currentText()
        setup.perturb_method = self.perturb_method_edits.currentText()

        setup.fast_ballistic = self.fast_ballistic_edits.currentText()
        setup.fit_type = self.fit_type_edits.text()
        setup.n_theta = self.n_theta_edits.text()
        setup.n_phi = self.n_phi_edits.text()
        setup.angle_precision = self.angle_precision_edits.text()
        setup.angle_error_tol = self.angle_error_tol_edits.text()

        setup.maxiter = self.maxiter_edits.text()
        setup.swarmsize = self.swarmsize_edits.text()
        setup.run_times = self.run_times_edits.text()
        setup.minfunc = self.minfunc_edits.text()
        setup.minstep = self.minstep_edits.text()
        setup.phip = self.phip_edits.text()
        setup.phig = self.phig_edits.text()
        setup.omega = self.omega_edits.text()
        setup.pso_debug = self.pso_debug_edits.currentText()

        setup.plot_all_stations = self.plot_all_stations_edits.currentText()
        setup.colortoggle = self.color_toggle_edits.currentText()
        setup.dot_tol = self.dot_tol_edits.text()
        setup.contour_res = self.contour_res_edits.text()
        setup.high_f = self.high_f_edits.text()
        setup.high_b = self.high_b_edits.text()
        setup.rm_stat = self.rm_stat_edits.text()
        setup.img_dim = self.img_dim_edits.text()
        setup.reported_points = self.fromTable(self.reported_points)

        setup.stations = self.fromTable(self.extra_point)

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

    def scatterPlot(self, setup, results, n_stations, xstn, s_name, dataset):

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

        # colorbars
        b = plt.colorbar(res, ax=ax)
        b.set_label("Station Residuals (s)")

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
        ax.mouse_init()
        SolutionGUI.update(self)

    def residPlot(self, results_arr, s_name, xstn, output_name, n_stations):
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
        SolutionGUI.update(self)

    def tabChange(self):

        setup = self.saveINI(False)

        # manual Supracenter
        self.lat_edit.setText(setup.manual_fragmentation_search[0])
        self.lon_edit.setText(setup.manual_fragmentation_search[1])
        self.elev_edit.setText(setup.manual_fragmentation_search[2])
        self.time_edit.setText(setup.manual_fragmentation_search[3])

        self.ref_edit.setDateTime(setup.start_datetime)
        self.picks_file_edit.setText(setup.station_picks_file)
        self.atmospheric_file_edit.setText(setup.sounding_file)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    gui = SolutionGUI()

    w = 1280; h = 1024
    gui.resize(w, h)
    gui.show()

    app.exec_()