from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *

from functools import partial

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

import pyqtgraph as pg

from supra.GUI.GUITools import *

def addStationsWidgets(obj):
    station_tab = QWidget()

    obj.station_layout = QVBoxLayout()
    station_tab.setLayout(obj.station_layout)
    obj.tab_widget.addTab(station_tab, "Stations")

    obj.station_label = QLabel("Stations:")
    obj.station_layout.addWidget(obj.station_label)

    obj.station_table = QTableWidget()
    obj.station_layout.addWidget(obj.station_table)

    obj.station_button = QPushButton("Get Station Data")
    obj.station_layout.addWidget(obj.station_button)
    obj.station_button.clicked.connect(obj.getStations)

def addPicksReadWidgets(obj):
    picks_read_tab = QWidget()
    picks_read_tab_content = QGridLayout()
    picks_read_tab.setLayout(picks_read_tab_content)

    obj.tab_widget.addTab(picks_read_tab, "Picks Read")

    obj.csv_table = QTableWidget(0, 9)
    picks_read_tab_content.addWidget(obj.csv_table, 1, 1, 1, 4)
    obj.csv_table.setHorizontalHeaderLabels(['Pick Group', 'Network', 'Code', 'Latitude', 'Longitude', 'Elevation', 'Pick JD', 'Pick Time', 'station_number'])
    header = obj.csv_table.horizontalHeader()
    header.setSectionResizeMode(QHeaderView.Stretch)

    obj.csv_table_add = QPushButton("+")
    picks_read_tab_content.addWidget(obj.csv_table_add, 2, 2, 1, 1)
    obj.csv_table_add.clicked.connect(partial(obj.changeRows, obj.csv_table, 1))
    obj.csv_table_add.setToolTip("Add row")

    obj.csv_table_min = QPushButton("-")
    picks_read_tab_content.addWidget(obj.csv_table_min, 2, 1, 1, 1)
    obj.csv_table_min.clicked.connect(partial(obj.changeRows, obj.csv_table, -1))
    obj.csv_table_min.setToolTip("Remove row")

    obj.csv_table_load = QPushButton("Load")
    picks_read_tab_content.addWidget(obj.csv_table_load, 2, 3, 1, 1)
    obj.csv_table_load.clicked.connect(obj.csvLoad)

    obj.csv_table_save = QPushButton("Save")
    picks_read_tab_content.addWidget(obj.csv_table_save, 2, 4, 1, 1)
    obj.csv_table_save.clicked.connect(obj.csvSave)

def addSupraWidgets(obj):

    supra_tab = QWidget()
    obj.master_supra = QGridLayout()
    obj.supra_tab_content = QVBoxLayout()
    obj.plots = QVBoxLayout()

    obj.two_canvas = FigureCanvas(Figure(figsize=(0, 0)))
    obj.two_canvas.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
    obj.plots.addWidget(obj.two_canvas)

    obj.three_canvas = FigureCanvas(Figure(figsize=(0, 0)))
    obj.three_canvas.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
    obj.plots.addWidget(obj.three_canvas)

    obj.results_label = QLabel("Results: ")
    obj.supra_tab_content.addWidget(obj.results_label)

    obj.tableWidget = QTableWidget(0, 0)
    obj.supra_tab_content.addWidget(obj.tableWidget)

    obj.search_button = QPushButton('Search')
    obj.supra_tab_content.addWidget(obj.search_button)
    obj.search_button.clicked.connect(obj.supraSearch)

    obj.master_supra.addLayout(obj.supra_tab_content, 1, 1, 1, 100)
    obj.master_supra.addLayout(obj.plots, 1, 101)

    supra_tab.setLayout(obj.master_supra)
    obj.tab_widget.addTab(supra_tab, "Supracenter Manual Search")

def addSupWidgets(obj):

    sup_tab = QWidget()
    obj.master_sup = QHBoxLayout()
    obj.sup_tab_content = QGridLayout()
    obj.sup_plots = QVBoxLayout()

    obj.sup_two_canvas = FigureCanvas(Figure(figsize=(0, 0)))
    obj.sup_two_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
    obj.sup_plots.addWidget(obj.sup_two_canvas)
    
    obj.sup_three_canvas = FigureCanvas(Figure(figsize=(0, 0)))
    obj.sup_three_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
    obj.sup_plots.addWidget(obj.sup_three_canvas)

    obj.sup_search_button = QPushButton('Search')
    obj.sup_tab_content.addWidget(obj.sup_search_button, 5, 1, 1, 75)
    obj.sup_search_button.clicked.connect(obj.supSearch)

    obj.sup_results_label = QLabel("Results: ")
    obj.sup_tab_content.addWidget(obj.sup_results_label, 3, 1, 1, 1)

    obj.sup_results_table = QTableWidget(0, 0)
    obj.sup_tab_content.addWidget(obj.sup_results_table, 4, 1, 1, 75)

    obj.master_sup.addLayout(obj.sup_tab_content)
    obj.master_sup.addLayout(obj.sup_plots)

    sup_tab.setLayout(obj.master_sup)
    obj.tab_widget.addTab(sup_tab, "Supracenter PSO Search")

def addMakePicksWidgets(obj):
    make_picks_master_tab = QWidget()
    make_picks_master = QVBoxLayout()
    make_picks_master_tab.setLayout(make_picks_master)

    obj.make_picks_top_graphs = QHBoxLayout()
    obj.make_picks_bottom_graphs = QVBoxLayout()
    make_picks_master.addLayout(obj.make_picks_top_graphs)
    make_picks_master.addLayout(obj.make_picks_bottom_graphs)

    obj.make_picks_station_graph_view = pg.GraphicsLayoutWidget()
    obj.make_picks_station_graph_canvas = obj.make_picks_station_graph_view.addPlot()
    obj.make_picks_top_graphs.addWidget(obj.make_picks_station_graph_view)
    obj.make_picks_station_graph_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

    obj.make_picks_map_graph_view = pg.GraphicsLayoutWidget()
    obj.make_picks_map_graph_canvas = obj.make_picks_map_graph_view.addPlot()
    obj.make_picks_top_graphs.addWidget(obj.make_picks_map_graph_view)
    obj.make_picks_map_graph_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

    obj.make_picks_waveform_view = pg.GraphicsLayoutWidget()
    obj.make_picks_waveform_canvas = obj.make_picks_waveform_view.addPlot()
    obj.make_picks_bottom_graphs.addWidget(obj.make_picks_waveform_view)
    obj.make_picks_waveform_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

    make_picks_control_panel = QHBoxLayout()
    obj.make_picks_bottom_graphs.addLayout(make_picks_control_panel)

    make_picks_station_group = QGroupBox("Station Navigation")
    make_picks_control_panel.addWidget(make_picks_station_group)

    station_group_layout = QGridLayout()
    make_picks_station_group.setLayout(station_group_layout)

    obj.make_picks_station_choice = QComboBox()
    station_group_layout.addWidget(obj.make_picks_station_choice, 0, 0, 1, 2)

    obj.make_picks_channel_choice = QComboBox()
    station_group_layout.addWidget(obj.make_picks_channel_choice, 1, 0, 1, 2)
    obj.make_picks_channel_choice.currentIndexChanged.connect(partial(obj.drawWaveform, 1))

    obj.prev_stat = QPushButton('Prev')
    station_group_layout.addWidget(obj.prev_stat, 2, 0, 1, 1)

    obj.next_stat = QPushButton('Next')
    station_group_layout.addWidget(obj.next_stat, 2, 1, 1, 1)

    launch = QPushButton('Load Station Data')
    station_group_layout.addWidget(launch, 3, 0, 1, 2)
    launch.clicked.connect(obj.makePicks)

    make_picks_filter_group = QGroupBox("Waveform Filtering")
    make_picks_control_panel.addWidget(make_picks_filter_group)

    filter_group_layout = QGridLayout()
    make_picks_filter_group.setLayout(filter_group_layout)

    obj.low_bandpass_label = QLabel('Low: 2 Hz')
    filter_group_layout.addWidget(obj.low_bandpass_label, 0, 0, 1, 1)

    obj.low_bandpass_slider = QSlider(Qt.Horizontal)
    filter_group_layout.addWidget(obj.low_bandpass_slider, 0, 1, 1, 1)

    obj.high_bandpass_label = QLabel("High: 8 Hz")
    filter_group_layout.addWidget(obj.high_bandpass_label, 1, 0, 1, 1)

    obj.high_bandpass_slider = QSlider(Qt.Horizontal)
    filter_group_layout.addWidget(obj.high_bandpass_slider, 1, 1, 1, 1)

    obj.low_bandpass_slider.setMinimum(0/obj.bandpass_scale)
    obj.low_bandpass_slider.setMaximum(5/obj.bandpass_scale)
    obj.low_bandpass_slider.setValue(2/obj.bandpass_scale)
    obj.low_bandpass_slider.setTickInterval(0.5)
    obj.low_bandpass_slider.valueChanged.connect(partial(obj.makeValueChange, obj.low_bandpass_label, obj.low_bandpass_slider))

    obj.high_bandpass_slider.setMinimum(3/obj.bandpass_scale)
    obj.high_bandpass_slider.setMaximum(40/obj.bandpass_scale)
    obj.high_bandpass_slider.setValue(8/obj.bandpass_scale)
    obj.high_bandpass_slider.setTickInterval(0.5)
    obj.high_bandpass_slider.valueChanged.connect(partial(obj.makeValueChange, obj.high_bandpass_label, obj.high_bandpass_slider))


    obj.filter_combo_box = QComboBox()
    filter_group_layout.addWidget(obj.filter_combo_box, 2, 0, 1, 2)

    make_picks_picks_group = QGroupBox("Arrival Picks")
    make_picks_control_panel.addWidget(make_picks_picks_group)

    pick_group_layout = QGridLayout()
    make_picks_picks_group.setLayout(pick_group_layout)

    obj.export_to_csv = QPushButton('Export to CSV')
    pick_group_layout.addWidget(obj.export_to_csv)

    obj.export_to_all_times = QPushButton('Export All Times')
    pick_group_layout.addWidget(obj.export_to_all_times)

    obj.tab_widget.addTab(make_picks_master_tab, 'Make Picks')  

def addSeisTrajWidgets(obj):

    seis_tab = QWidget()
    obj.master_seis = QHBoxLayout()
    obj.seis_tab_input = QVBoxLayout()
    obj.seis_tab_output = QVBoxLayout()

    seis_tab.setLayout(obj.master_seis)
    obj.tab_widget.addTab(seis_tab, "Seismic Trajectory")

    obj.master_seis.addLayout(obj.seis_tab_input)
    obj.master_seis.addLayout(obj.seis_tab_output)

    table_group = QGridLayout()
    obj.seis_tab_input.addLayout(table_group)

    tab_layout = QGridLayout()
    obj.seis_tab_input.addLayout(tab_layout)

    obj.seis_search = QPushButton('Search')
    tab_layout.addWidget(obj.seis_search, 1, 1, 1, 4)
    obj.seis_search.clicked.connect(obj.seisSearch)
    
    obj.seis_table = QTableWidget()
    table_group.addWidget(obj.seis_table, 1, 1, 1, 4)

    obj.seis_resids = QTableWidget()
    table_group.addWidget(obj.seis_resids, 2, 1, 1, 4)

    obj.seis_three_canvas = FigureCanvas(Figure(figsize=(5, 5)))
    obj.seis_three_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
    obj.seis_tab_output.addWidget(obj.seis_three_canvas)

    two_graphs = QGridLayout()
    obj.seis_tab_output.addLayout(two_graphs)

    obj.seis_two_lat_view = pg.GraphicsLayoutWidget()
    obj.seis_two_lat_canvas = obj.seis_two_lat_view.addPlot()
    two_graphs.addWidget(obj.seis_two_lat_view, 1, 1)
    obj.seis_two_lat_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)
    
    obj.seis_two_time_view = pg.GraphicsLayoutWidget()
    obj.seis_two_time_canvas = obj.seis_two_time_view.addPlot()
    two_graphs.addWidget(obj.seis_two_time_view, 1, 2)
    obj.seis_two_time_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)
    
    obj.seis_two_angle_view = pg.GraphicsLayoutWidget()
    obj.seis_two_angle_canvas = obj.seis_two_angle_view.addPlot()
    two_graphs.addWidget(obj.seis_two_angle_view, 2, 1)
    obj.seis_two_angle_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)
    
    obj.seis_two_plot_view = pg.GraphicsLayoutWidget()
    obj.seis_two_plot_canvas = obj.seis_two_plot_view.addPlot()
    two_graphs.addWidget(obj.seis_two_plot_view, 2, 2)
    obj.seis_two_plot_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

def addFetchATMWidgets(obj):
    fetch = QWidget()
    fetch_master = QHBoxLayout()
    fetch_content = QGridLayout()
    fetch_plots = QVBoxLayout()
    fetch_master.addLayout(fetch_plots)
    fetch_master.addLayout(fetch_content)
    fetch.setLayout(fetch_master)
    obj.tab_widget.addTab(fetch, "Fetch Atmosphere")

    obj.fatm_view = pg.GraphicsLayoutWidget()
    obj.fatm_canvas = obj.fatm_view.addPlot()
    fetch_plots.addWidget(obj.fatm_view)
    obj.fatm_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

    obj.fatm_variable_combo = QComboBox()
    fetch_plots.addWidget(obj.fatm_variable_combo)
    obj.fatm_variable_combo.currentTextChanged.connect(obj.fatmPlot)

    obj.fatm_name_label, obj.fatm_name_edits = createLabelEditObj('Name:', fetch_content, 1)

    obj.fatm_button = QPushButton('Browse')
    fetch_content.addWidget(obj.fatm_button, 1, 3)
    obj.fatm_button.clicked.connect(partial(obj.folderSearch, obj.fatm_name_edits))

    obj.fatm_datetime_label = QLabel("Time of Profile:")
    obj.fatm_datetime_edits = QDateTimeEdit()
    fetch_content.addWidget(obj.fatm_datetime_label, 2, 1)
    fetch_content.addWidget(obj.fatm_datetime_edits, 2, 2)
    obj.fatm_datetime_edits.setCalendarPopup(True)

    obj.fatm_variable_group = QGroupBox("Variables")
    fetch_content.addWidget(obj.fatm_variable_group, 4, 1, 1, 2)

    #############################
    group_box = QGridLayout()
    obj.fatm_variable_group.setLayout(group_box)

    obj.fatm_temp = QCheckBox('Temperature')
    group_box.addWidget(obj.fatm_temp, 1, 1)

    obj.fatm_u_wind = QCheckBox('U Wind')
    group_box.addWidget(obj.fatm_u_wind, 1, 2)

    obj.fatm_v_wind = QCheckBox('V Wind')
    group_box.addWidget(obj.fatm_v_wind, 1, 3)
    ##############################

    obj.fatm_fetch = QPushButton("Download")
    fetch_content.addWidget(obj.fatm_fetch, 5, 1, 1, 2)
    obj.fatm_fetch.clicked.connect(partial(obj.fatmFetch, True))

    obj.fatm_open = QPushButton("Open")
    fetch_content.addWidget(obj.fatm_open, 6, 1, 1, 2)
    obj.fatm_open.clicked.connect(partial(obj.fatmFetch, False))

    obj.fatm_print = QPushButton("Print")
    fetch_content.addWidget(obj.fatm_print, 9, 1, 1, 2)
    obj.fatm_print.clicked.connect(obj.fatmPrint)

    obj.fatm_lat_label = QLabel("Latitude: 0")
    obj.fatm_lat_slide = QSlider(Qt.Horizontal)
    fetch_content.addWidget(obj.fatm_lat_label, 7, 1)
    fetch_content.addWidget(obj.fatm_lat_slide, 7, 2)
    obj.fatm_lat_slide.setMinimum(-90/obj.slider_scale)
    obj.fatm_lat_slide.setMaximum(90/obj.slider_scale)
    obj.fatm_lat_slide.setValue(0)
    obj.fatm_lat_slide.setTickInterval(0.5)
    obj.fatm_lat_slide.valueChanged.connect(partial(obj.fatmValueChange, obj.fatm_lat_label, obj.fatm_lat_slide))

    obj.fatm_lon_label = QLabel("Longitude: 0")
    obj.fatm_lon_slide = QSlider(Qt.Horizontal)
    fetch_content.addWidget(obj.fatm_lon_label, 8, 1)
    fetch_content.addWidget(obj.fatm_lon_slide, 8, 2)
    obj.fatm_lon_slide.setMinimum(-180/obj.slider_scale)
    obj.fatm_lon_slide.setMaximum(180/obj.slider_scale)
    obj.fatm_lon_slide.setValue(0)
    obj.fatm_lon_slide.setTickInterval(0.5)
    obj.fatm_lon_slide.valueChanged.connect(partial(obj.fatmValueChange, obj.fatm_lon_label, obj.fatm_lon_slide))

def addProfileWidgets(obj):
    profile_tab = QWidget()
    profile_master = QHBoxLayout()
    profile_tab_content = QGridLayout()
    profile_tab_content_graph = QVBoxLayout()
    profile_tab.setLayout(profile_master)
    profile_master.addLayout(profile_tab_content_graph)
    profile_master.addLayout(profile_tab_content)

    obj.atm_lat_label = QLabel("Latitude: 0")
    obj.atm_lat_slide = QSlider(Qt.Horizontal)
    profile_tab_content.addWidget(obj.atm_lat_label, 6, 1)
    profile_tab_content.addWidget(obj.atm_lat_slide, 6, 2, 1, 3)
    obj.atm_lat_slide.setMinimum(-90/obj.slider_scale)
    obj.atm_lat_slide.setMaximum(90/obj.slider_scale)
    obj.atm_lat_slide.setValue(0)
    obj.atm_lat_slide.setTickInterval(0.5)
    obj.atm_lat_slide.valueChanged.connect(partial(obj.atmValueChange, obj.atm_lat_label, obj.atm_lat_slide))

    obj.atm_lon_label = QLabel("Longitude: 0")
    obj.atm_lon_slide = QSlider(Qt.Horizontal)
    profile_tab_content.addWidget(obj.atm_lon_label, 7, 1)
    profile_tab_content.addWidget(obj.atm_lon_slide, 7, 2, 1, 3)
    obj.atm_lon_slide.setMinimum(-180/obj.slider_scale)
    obj.atm_lon_slide.setMaximum(180/obj.slider_scale)
    obj.atm_lon_slide.setValue(0)
    obj.atm_lon_slide.setTickInterval(0.5)
    obj.atm_lon_slide.valueChanged.connect(partial(obj.atmValueChange, obj.atm_lon_label, obj.atm_lon_slide))

    obj.atm_view = pg.GraphicsLayoutWidget()
    obj.atm_canvas = obj.atm_view.addPlot()

    profile_tab_content_graph.addWidget(obj.atm_view)
    obj.atm_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

    obj.atm_T_button = QPushButton('Temperature')
    profile_tab_content_graph.addWidget(obj.atm_T_button)
    obj.atm_T_button.clicked.connect(partial(obj.atmPlotProfile, obj.atm_lat_slide.value()*obj.slider_scale, obj.atm_lon_slide.value()*obj.slider_scale, var_typ='t'))

    obj.atm_mag_button = QPushButton('Wind Magnitude')
    profile_tab_content_graph.addWidget(obj.atm_mag_button)
    obj.atm_mag_button.clicked.connect(partial(obj.atmPlotProfile, obj.atm_lat_slide.value()*obj.slider_scale, obj.atm_lon_slide.value()*obj.slider_scale, var_typ='m'))

    obj.atm_dir_button = QPushButton('Wind Direction')
    profile_tab_content_graph.addWidget(obj.atm_dir_button)
    obj.atm_dir_button.clicked.connect(partial(obj.atmPlotProfile, obj.atm_lat_slide.value()*obj.slider_scale, obj.atm_lon_slide.value()*obj.slider_scale, var_typ='d'))

    obj.tab_widget.addTab(profile_tab, "Atmospheric Profile")

def addRayTracerWidgets(obj):
    ray_tab = QWidget()

    obj.master_ray = QVBoxLayout()
    obj.ray_graphs = QHBoxLayout()
    obj.ray_control = QGridLayout()

    obj.master_ray.addLayout(obj.ray_graphs)
    obj.master_ray.addLayout(obj.ray_control)

    ray_tab.setLayout(obj.master_ray)
    obj.tab_widget.addTab(ray_tab, "Ray Tracer")

    obj.ray_view = pg.GraphicsLayoutWidget()
    obj.ray_canvas = obj.ray_view.addPlot()
    obj.ray_graphs.addWidget(obj.ray_view)
    obj.ray_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)
    obj.ray_canvas.setLabel('bottom', "Longitude", units='deg E')
    obj.ray_canvas.setLabel('left', "Latitude", units='deg N')

    obj.ray_line_canvas = FigureCanvas(Figure(figsize=(0, 0)))
    obj.ray_line_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
    obj.ray_graphs.addWidget(obj.ray_line_canvas)

    obj.ray_label = QLabel('Starting Point')
    obj.ray_control.addWidget(obj.ray_label, 1, 0)
    
    obj.ray_label2 = QLabel('Ending Point')
    obj.ray_control.addWidget(obj.ray_label2, 3, 0, 1, 1)

    obj.ray_height_label, obj.ray_height_edits = createLabelEditObj("Height", obj.ray_control, 1, width=3)
    
    obj.ray_button = QPushButton('Solve for Lat/Lon')
    obj.ray_control.addWidget(obj.ray_button, 7, 3, 4, 1)
    obj.ray_button.clicked.connect(obj.trajSolver)

    obj.ray_button = QPushButton('Load Data')
    obj.ray_control.addWidget(obj.ray_button, 7, 1, 1, 1)
    obj.ray_button.clicked.connect(obj.loadRayGraph)

    obj.ray_lat_label, obj.ray_lat_edits = createLabelEditObj("Lat", obj.ray_control, 2)
    obj.ray_lon_label, obj.ray_lon_edits = createLabelEditObj("Lon", obj.ray_control, 2, h_shift=2)

    obj.ray_pick_label = QLabel('')
    obj.ray_control.addWidget(obj.ray_pick_label, 3, 1, 1, 5)

    obj.ray_enable_windfield = QCheckBox('Enable Wind Field')
    obj.ray_control.addWidget(obj.ray_enable_windfield, 1, 5)
    obj.ray_enable_windfield.stateChanged.connect(obj.trajSolver)

    obj.ray_enable_perts = QCheckBox('Enable Perturbations')
    obj.ray_control.addWidget(obj.ray_enable_perts, 2, 5)
    obj.ray_enable_perts.stateChanged.connect(obj.trajSolver)

    obj.ray_canvas.scene().sigMouseClicked.connect(obj.rayMouseClicked)