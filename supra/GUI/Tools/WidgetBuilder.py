import os

from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5.QtWebEngineWidgets import QWebEngineView as QWebView

from functools import partial

from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl

from supra.GUI.Tools.GUITools import *
from supra.GUI.Tools.CustomWidgets import *
from supra.Files.SaveLoad import save, load
from supra.GUI.Tools.Theme import theme

from supra.GUI.Tabs.StationInfo import *
from supra.GUI.Tabs.SourcesInfo import *


def initTabs(obj):
    
    addDocWidgets(obj)
    addStationsWidgets(obj)
    addSourceWidgets(obj)
    addPicksReadWidgets(obj)
    addMakePicksWidgets(obj)
    addSupraWidgets(obj)
    addSeisTrajWidgets(obj)
    addFetchATMWidgets(obj)
    # addProfileWidgets(obj)
    # addRayTracerWidgets(obj)
    # addRefineWidgets(obj)


def initMenuBar(obj, layout):
    menu_bar = obj.menuBar() 
    layout.addWidget(menu_bar, 0, 1)
    file_menu = menu_bar.addMenu('&File')
    about_menu = menu_bar.addMenu('&About')
    view_menu = menu_bar.addMenu('&View')
    setting_menu = menu_bar.addMenu('&Settings')
    tools_menu = menu_bar.addMenu('&Tools')

    file_qsave = QAction("Quick Save", obj)
    file_qsave.setShortcut('Ctrl+S')
    file_qsave.setStatusTip('Saves setup file')
    file_qsave.triggered.connect(partial(save, obj))
    file_menu.addAction(file_qsave)

    file_save = QAction("Save", obj)
    file_save.setShortcut('Ctrl+Shift+S')
    file_save.setStatusTip('Saves setup file')
    file_save.triggered.connect(partial(save, obj))
    file_menu.addAction(file_save)

    file_rep = QAction("Generate Report", obj)
    file_rep.setStatusTip('Opens report dialog')
    file_rep.triggered.connect(obj.genReport)
    file_menu.addAction(file_rep)    

    file_load = QAction("Load", obj)
    file_load.setShortcut('Ctrl+L')
    file_load.setStatusTip('Loads setup file')
    file_load.triggered.connect(partial(load, obj))
    file_menu.addAction(file_load)

    file_exit = QAction("Exit", obj)
    file_exit.setShortcut('Ctrl+Q')
    file_exit.setStatusTip('Exit application')
    file_exit.triggered.connect(obj.quitApp)
    file_menu.addAction(file_exit)

    about_github = QAction("GitHub", obj)
    about_github.triggered.connect(obj.openGit)
    about_menu.addAction(about_github)

    about_docs = QAction("Documentation", obj)
    about_docs.setShortcut('Ctrl+D')
    about_docs.triggered.connect(obj.openDocs)
    about_menu.addAction(about_docs)

    view_vartools = QAction("Show/Hide Toolbar", obj)
    view_vartools.setShortcut('V')
    view_vartools.setStatusTip('Toggle if the variable toolbar is visible')
    view_vartools.triggered.connect(obj.viewToolbar)
    view_menu.addAction(view_vartools)

    view_fullscreen = QAction("Fullscreen", obj)
    view_fullscreen.setShortcut('F11')
    view_fullscreen.setStatusTip('Toggles fullscreen')
    view_fullscreen.triggered.connect(obj.viewFullscreen)
    view_menu.addAction(view_fullscreen)

    preferences_menu = QAction("Preferences", obj)
    preferences_menu.setStatusTip('Opens preferences dialog')
    preferences_menu.triggered.connect(obj.preferencesDialog)
    setting_menu.addAction(preferences_menu)

    stndownload_menu = QAction("Station Download Sources", obj)
    stndownload_menu.triggered.connect(obj.stndownloadDialog)
    setting_menu.addAction(stndownload_menu)

    traj_interp = QAction("Trajectory to Points Wizard", obj)
    traj_interp.triggered.connect(obj.trajInterpDialog)
    tools_menu.addAction(traj_interp)

def initMainGUI(obj):

    obj._main = QWidget()
    obj.setCentralWidget(obj._main)
    layout = QGridLayout(obj._main)
    initMenuBar(obj, layout)

    obj.tab_widget = QTabWidget()
    obj.tab_widget.blockSignals(True)

    obj.tab_widget.blockSignals(False)
    layout.addWidget(obj.tab_widget, 1, 1)

    obj.doc_file = os.path.join('supra', 'Docs', 'index.html')
       
    obj.contour_data = None

    obj.slider_scale = 0.25
    obj.bandpass_scale = 0.1
    
    obj.inverted = False
    obj.showtitled = False

    obj.ini_dock = QDockWidget("Variables", obj)
    obj.addDockWidget(Qt.LeftDockWidgetArea, obj.ini_dock)
    obj.ini_dock.setFeatures(QtGui.QDockWidget.DockWidgetFloatable | QtGui.QDockWidget.DockWidgetMovable)
    
    obj.group_no = 0
    obj.position = []

    obj.contour_data_squares = None

    obj.var_typ = 't'

    initTabs(obj)

    obj.doc_view.load(QUrl.fromLocalFile(os.path.abspath(obj.doc_file)))

    pg.setConfigOptions(antialias=True)
    obj.ray_pick = pg.ScatterPlotItem()
    obj.ray_pick_traj = pg.ScatterPlotItem()
    obj.ray_pick_point = [0, 0, 0]
    obj.ctrl_pressed = False
    obj.shift_pressed = False
    obj.alt_pressed = False

    

def initMainGUICosmetic(obj):

    obj.setWindowTitle('Bolide Acoustic Modelling')
    app_icon = QtGui.QIcon()
    app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM.png'), QtCore.QSize(16, 16))
    obj.setWindowIcon(app_icon)

    p = obj.palette()
    p.setColor(obj.backgroundRole(), Qt.black)
    obj.setPalette(p)

    obj.colors = [(0, 255, 26), (3, 252, 219), (252, 3, 3), (223, 252, 3), (255, 133, 3),
                  (149, 0, 255), (76, 128, 4), (82, 27, 27), (101, 128, 125), (255, 230, 249)]


    theme(obj)

def addStationsWidgets(obj):
    station_tab = QWidget()

    obj.station_layout = QVBoxLayout()
    station_tab.setLayout(obj.station_layout)
    obj.tab_widget.addTab(station_tab, "Stations")

    obj.station_label = QLabel("Stations:")
    obj.station_layout.addWidget(obj.station_label)

    obj.station_table = QScrollArea()
    obj.station_layout.addWidget(obj.station_table)
    obj.station_table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    obj.station_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

    obj.station_table.setWidgetResizable(True)

    container = QWidget()
    container.setStyleSheet("""
        QWidget {
            background-color: rgb(0, 0, 0);
            }
        """)
    obj.station_table.setWidget(container)
    obj.station_table_layout = QVBoxLayout(container)
    obj.station_table_layout.setSpacing(10)
    # obj.station_table.setLayout(obj.station_table_layout)


    obj.station_control_layout = QGridLayout()

    obj.station_button = createButton("Download Station Data", obj.station_control_layout, 1, 1, \
                                                        getStations, args=[obj])

    obj.station_button = createButton("Save Station Data", obj.station_control_layout, 1, 3, \
                                                        saveStations, args=[obj])  

    obj.station_button = createButton("Load Station Data", obj.station_control_layout, 1, 2, \
                                                        loadStations, args=[obj])

    obj.station_button = createButton("Clear Station Data", obj.station_control_layout, 2, 1, \
                                                        clearStationWidgets, args=[obj])

    obj.station_button = createButton("Add Station", obj.station_control_layout, 2, 2, \
                                                        addStation, args=[obj])

    obj.station_button = createButton("Refresh Station", obj.station_control_layout, 2, 3, \
                                                        refreshStation, args=[obj])

    obj.station_button = createButton("Output Station Metadata", obj.station_control_layout, 3, 1, \
                                                        outputStation, args=[obj])

    obj.station_button = createButton("Count Stations", obj.station_control_layout, 3, 2, \
                                                        countStation, args=[obj])    

    obj.station_layout.addLayout(obj.station_control_layout)


def addSourceWidgets(obj):
    sources_tab = QWidget()

    obj.sources_layout = QVBoxLayout()
    sources_tab.setLayout(obj.sources_layout)
    obj.tab_widget.addTab(sources_tab, "Sources")

    obj.sources_label = QLabel("Sources:")
    obj.sources_layout.addWidget(obj.sources_label)

    obj.sources_table = QScrollArea()
    obj.sources_layout.addWidget(obj.sources_table)
    obj.sources_table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOn)
    obj.sources_table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

    obj.sources_table.setWidgetResizable(True)

    container = QWidget()
    container.setStyleSheet("""
        QWidget {
            background-color: rgb(0, 0, 0);
            }
        """)
    obj.sources_table.setWidget(container)
    obj.sources_table_layout = QVBoxLayout(container)
    obj.sources_table_layout.setSpacing(10)
    # obj.station_table.setLayout(obj.station_table_layout)


    obj.sources_control_layout = QGridLayout()


    obj.sources_save_button = createButton("Load Sources", obj.sources_control_layout, 1, 2, \
                                                loadSource, args=[obj])

    obj.sources_add_button = createButton("Add Source", obj.sources_control_layout, 2, 2, \
                                                addSource, args=[obj])

    obj.sources_save_button = createButton("Save Sources", obj.sources_control_layout, 1, 1, \
                                                saveSource, args=[obj])

    obj.sources_del_button = createButton("Delete Sources", obj.sources_control_layout, 2, 1, \
                                                delSource, args=[obj])


    obj.sources_layout.addLayout(obj.sources_control_layout)

def addPicksReadWidgets(obj):
    picks_read_tab = QWidget()
    picks_read_tab_content = QGridLayout()
    picks_read_tab.setLayout(picks_read_tab_content)

    obj.tab_widget.addTab(picks_read_tab, "Picks Read")

    obj.csv_table = QTableWidget()
    picks_read_tab_content.addWidget(obj.csv_table, 1, 1, 1, 4)

    obj.csv_table_add = QPushButton("+")
    picks_read_tab_content.addWidget(obj.csv_table_add, 2, 2, 1, 1)
    obj.csv_table_add.clicked.connect(partial(changeRows, obj.csv_table, 1))
    obj.csv_table_add.setToolTip("Add row")

    obj.csv_table_min = QPushButton("-")
    picks_read_tab_content.addWidget(obj.csv_table_min, 2, 1, 1, 1)
    obj.csv_table_min.clicked.connect(partial(changeRows, obj.csv_table, -1))
    obj.csv_table_min.setToolTip("Remove row")

    obj.csv_table_load = QPushButton("Load")
    picks_read_tab_content.addWidget(obj.csv_table_load, 2, 3, 1, 1)
    obj.csv_table_load.clicked.connect(partial(obj.csvLoad, obj.csv_table))

    obj.csv_table_save = QPushButton("Save")
    picks_read_tab_content.addWidget(obj.csv_table_save, 2, 4, 1, 1)
    obj.csv_table_save.clicked.connect(partial(obj.csvSave, obj.csv_table))

def addSupraWidgets(obj):

    supra_tab = QWidget()
    obj.master_supra = QGridLayout()
    obj.supra_tab_content = QVBoxLayout()
    obj.plots = QVBoxLayout()

    obj.residfig = FigureCanvas(Figure())
    obj.residfig.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
    obj.plots.addWidget(obj.residfig)

    obj.suprafig = FigureCanvas(Figure())
    obj.suprafig.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
    obj.plots.addWidget(obj.suprafig)

    obj.results_label = QLabel("Results: ")
    obj.supra_tab_content.addWidget(obj.results_label)

    obj.supra_res_table = QTableWidget(0, 0)
    obj.supra_tab_content.addWidget(obj.supra_res_table)

    obj.search_button = QPushButton('Manual Search')
    obj.supra_tab_content.addWidget(obj.search_button)
    obj.search_button.clicked.connect(partial(obj.supSearchSetup, True))

    obj.search_p_button = QPushButton('PSO Search')
    obj.supra_tab_content.addWidget(obj.search_p_button)
    obj.search_p_button.clicked.connect(partial(obj.supSearchSetup, False))

    obj.master_supra.addLayout(obj.supra_tab_content, 1, 1, 1, 100)
    obj.master_supra.addLayout(obj.plots, 1, 101)

    supra_tab.setLayout(obj.master_supra)
    obj.tab_widget.addTab(supra_tab, "Supracenter Search")

def addRefineWidgets(obj):
    refine_tab = QWidget()
    obj.master_ref = QHBoxLayout()
    refine_tab.setLayout(obj.master_ref)

    obj.ref_control = QVBoxLayout()
    obj.ref_wave_layout = QVBoxLayout()
    obj.ref_wave_control = QHBoxLayout()

    obj.ref_waveforms = QScrollArea()

    obj.ref_traj = QGridLayout()
    obj.ref_frags = QVBoxLayout()
    obj.ref_buttons = QHBoxLayout()

    trajLabel = QLabel("Trajectory")
    obj.ref_traj.addWidget(trajLabel, 0, 0)

    latlabel = QLabel('Lat: ')
    obj.ref_traj.addWidget(latlabel, 1, 0)
    obj.latedit = QLineEdit('')
    obj.ref_traj.addWidget(obj.latedit, 1, 1)

    lonlabel = QLabel('Lon: ')
    obj.ref_traj.addWidget(lonlabel, 1, 2)
    obj.lonedit = QLineEdit('')
    obj.ref_traj.addWidget(obj.lonedit, 1, 3)

    azlabel = QLabel('Az: ')
    obj.ref_traj.addWidget(azlabel, 2, 0)
    obj.azedit = QLineEdit('')
    obj.ref_traj.addWidget(obj.azedit, 2, 1)

    zelabel = QLabel('Ze: ')
    obj.ref_traj.addWidget(zelabel, 2, 2)
    obj.zeedit = QLineEdit('')
    obj.ref_traj.addWidget(obj.zeedit, 2, 3)

    velabel = QLabel('V: ')
    obj.ref_traj.addWidget(velabel, 3, 0)
    obj.veedit = QLineEdit('')
    obj.ref_traj.addWidget(obj.veedit, 3, 1)

    vflabel = QLabel('Vf: ')
    obj.ref_traj.addWidget(vflabel, 4, 0)
    obj.vfedit = QLineEdit('')
    obj.ref_traj.addWidget(obj.vfedit, 4, 1)

    tilabel = QLabel('T: ')
    obj.ref_traj.addWidget(tilabel, 3, 2)
    obj.tiedit = QLineEdit('')
    obj.ref_traj.addWidget(obj.tiedit, 3, 3)

    fragLabel = QLabel("Fragmentations")
    obj.ref_frags.addWidget(fragLabel)

    obj.ref_table = QTableWidget()
    obj.ref_frags.addWidget(obj.ref_table)
    defTable(obj.ref_table, 1, 5, headers=['Height', 'Latitude', 'Longitude', 'Time', 'δ Height'])
    obj.ref_table.setStyleSheet('')

    table_control = QHBoxLayout()
    obj.ref_frags.addLayout(table_control)

    frag_table_sub = QPushButton("-")
    table_control.addWidget(frag_table_sub)
    frag_table_sub.clicked.connect(partial(changeRows, obj.ref_table, -1))
    frag_table_sub.setToolTip("Remove row")

    frag_table_add = QPushButton("+")
    table_control.addWidget(frag_table_add)
    frag_table_add.clicked.connect(partial(changeRows, obj.ref_table, 1))
    frag_table_add.setToolTip("Add row")

    frag_table_save = QPushButton("Save")
    table_control.addWidget(frag_table_save)
    frag_table_save.clicked.connect(obj.refSave)

    frag_table_load = QPushButton("Load")
    table_control.addWidget(frag_table_load)
    frag_table_load.clicked.connect(obj.refLoad)

    load_stations = QPushButton("Load All Stations")
    obj.ref_buttons.addWidget(load_stations)
    load_stations.clicked.connect(obj.refLoadStations)

    load_traj = QPushButton("Load Trajectory")
    obj.ref_buttons.addWidget(load_traj)
    load_traj.clicked.connect(obj.refLoadTrajectory)

    high_stats = QPushButton("Highlight Stations")
    obj.ref_buttons.addWidget(high_stats)
    high_stats.clicked.connect(obj.refHighStats)

    sync_heights = QPushButton("Sync Heights to Trajectory")
    obj.ref_buttons.addWidget(sync_heights)
    sync_heights.clicked.connect(obj.refSyncHeights) 

    plot_heights = QPushButton("Plot Fragmentation Times")
    obj.ref_buttons.addWidget(plot_heights)
    plot_heights.clicked.connect(obj.refPlotHeights) 

    plot_traj = QPushButton("Plot Trajectory Times")
    obj.ref_buttons.addWidget(plot_traj)
    plot_traj.clicked.connect(obj.refPlotTraj)

    clear_stations = QPushButton("Clear Stations")
    obj.ref_wave_control.addWidget(clear_stations)
    clear_stations.clicked.connect(obj.refClearStations)

    ball_stations = QPushButton("Show Ballistic")
    obj.ref_wave_control.addWidget(ball_stations)
    ball_stations.clicked.connect(obj.refBallStations)

    frag_stations = QPushButton("Show Fragmentations")
    obj.ref_wave_control.addWidget(frag_stations)
    frag_stations.clicked.connect(obj.refFragStations)

    both_stations = QPushButton("Show Fragmentation and Ballistic")
    obj.ref_wave_control.addWidget(both_stations)
    both_stations.clicked.connect(obj.refBothStations)


    obj.ref_control.addLayout(obj.ref_traj)
    obj.ref_control.addLayout(obj.ref_frags)
    obj.ref_control.addLayout(obj.ref_buttons)

    obj.master_ref.addLayout(obj.ref_control)
    obj.master_ref.addLayout(obj.ref_wave_layout)
    
    obj.ref_wave_layout.addWidget(obj.ref_waveforms)
    obj.ref_wave_layout.addLayout(obj.ref_wave_control)

    obj.tab_widget.addTab(refine_tab, "Refine")

def addDocWidgets(obj):
    doc_master_tab = QWidget()
    doc_master = QVBoxLayout()
    doc_master_tab.setLayout(doc_master)

    obj.doc_view = QWebView()
    doc_master.addWidget(obj.doc_view)
    obj.doc_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

    obj.tab_widget.addTab(doc_master_tab, 'Documentation')  

def addMakePicksWidgets(obj):
    
    make_picks_master_tab = QWidget()
    make_picks_master = QVBoxLayout()
    make_picks_master_tab.setLayout(make_picks_master)

    obj.make_picks_top_graphs = QHBoxLayout()
    obj.make_picks_bottom_graphs = QGridLayout()
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
    obj.make_picks_bottom_graphs.addWidget(obj.make_picks_waveform_view, 1, 2)
    obj.make_picks_waveform_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)


    toggle_button_array = QVBoxLayout()
    toggle_button_array.setSpacing(0)
    obj.make_picks_bottom_graphs.addLayout(toggle_button_array, 1, 1)

    obj.tog_picks = ToggleButton(False, 1)
    obj.tog_picks.setToolTip("Click the waveform to make a pick")
    obj.tog_picks.clicked.connect(obj.tog_picks.clickedEvt)
    toggle_button_array.addWidget(obj.tog_picks)

    obj.tog_rm_picks = ToggleButton(False, 4)
    obj.tog_rm_picks.setToolTip("Click the waveform to remove a pick")
    obj.tog_rm_picks.clicked.connect(obj.tog_rm_picks.clickedEvt)
    toggle_button_array.addWidget(obj.tog_rm_picks)

    # obj.annote_picks = ToggleButton(False, 2)
    # obj.annote_picks.setToolTip("Click the waveform to make an annotation")
    # obj.annote_picks.clicked.connect(obj.annote_picks.clickedEvt)
    # toggle_button_array.addWidget(obj.annote_picks)

    # obj.gnd_mot_picks = ToggleButton(False, 3)
    # obj.gnd_mot_picks.setToolTip("Click the waveform to get the ground motion")
    # obj.gnd_mot_picks.clicked.connect(obj.gnd_mot_picks.clickedEvt)
    # toggle_button_array.addWidget(obj.gnd_mot_picks)

    # obj.bandpass_picks = ToggleButton(False, 6)
    # obj.bandpass_picks.setToolTip("Click the waveform to open optimal bandpass dialog")
    # obj.bandpass_picks.clicked.connect(obj.bandpass_picks.clickedEvt)
    # toggle_button_array.addWidget(obj.bandpass_picks)

    # obj.polmap_picks = ToggleButton(False, 7)
    # obj.polmap_picks.setToolTip("Click the waveform to open polarization heat map dialog")
    # obj.polmap_picks.clicked.connect(obj.polmap_picks.clickedEvt)
    # toggle_button_array.addWidget(obj.polmap_picks)

    toggle_button_array.insertStretch(-1, 0)

    make_picks_control_panel = QHBoxLayout()
    obj.make_picks_bottom_graphs.addLayout(make_picks_control_panel, 2, 1, 1, 2)

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

    obj.low_bandpass_label, obj.low_bandpass_edits = createLabelEditObj('Low: ', filter_group_layout, 0)
    obj.high_bandpass_label, obj.high_bandpass_edits = createLabelEditObj('High: ', filter_group_layout, 1)

    obj.low_bandpass_edits.setText('2')
    obj.high_bandpass_edits.setText('8')

    obj.low_bandpass_edits.textChanged.connect(obj.updatePlot)
    obj.high_bandpass_edits.textChanged.connect(obj.updatePlot)

    obj.filter_combo_box = QComboBox()
    filter_group_layout.addWidget(obj.filter_combo_box, 2, 0, 1, 4)

    make_picks_picks_group = QGroupBox("Arrival Picks")
    make_picks_control_panel.addWidget(make_picks_picks_group)

    pick_group_layout = QGridLayout()
    make_picks_picks_group.setLayout(pick_group_layout)

    obj.export_to_csv = QPushButton('Export Picks to CSV')
    pick_group_layout.addWidget(obj.export_to_csv)

    obj.export_to_all_times = QPushButton('Export Arrival Times')
    pick_group_layout.addWidget(obj.export_to_all_times)

    obj.export_to_image = QPushButton('Export Image')
    pick_group_layout.addWidget(obj.export_to_image)
    obj.export_to_image.clicked.connect(obj.exportImage)

    obj.W_est = QPushButton('Yield Estimate')
    pick_group_layout.addWidget(obj.W_est)
    obj.W_est.clicked.connect(obj.W_estGUI)

    make_picks_check_group = QGroupBox("Toggles")
    make_picks_control_panel.addWidget(make_picks_check_group)

    check_group_layout = QVBoxLayout()
    make_picks_check_group.setLayout(check_group_layout)

    obj.show_frags = QCheckBox('Show Fragmentations')
    check_group_layout.addWidget(obj.show_frags)
    obj.show_frags.stateChanged.connect(partial(obj.updatePlot, True))

    obj.show_ball = QCheckBox('Show Ballistic')
    check_group_layout.addWidget(obj.show_ball)
    obj.show_ball.stateChanged.connect(partial(obj.updatePlot, True))

    # obj.show_prec = QCheckBox('Show Precursors')
    # check_group_layout.addWidget(obj.show_prec)
    # obj.show_prec.stateChanged.connect(partial(obj.updatePlot, True))

    # obj.rm_resp = QCheckBox('Remove Response (EXPERIMENTAL)')
    # check_group_layout.addWidget(obj.rm_resp)
    # obj.rm_resp.stateChanged.connect(partial(obj.updatePlot, True))

    obj.show_perts = QCheckBox('Show Perturbations')
    check_group_layout.addWidget(obj.show_perts)
    obj.show_perts.stateChanged.connect(partial(obj.updatePlot, True))

    obj.show_height = QCheckBox('Show Height Prediction')
    check_group_layout.addWidget(obj.show_height)

    # obj.psd = QCheckBox('[EXPERIMENTAL] PSD')
    # check_group_layout.addWidget(obj.psd)
    # obj.psd.stateChanged.connect(obj.psdPlot)

    # obj.solve_height = QCheckBox('Solve Heights')
    # check_group_layout.addWidget(obj.solve_height)

    make_picks_plot_tools_group = QGroupBox("Plot Tools")
    make_picks_control_panel.addWidget(make_picks_plot_tools_group)

    plot_tweaks_layout = QVBoxLayout()
    make_picks_plot_tools_group.setLayout(plot_tweaks_layout)

    obj.invert = QCheckBox('Invert') 
    plot_tweaks_layout.addWidget(obj.invert)
    obj.invert.stateChanged.connect(obj.invertGraph)

    # obj.show_title = QCheckBox('Show Title')
    # plot_tweaks_layout.addWidget(obj.show_title)
    # obj.show_title.stateChanged.connect(obj.showTitle)

    obj.show_contour = QPushButton('Show Ballistic Contour')
    plot_tweaks_layout.addWidget(obj.show_contour)
    obj.show_contour.clicked.connect(partial(obj.showContour, 'ballistic'))

    obj.show_f_contour = QPushButton('Show Fragmentation Contour')
    plot_tweaks_layout.addWidget(obj.show_f_contour)
    obj.show_f_contour.clicked.connect(partial(obj.showContour, 'fragmentation'))

    # obj.save_contour = QPushButton('Save Contour')
    # plot_tweaks_layout.addWidget(obj.save_contour)
    # obj.save_contour.clicked.connect(obj.saveContour)

    # obj.load_contour = QPushButton('Load Contour')
    # plot_tweaks_layout.addWidget(obj.load_contour)
    # obj.load_contour.clicked.connect(obj.loadContour)

    # obj.clear_contour = QPushButton('Clear Contour')
    # plot_tweaks_layout.addWidget(obj.clear_contour)
    # obj.clear_contour.clicked.connect(obj.clearContour)

    obj.tab_widget.addTab(make_picks_master_tab, 'Make Picks')  

def addSeisTrajWidgets(obj):

    seis_tab = QWidget()
    obj.master_seis = QHBoxLayout()
    obj.seis_tab_input = QVBoxLayout()
    obj.seis_tab_output = QGridLayout()

    seis_tab.setLayout(obj.master_seis)
    obj.tab_widget.addTab(seis_tab, "Seismic Trajectory")

    obj.master_seis.addLayout(obj.seis_tab_input)
    obj.master_seis.addLayout(obj.seis_tab_output)

    table_group = QGridLayout()
    obj.seis_tab_input.addLayout(table_group)

    tab_layout = QGridLayout()
    obj.seis_tab_input.addLayout(tab_layout)

    obj.seis_search = QPushButton('Search')
    tab_layout.addWidget(obj.seis_search, 0, 1, 1, 100)
    obj.seis_search.clicked.connect(obj.trajSearchSetup)
    #obj.seis_search.clicked.connect(obj.seisSearch)
    
    obj.seis_table = QTableWidget()
    table_group.addWidget(obj.seis_table, 1, 1, 1, 100)

    obj.seis_resids = QTableWidget()
    table_group.addWidget(obj.seis_resids, 2, 1, 1, 100)

    # obj.seis_plane = QTableWidget()
    # table_group.addWidget(obj.seis_plane, 3, 1, 1, 100)
    # defTable(obj.seis_plane, 3, 3, headers=['Latitude', 'Longitude', 'Height'])
    #obj.seis_three_canvas = FigureCanvas(Figure(figsize=(5, 5)))
    # obj.seis_three_canvas.setSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)

    obj.seis_view = gl.GLViewWidget()
    #view.setMinimumSize(384,360)
    xgrid = gl.GLGridItem()
    xgrid.setSize(x=100000, y=100000, z=100000)
    xgrid.setSpacing(x=10000, y=10000, z=10000)
    obj.seis_view.addItem(xgrid)
    obj.seis_tab_output.addWidget(obj.seis_view, 1, 1, 50, 1)

    two_graphs = QGridLayout()
    obj.seis_tab_output.addLayout(two_graphs, 51, 1, 50, 1)    

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

    obj.fatm_variable_combo.addItem('Sound Speed')
    obj.fatm_variable_combo.addItem('Effective Sound Speed')
    obj.fatm_variable_combo.addItem('U-Component of Wind')
    obj.fatm_variable_combo.addItem('V-Component of Wind')
    obj.fatm_variable_combo.addItem('Wind Magnitude')
    obj.fatm_variable_combo.addItem('Wind Direction')

    
    try:
        obj.fatm_variable_combo.currentTextChanged.connect(obj.fatmLoadAtm)
    except:
        pass

    obj.fatm_name_label, obj.fatm_name_edits = createLabelEditObj('Name:', fetch_content, 1)

    obj.fatm_button = QPushButton('Browse')
    fetch_content.addWidget(obj.fatm_button, 1, 3)
    obj.fatm_button.clicked.connect(partial(fileSearch, ['NetCDF (*.nc)', 'HDF (*.HDF)', 'CSV (*.csv)', 'TXT (*.txt)'], obj.fatm_name_edits))

    obj.fatm_datetime_label = QLabel("Time of Profile:")
    obj.fatm_datetime_edits = QDateTimeEdit()
    fetch_content.addWidget(obj.fatm_datetime_label, 2, 1)
    fetch_content.addWidget(obj.fatm_datetime_edits, 2, 2)
    obj.fatm_datetime_edits.setCalendarPopup(True)

    obj.fatm_source_type = QComboBox()
    obj.fatm_source_type.addItem("Copernicus Climate Change Service (ECMWF)")
    obj.fatm_source_type.addItem("Copernicus Climate Change Service (ECMWF) - Spread File")
    obj.fatm_source_type.addItem("Radiosonde")
    fetch_content.addWidget(obj.fatm_source_type, 3, 1, 1, 2)

    save_into_bam = QPushButton('Save into BAM')
    fetch_content.addWidget(save_into_bam, 4, 1, 1, 1)
    save_into_bam.clicked.connect(obj.fatmSaveAtm)

    load_from_bam = QPushButton('Load from BAM')
    fetch_content.addWidget(load_from_bam, 4, 2, 1, 1)
    load_from_bam.clicked.connect(obj.fatmLoadAtm)

    #############################

    # indep_group = QGroupBox('Inependent Variables (x)')
    # fetch_content.addWidget(indep_group, 5, 1, 1, 2)

    # group_box = QGridLayout()
    # indep_group.setLayout(group_box)

    # obj.fatm_temp = QCheckBox('Temperature')
    # group_box.addWidget(obj.fatm_temp, 1, 1)

    # obj.fatm_u_wind = QCheckBox('U Wind')
    # group_box.addWidget(obj.fatm_u_wind, 1, 2)

    # obj.fatm_v_wind = QCheckBox('V Wind')
    # group_box.addWidget(obj.fatm_v_wind, 1, 3)

    # ##############################

    # dep_group = QGroupBox('Dependent Variables (y)')
    # fetch_content.addWidget(dep_group, 6, 1, 1, 2)

    # dgroup_box = QGridLayout()
    # dep_group.setLayout(dgroup_box)

    # obj.fatm_geo_height = QCheckBox('Geopotential Height')
    # dgroup_box.addWidget(obj.fatm_geo_height, 1, 1)

    # obj.fatm_pressure = QCheckBox('Pressure')
    # dgroup_box.addWidget(obj.fatm_pressure, 1, 2)

    ###############################
    # op_group = QGroupBox('Options')
    # fetch_content.addWidget(op_group, 7, 1, 1, 2)

    # opgroup_box = QGridLayout()
    # op_group.setLayout(opgroup_box)

    # obj.fatm_perts = QCheckBox('Show Perturbations')
    # opgroup_box.addWidget(obj.fatm_perts, 2, 1)
    #############################################

    obj.fatm_fetch = QPushButton("Download")
    fetch_content.addWidget(obj.fatm_fetch, 8, 1, 1, 2)
    obj.fatm_fetch.clicked.connect(partial(obj.fatmFetch, True))

    # obj.fatm_open = QPushButton("Open")
    # fetch_content.addWidget(obj.fatm_open, 9, 1, 1, 2)
    # obj.fatm_open.clicked.connect(partial(obj.fatmFetch, False))

    obj.fatm_print = QPushButton("Print")
    fetch_content.addWidget(obj.fatm_print, 16, 1, 1, 2)
    obj.fatm_print.clicked.connect(obj.fatmPrint)

    obj.fatm_start_label = QLabel("Start lat/lon/elev")
    obj.fatm_start_lat = QLineEdit()
    obj.fatm_start_lon = QLineEdit()
    obj.fatm_start_elev = QLineEdit()
    fetch_content.addWidget(obj.fatm_start_label, 10, 1)
    fetch_content.addWidget(obj.fatm_start_lat, 10, 2)
    fetch_content.addWidget(obj.fatm_start_lon, 11, 2)
    fetch_content.addWidget(obj.fatm_start_elev, 12, 2)

    obj.fatm_end_label = QLabel("End lat/lon/elev")
    obj.fatm_end_lat = QLineEdit()
    obj.fatm_end_lon = QLineEdit()
    obj.fatm_end_elev = QLineEdit()
    fetch_content.addWidget(obj.fatm_end_label, 13, 1)
    fetch_content.addWidget(obj.fatm_end_lat, 13, 2)
    fetch_content.addWidget(obj.fatm_end_lon, 14, 2)
    fetch_content.addWidget(obj.fatm_end_elev,15, 2)

    # self.LineEdit.setValidator(QIntValidator())

    # obj.fatm_lat_label = QLabel("Latitude: 0")
    # obj.fatm_lat_slide = QSlider(Qt.Horizontal)
    # fetch_content.addWidget(obj.fatm_lat_label, 7, 1)
    # fetch_content.addWidget(obj.fatm_lat_slide, 7, 2)
    # obj.fatm_lat_slide.setMinimum(-90/obj.slider_scale)
    # obj.fatm_lat_slide.setMaximum(90/obj.slider_scale)
    # obj.fatm_lat_slide.setValue(0)
    # obj.fatm_lat_slide.setTickInterval(0.5)
    # obj.fatm_lat_slide.valueChanged.connect(partial(obj.fatmValueChange, obj.fatm_lat_label, obj.fatm_lat_slide))

    # obj.fatm_lon_label = QLabel("Longitude: 0")
    # obj.fatm_lon_slide = QSlider(Qt.Horizontal)
    # fetch_content.addWidget(obj.fatm_lon_label, 8, 1)
    # fetch_content.addWidget(obj.fatm_lon_slide, 8, 2)
    # obj.fatm_lon_slide.setMinimum(-180/obj.slider_scale)
    # obj.fatm_lon_slide.setMaximum(180/obj.slider_scale)
    # obj.fatm_lon_slide.setValue(0)
    # obj.fatm_lon_slide.setTickInterval(0.5)
    # obj.fatm_lon_slide.valueChanged.connect(partial(obj.fatmValueChange, obj.fatm_lon_label, obj.fatm_lon_slide))

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
    obj.atm_view.setBackground((0, 0, 0))
    obj.atm_canvas.getAxis('bottom').setPen((255, 255, 255)) 
    obj.atm_canvas.getAxis('left').setPen((255, 255, 255))

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

    obj.ray_line_canvas = FigureCanvas(Figure())
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

    obj.ray_enable_vars = QCheckBox('Enable Ray Grid')
    obj.ray_control.addWidget(obj.ray_enable_vars, 3, 5)
    obj.ray_enable_vars.stateChanged.connect(obj.trajSolver)

    obj.ray_canvas.scene().sigMouseClicked.connect(obj.rayMouseClicked)

