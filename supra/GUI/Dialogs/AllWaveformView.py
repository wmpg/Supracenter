import numpy as np
import warnings

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

import obspy
import os
from functools import partial

import scipy.signal
import pyqtgraph.exporters


import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})


from supra.Supracenter.cyscan5 import cyscan

from supra.GUI.Tools.GUITools import *
from supra.GUI.Tools.WidgetBuilder import *

from supra.Utils.Formatting import *
from supra.Utils.TryObj import *

from supra.GUI.Dialogs.ExportWindow import ExportWindow
from supra.GUI.Tools.Theme import theme

from supra.Atmosphere.Parse import parseWeather

HEIGHT_SOLVER_DIV = 100
PEN = [(0     *255, 0.4470*255, 0.7410*255),        
       (0.8500*255, 0.3250*255, 0.0980*255),            
       (0.9290*255, 0.6940*255, 0.1250*255),          
       (0.4940*255, 0.1840*255, 0.5560*255),                
       (0.4660*255, 0.6740*255, 0.1880*255),                
       (0.3010*255, 0.7450*255, 0.9330*255),                
       (0.6350*255, 0.0780*255, 0.1840*255)]

class AllWaveformViewer(QScrollArea):

    def __init__(self, setup, stn_list, position, pert_idx):

        QWidget.__init__(self)
        self.setWindowTitle('Height Solver')
        self.setup = setup
        self.stn_list = stn_list
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        menu_bar = QMenuBar(self) 
        #layout.addWidget(menu_bar, 0, 1)
        file_menu = menu_bar.addMenu('&File')
        select_menu = menu_bar.addMenu('&Select')

        self.selected = False

        file_export = QAction("Export Selected Menu", self)
        file_export.setShortcut('Ctrl+E')
        file_export.setStatusTip('Brings up the Export Menu with the selected stations')
        file_export.triggered.connect(self.export)
        file_menu.addAction(file_export)

        select_all = QAction("Select All", self)
        select_all.setShortcut('Ctrl+A')
        select_all.setStatusTip('Selects all')
        select_all.triggered.connect(self.select)
        select_menu.addAction(select_all)

        theme(self)

        self.position = position
        point = position
        widget = QWidget()
        self.layout = QVBoxLayout(widget)
        self.layout.setAlignment(Qt.AlignTop)
        for point in position:
            print('Best Position: {:}'.format(point))

        ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)
        count = 0
        max_steps = len(stn_list)*len(position)*self.setup.perturb_times
        dataset = parseWeather(self.setup)
        self.save_button = []
        self.stn_view = []
        self.stn_canvas = []
        self.invert = []
        self.waveform_data = [None]*len(stn_list)


        A = self.setup.trajectory.pos_i
        B = self.setup.trajectory.pos_f

        A.pos_loc(B)
        B.pos_loc(B)

        # Get prediction of time of the meteor, so the timing of each fragmentation can be known
        length_of_meteor = np.sqrt((A.x - B.x)**2 + (A.y - B.y)**2 + (A.z - B.z)**2)
        time_of_meteor = length_of_meteor/self.setup.trajectory.v
        for index in range(len(stn_list)):
        
            stn = stn_list[index]
            station_layout = QGridLayout()
            self.layout.addLayout(station_layout)
            label_layout = QVBoxLayout()
            control_layout = QGridLayout()
            waveform_layout = QVBoxLayout()
            station_layout.addLayout(label_layout, 0, 100, 1, 100)
            station_layout.addLayout(control_layout, 0, 0, 1, 100)
            station_layout.addLayout(waveform_layout, 0, 200, 1, 100)

            label_layout.addWidget(QLabel('Station: {:} - {:}'.format(stn_list[index].network, stn_list[index].code)))
            label_layout.addWidget(QLabel('Channel: {:}'.format(stn_list[index].channel)))
            label_layout.addWidget(QLabel('Name: {:}'.format(stn_list[index].name)))
            label_layout.addWidget(QLabel('Lat: {:}'.format(stn_list[index].position.lat)))
            label_layout.addWidget(QLabel('Lon: {:}'.format(stn_list[index].position.lon)))
            label_layout.addWidget(QLabel('Elev: {:}'.format(stn_list[index].position.elev)))

            self.save_button.append(QPushButton('Select'))
            self.save_button[index].clicked.connect(partial(self.saveButton, index))
            control_layout.addWidget(self.save_button[index], 0, 0)
            control_layout.addWidget(QPushButton('-'), 0, 1)
            control_layout.addWidget(QPushButton('-'), 1, 0)
            control_layout.addWidget(QPushButton('-'), 1, 1)

            self.stn_view.append(pg.GraphicsLayoutWidget())
            self.stn_canvas.append(self.stn_view[index].addPlot())
            waveform_layout.addWidget(self.stn_view[index])
            self.invert.append(False)

            min_point, max_point = self.discountDrawWaveform(setup, index, self.stn_canvas[index])
            # for ptb_n in range(self.setup.perturb_times):
            for p, point in enumerate(position): 
                for ptb_n in range(self.setup.perturb_times):               
                    dataset = self.perturbGenerate(ptb_n, dataset, self.perturbSetup())
                    zProfile, _ = getWeather(np.array([point.lat, point.lon, point.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), self.setup.weather_type, \
                                    ref_pos, dataset, convert=True)
                    point.pos_loc(ref_pos)
                    stn.position.pos_loc(ref_pos)
                    f_time, _, _, _ = cyscan(np.array([point.x, point.y, point.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True, \
                            n_theta=self.setup.n_theta, n_phi=self.setup.n_phi, h_tol=self.setup.h_tol, v_tol=self.setup.v_tol)


                    correction = time_of_meteor - A.z/self.setup.trajectory.pos_i.elev*(time_of_meteor)

                    nom_time = f_time + correction
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        try:
                            self.stn_canvas[index].plot(x=[f_time, f_time], y=[min_point, max_point], pen=PEN[p%7], brush=PEN[p%7])
                        except:
                            pass
                        count += 1
                        loadingBar("Generating Plots", count, max_steps)
                    try:
                        self.stn_canvas[index].setXRange(nom_time-25, nom_time+25, padding=1)
                    except:
                        avg_time = stn.position.pos_distance(point)/330
                        self.stn_canvas[index].setXRange(avg_time-25, avg_time+25, padding=1)

        self.setWidget(widget)
        self.setWidgetResizable(True)

    def yesSelect(self, index):
        self.stn_view[index].setBackground((255, 255, 255))
        self.waveform_data[index].setPen((0, 0, 0))
        self.stn_canvas[index].getAxis('bottom').setPen((0, 0, 0)) 
        self.stn_canvas[index].getAxis('left').setPen((0, 0, 0)) 
        self.invert[index] = True

    def noSelect(self, index):
        self.stn_view[index].setBackground((0, 0, 0))
        self.waveform_data[index].setPen((255, 255, 255))
        self.stn_canvas[index].getAxis('bottom').setPen((255, 255, 255)) 
        self.stn_canvas[index].getAxis('left').setPen((255, 255, 255)) 
        self.invert[index] = False


    def select(self):

        for i in range(len(self.waveform_data)):

            if self.selected:
                self.noSelect(i)

            else:
                self.yesSelect(i)


        self.selected = not self.selected


    def saveButton(self, index):

        if self.invert[index]:
            #Make normal
            self.noSelect(index)
        else:
            #Make inverted
            self.yesSelect(index)


    def discountDrawWaveform(self, setup, station_no, canvas):
        # Extract current station
        stn = self.stn_list[station_no]

        # Get the miniSEED file path
        mseed_file_path = os.path.join(setup.working_directory, setup.fireball_name, stn.file_name)

        # Try reading the mseed file, if it doesn't work, skip to the next frame
        try:
            mseed = obspy.read(mseed_file_path)

        except TypeError:
            if setup.debug:
                print('mseed file could not be read:', mseed_file_path)
            return None

        # Unpact miniSEED data
        delta = mseed[0].stats.delta
        start_datetime = mseed[0].stats.starttime.datetime
        end_datetime = mseed[0].stats.endtime.datetime

        stn.offset = (start_datetime - setup.fireball_datetime).total_seconds()

        waveform_data = mseed[0].data

        # Store raw data for bookkeeping on first open
        self.current_waveform_raw = waveform_data

        self.current_waveform_delta = delta
        self.current_waveform_time = np.arange(0, (end_datetime - start_datetime).total_seconds() + delta, \
            delta)

        # Construct time array, 0 is at start_datetime
        time_data = np.copy(self.current_waveform_time)

        # Cut the waveform data length to match the time data
        waveform_data = waveform_data[:len(time_data)]
        time_data = time_data[:len(waveform_data)] + stn.offset

        # Get bandpass filter values
        bandpass_low = float(2)
        bandpass_high = float(8)


        # Init the butterworth bandpass filter
        butter_b, butter_a = butterworthBandpassFilter(bandpass_low, bandpass_high, \
            1.0/self.current_waveform_delta, order=6)

        # Filter the data
        waveform_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(self.current_waveform_raw))

        # Plot the waveform
        self.waveform_data[station_no] = pg.PlotDataItem(x=time_data, y=waveform_data, pen='w')
        canvas.addItem(self.waveform_data[station_no])
        #canvas.setXRange(t_arrival-100, t_arrival+100, padding=1)
        #canvas.setLabel('bottom', "Time after {:}".format(setup.fireball_datetime), units='s')
        canvas.setLabel('left', "Signal Response")

        return np.min(waveform_data), np.max(waveform_data)

    def export(self):
        self.w = ExportWindow(self.invert, self.setup, self.stn_list, self.position)
        self.w.setGeometry(QRect(100, 100, 900, 900))
        self.w.show()

    def perturbSetup(self):
        """ Pulls the correct file names for the perturbing function to read, depending on the perturbation type
        """

        if self.setup.perturb_method == 'temporal':

            # sounding data one hour later
            sounding_u = parseWeather(self.setup, time= 1)

            # sounding data one hour earlier
            sounding_l = parseWeather(self.setup, time=-1)

        else:
            sounding_u = []
            sounding_l = []

        if self.setup.perturb_method == 'ensemble':
            ensemble_file = self.setup.perturbation_spread_file
        else:
            ensemble_file = ''

        if self.setup.perturb_times == 0: self.setup.perturb_times = 1

        if not self.setup.perturb:
            self.setup.perturb_times = 1

        return np.array([sounding_l, sounding_u, ensemble_file])

    def perturbGenerate(self, ptb_n, dataset, perturb_data, line=False):
        """ Generates a perturbed cubic atmospheric profile (similar to 'dataset') based off of the perturb data 
        """

        sounding_l, sounding_u, ensemble_file = perturb_data[0], perturb_data[1], perturb_data[2]

        # Perturbed soundings
        if ptb_n > 0:
            

            sounding_p = perturbation_method(self.setup, dataset, self.setup.perturb_method, \
                sounding_u=sounding_u, sounding_l=sounding_l, \
                spread_file=self.setup.perturbation_spread_file, lat=self.setup.lat_centre, lon=self.setup.lon_centre, \
                ensemble_file=ensemble_file, ensemble_no=ptb_n, line=line)

        # Nominal sounding
        else:
            sounding_p = dataset


        return sounding_p
if __name__ == '__main__':

    pass