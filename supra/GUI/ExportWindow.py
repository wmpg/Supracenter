import numpy as np

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.WidgetBuilder import theme

class ExportWindow(QScrollArea):
    def __init__(self, invert, setup, stn_list, position):
        self.link = False
        self.grid = False
        QWidget.__init__(self)
        self.selected_stn_view = pg.GraphicsLayoutWidget()
        self.stn_canvas = []
        self.stn_list = stn_list
        self.setup = setup
        self.waveform_data = [None]*len(stn_list)
        ref_pos = Position(self.setup.lat_centre, self.setup.lon_centre, 0)
        count = 0
        max_steps = len(position)*len(stn_list)*setup.perturb_times
        dataset = parseWeather(self.setup)

        theme(self)
        for ii, index in enumerate(invert):

            if index:
                stn = stn_list[ii]
                self.stn_canvas.append(self.selected_stn_view.addPlot())
                self.selected_stn_view.nextRow()
                min_point, max_point = self.discountDrawWaveform(setup, ii, self.stn_canvas[-1])
                self.stn_canvas[-1].getAxis('bottom').setPen((0, 0, 0)) 
                self.stn_canvas[-1].getAxis('left').setPen((0, 0, 0)) 
                self.waveform_data[ii].setPen((0, 0, 0))
                #self.stn_canvas[-1].addItem(pg.LabelItem(text="{:}-{:}".format(stn.network, stn.code), color=(0, 0, 0)))
                self.stn_canvas[-1].setTitle("{:}-{:}".format(stn.network, stn.code), color=(0, 0, 0))

                for p, point in enumerate(position):
                    for ptb_n in range(self.setup.perturb_times):             
                        dataset = SolutionGUI.perturbGenerate(self, ptb_n, dataset, SolutionGUI.perturbSetup(self))
                        zProfile, _ = getWeather(np.array([point.lat, point.lon, point.elev]), np.array([stn.position.lat, stn.position.lon, stn.position.elev]), self.setup.weather_type, \
                                            [ref_pos.lat, ref_pos.lon, ref_pos.elev], dataset, convert=False)
                        point.pos_loc(ref_pos)
                        stn.position.pos_loc(ref_pos)
                        f_time, _, _, _ = cyscan(np.array([point.x, point.y, point.z]), np.array([stn.position.x, stn.position.y, stn.position.z]), zProfile, wind=True, \
                            n_theta=self.setup.n_theta, n_phi=self.setup.n_phi, h_tol=self.setup.h_tol, v_tol=self.setup.v_tol)

                        A = self.setup.trajectory.pos_i
                        B = self.setup.trajectory.pos_f

                        A.pos_loc(B)
                        B.pos_loc(B)

                        # Get prediction of time of the meteor, so the timing of each fragmentation can be known
                        length_of_meteor = np.sqrt((A.x - B.x)**2 + (A.y - B.y)**2 + (A.z - B.z)**2)
                        time_of_meteor = length_of_meteor/self.setup.trajectory.v

                        correction = time_of_meteor - A.z/self.setup.trajectory.pos_i.elev*(time_of_meteor)

                        nom_time = f_time + correction
                        with warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            try:
                                self.stn_canvas[-1].plot(x=[f_time, f_time], y=[min_point, max_point], pen=PEN[p%7], brush=PEN[p%7])
                            except:
                                pass
                        count += 1
                        loadingBar("Generating Plots", count, max_steps)
                        try:
                            self.stn_canvas[-1].setXRange(nom_time-25, nom_time+25, padding=1)
                        except:
                            pass

        self.selected_stn_view.setBackground((255, 255, 255))


        layout = QVBoxLayout()
        layout.addWidget(self.selected_stn_view)
        self.setLayout(layout)

        toggle_group = QGroupBox("Toggles")
        layout.addWidget(toggle_group)

        toggle_group_layout = QVBoxLayout()
        toggle_group.setLayout(toggle_group_layout)

        sync_button = QCheckBox('Sync')
        toggle_group_layout.addWidget(sync_button)
        sync_button.stateChanged.connect(self.sync)

        grid_button = QCheckBox('Grid')
        toggle_group_layout.addWidget(grid_button)
        grid_button.stateChanged.connect(self.grid_toggle)

        export_button = QPushButton("Export")
        export_button.clicked.connect(self.export)
        layout.addWidget(export_button)


    def export(self):

        # set export parameters if needed
        #exporter.parameters()['width'] = 1000   # (note this also affects height parameter)
        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = dlg[0]

        exporter = pg.exporters.SVGExporter(self.selected_stn_view.scene())

        file_name = file_name + '.svg'
        exporter.export(file_name)


    def sync(self):
        if self.link:
            for i in range(len(self.stn_canvas)):
                self.stn_canvas[i].setXLink(None)
                self.stn_canvas[i].showAxis('bottom') 
        else:
            for i in range(len(self.stn_canvas)):
                self.stn_canvas[i].setXLink(self.stn_canvas[0])
                self.stn_canvas[i].hideAxis('bottom')    
                self.stn_canvas[-1].showAxis('bottom')

        self.link = not self.link

    def grid_toggle(self):
        if self.grid:
            for i in range(len(self.stn_canvas)):
                self.stn_canvas[i].showGrid(x=False, y=False)
        else:
            for i in range(len(self.stn_canvas)):
                self.stn_canvas[i].showGrid(x=True, y=False)

        self.grid = not self.grid

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