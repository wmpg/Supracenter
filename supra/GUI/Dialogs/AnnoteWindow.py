
import numpy as np

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.GUITools import *
from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.WidgetBuilder import center
from supra.Utils.Classes import Annote
from supra.Stations.ProcessStation import *
from supra.Utils.Formatting import *
from supra.Files.SaveLoad import save
class AnnoteWindow(QWidget):

    def __init__(self, time, stn, bam, mode="new", an=None, current_channel=None):

        QWidget.__init__(self)
        
        self.bam = bam

        self.buildGUI(time)

        self.annote_color = QColor(0, 0, 255)

        self.stn = stn

        self.an = an

        self.mode = mode
        
        if mode == "edit":
            trace = self.stn.stream[0]

            start_time = self.an.time
            end_time = self.an.time + self.an.length

            if not hasattr(self.an, 'bandpass'):
                self.an.bandpass = None
            self.bandpass = self.an.bandpass

            clean = {"ref_datetime" : self.bam.setup.fireball_datetime, "resp" : self.stn.response, "bandpass" : None, "backup" : False}

            raw_waveform, raw_time = subTrace(trace, start_time, end_time, self.bam.setup.fireball_datetime, clean=clean)

            clean = {"ref_datetime" : self.bam.setup.fireball_datetime, "resp" : self.stn.response, "bandpass" : self.bandpass, "backup" : False}

            cut_waveform, cut_time = subTrace(trace, start_time, end_time, self.bam.setup.fireball_datetime, clean=clean)
            station_waveform = pg.PlotDataItem(x=cut_time, y=cut_waveform, pen='w')
            self.annote_waveform_canvas.addItem(station_waveform)

            start_time_noise = 0
            end_time_noise = 0 + self.an.length

            clean = {"ref_datetime" : self.bam.setup.fireball_datetime, "resp" : self.stn.response, "bandpass" : self.bandpass, "backup" : False}

            noise_waveform, cut_time = subTrace(trace, start_time_noise, end_time_noise, self.bam.setup.fireball_datetime, clean=clean)

            # cut_waveform =  [item for sublist in cut_waveform for item in sublist]
            # cut_time =      [item for sublist in cut_time for item in sublist]
            # noise_waveform =[item for sublist in noise_waveform for item in sublist]

            try:
                sampling_rate = stn.stream.select(channel=current_channel)[0].stats.sampling_rate
                freq, fas = genFFT(raw_waveform, sampling_rate)
                station_fft = pg.PlotDataItem(x=freq, y=fas, pen="w")
                freq_n, fas_n = genFFT(noise_waveform, sampling_rate)
                station_noise_fft = pg.PlotDataItem(x=freq_n, y=fas_n, pen="r")
                station_ratio_fft = pg.PlotDataItem(x=freq_n, y=fas/fas_n, pen="b")

                self.annote_fft_canvas.addItem(pg.InfiniteLine(pos=(0, 0), angle=0, pen=QColor(255, 255, 255)))
                self.annote_fft_canvas.addItem(station_fft)
                self.annote_fft_canvas.addItem(station_noise_fft)
                self.annote_fft_canvas.addItem(station_ratio_fft)
            except (IndexError, ValueError):
                print(printMessage("error"), " Not enough time data for FFT")



            self.loadAnnote(an)

    def color_picker(self):
        color = QColorDialog.getColor()
        
        self.color_button.setStyleSheet("background-color: %s" % color.name())
        self.annote_color = color

    def loadAnnote(self, an):
        self.title_edits.setText(an.title)
        self.time_edits.setText(str(an.time))
        self.length_edits.setText(str(an.length))
        comboSet(self.group_edits, str(an.group))
        comboSet(self.source_edits, str(an.source))
        self.height_edits.setText(str(an.height))
        self.notes_box.setPlainText(an.notes)

        if an.bandpass is not None:
            self.low_bandpass_edits.setText(str(an.bandpass[0]))
            self.high_bandpass_edits.setText(str(an.bandpass[1]))
        else:
            self.low_bandpass_edits.setText(str(2))
            self.high_bandpass_edits.setText(str(8))

    def addAnnote(self):
        
        title = self.title_edits.text()
        time = float(self.time_edits.text())
        length = float(self.length_edits.text())
        group = self.group_edits.currentText()
        source = self.source_edits.currentText()
        height = float(self.height_edits.text())
        notes = self.notes_box.toPlainText()
        color = self.annote_color
        bandpass = [float(self.low_bandpass_edits.text()), float(self.high_bandpass_edits.text())]

        an = Annote(title, time, length, group, source, height, notes, color, bandpass)

        if self.mode == 'new':
            self.stn.annotation.add(an)
        elif self.mode == 'edit':
            self.stn.annotation.overwrite(an)
        save(self.bam, file_check=False)

        self.close()

    def delAnnote(self):
        self.stn.annotation.remove(self.an)
        self.close()

    # def groupRefresh(self):
    #     for s in self.bam.source_list:
    #         if self.group_edits.currentText() == s.title:

    #             comboSet(self.source_edits, s.source_type)

    #             if s.source_type == "Fragmentation":
    #                 self.height_edits.setText(str(s.source.position.elev))

    #     AnnoteWindow.update(self)

    def buildGUI(self, time):
        self.setWindowTitle('Edit Annotation')
        
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)
     
        title_label = QLabel("Title: ")
        layout.addWidget(title_label, 1, 1, 1, 1)

        self.title_edits = QLineEdit()
        layout.addWidget(self.title_edits, 1, 2, 1, 1)

        self.color_button = QPushButton()
        layout.addWidget(self.color_button, 1, 3, 1, 1)
        self.color_button.clicked.connect(self.color_picker)

        group_label = QLabel("Group: ")
        layout.addWidget(group_label, 2, 1, 1, 1)

        self.group_edits = QComboBox()
        layout.addWidget(self.group_edits, 2, 2, 1, 1)

        _, self.low_bandpass_edits =  createLabelEditObj("Lowpass", layout, 3, width=1, h_shift=0, tool_tip='', validate='float', default_txt='2')
        _, self.high_bandpass_edits = createLabelEditObj("Highpass", layout, 4, width=1, h_shift=0, tool_tip='', validate='float', default_txt='8')

        self.group_edits.addItem("[Unknown]")

        if hasattr(self.bam, "source_list"):
            for s in self.bam.source_list:
                self.group_edits.addItem(s.title)

        time_label = QLabel("Time: ")
        layout.addWidget(time_label, 5, 1, 1, 1)

        self.time_edits = QLineEdit("{:.2f}".format(time))
        layout.addWidget(self.time_edits, 5, 2, 1, 1)
        self.time_edits.setValidator(QDoubleValidator())

        length_label = QLabel("Length: ")
        layout.addWidget(length_label, 6, 1, 1, 1)

        self.length_edits = QLineEdit('0')
        layout.addWidget(self.length_edits, 6, 2, 1, 1)
        self.length_edits.setValidator(QDoubleValidator())

        source_label = QLabel("Source: ")
        layout.addWidget(source_label, 7, 1, 1, 1)

        self.source_edits = QComboBox()
        layout.addWidget(self.source_edits, 7, 2, 1, 1)
        for s in ["[None]", "Fragmentation", "Ballisitc", "Indirect"]:
            self.source_edits.addItem(s)

        height_label = QLabel("Height: ")
        layout.addWidget(height_label, 8, 1, 1, 1)

        self.height_edits = QLineEdit('0')
        layout.addWidget(self.height_edits, 8, 2, 1, 1)
        self.height_edits.setValidator(QDoubleValidator())

        notes_label = QLabel("Notes: ")
        layout.addWidget(notes_label, 9, 1, 1, 1)

        self.notes_box = QPlainTextEdit()
        layout.addWidget(self.notes_box, 9, 2, 1, 2)

        # self.group_edits.currentIndexChanged.connect(self.groupRefresh)

        self.annote_waveform_view = pg.GraphicsLayoutWidget()
        self.annote_waveform_canvas = self.annote_waveform_view.addPlot()
        layout.addWidget(self.annote_waveform_view, 10, 1, 1, 4)
        self.annote_waveform_canvas.setLabel('bottom', "Time after Reference", units='s')
        self.annote_waveform_canvas.setLabel('left', "Amplitude")

        self.annote_fft_view = pg.GraphicsLayoutWidget()
        self.annote_fft_canvas = self.annote_fft_view.addPlot()
        layout.addWidget(self.annote_fft_view, 2, 3, 9, 4)
        self.annote_fft_canvas.setLogMode(True, True)
        self.annote_fft_canvas.setLabel('bottom', "Frequency", units='Hz')
        self.annote_fft_canvas.setLabel('left', "Amplitude")

        save_button = QPushButton('Save')
        layout.addWidget(save_button, 12, 2, 1, 1)
        save_button.clicked.connect(self.addAnnote)

        del_button = QPushButton('Delete')
        layout.addWidget(del_button, 12, 1, 1, 1)
        del_button.clicked.connect(self.delAnnote)

if __name__ == '__main__':

    pass