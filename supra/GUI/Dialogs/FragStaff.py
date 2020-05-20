import numpy as np

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.Utils.AngleConv import chauvenet
from supra.GUI.Tools.Theme import theme

class FragmentationStaff(QWidget):

    def __init__(self, setup, pack):

        QWidget.__init__(self)
        self.setWindowTitle('Fragmentation Staff')
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)
        
        theme(self)

        # Take important values from main window class
        stn, self.current_station, self.pick_list = pack

        nom_pick = self.pick_list[0]

        layout = QVBoxLayout()

        # Pass setup value
        self.setup = setup

        ###########
        # Build GUI
        ###########
        self.plot_view = pg.GraphicsLayoutWidget()
        self.height_canvas = self.plot_view.addPlot()
        self.plot_view.nextRow()
        self.angle_canvas = self.plot_view.addPlot()

        self.plot_view.setBackground((0,0,0))
        self.angle_canvas.getAxis('bottom').setPen((255, 255, 255)) 
        self.angle_canvas.getAxis('left').setPen((255, 255, 255))
        self.height_canvas.getAxis('bottom').setPen((255, 255, 255)) 
        self.height_canvas.getAxis('left').setPen((255, 255, 255)) 

        # Force x-axes to stay aligned
        self.angle_canvas.setXLink(self.height_canvas)   
        layout.addWidget(self.plot_view)
        self.plot_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        export_button = QPushButton('Export')
        layout.addWidget(export_button)
        export_button.clicked.connect(self.export)

        X = []
        Y = []
        Y_M = []

        A = self.setup.trajectory.pos_i
        B = self.setup.trajectory.pos_f

        A.pos_loc(B)
        B.pos_loc(B)

        #######################
        # Plot nominal points
        #######################
        base_points = pg.ScatterPlotItem()
        for i in range(len(self.setup.fragmentation_point)):
            f_time = stn.times.fragmentation[i][0][0]
            X = self.setup.fragmentation_point[i].position.elev
            Y = f_time - nom_pick.time
            base_points.addPoints(x=[X], y=[Y], pen=(255, 0, 238), brush=(255, 0, 238), symbol='o')


        ptb_colors = [(0, 255, 26, 150), (3, 252, 176, 150), (252, 3, 3, 150), (176, 252, 3, 150), (255, 133, 3, 150),
                      (149, 0, 255, 150), (76, 128, 4, 150), (82, 27, 27, 150), (101, 128, 125, 150), (5, 176, 249, 150)]
        
        
        base_points.setZValue(1)
        self.height_canvas.addItem(base_points, update=True)

        #########################
        # Perturbation points
        #########################
        prt_points = pg.ScatterPlotItem()
        for i in range(len(self.setup.fragmentation_point)):
            data, remove = self.obtainPerts(stn.times.fragmentation, i)
            for pt in data:
                X = self.setup.fragmentation_point[i].position.elev
                Y = pt - nom_pick.time
                prt_points.addPoints(x=[X], y=[Y], pen=(255, 0, 238, 150), brush=(255, 0, 238, 150), symbol='o')

        self.height_canvas.addItem(prt_points, update=True)
 

        for pick in self.pick_list:
            if pick.group == 0:
                self.height_canvas.plot(x=[np.nanmin(X), np.nanmax(X)], y=[pick.time - nom_pick.time, pick.time - nom_pick.time], pen='g')
            else:
                self.height_canvas.plot(x=[np.nanmin(X), np.nanmax(X)], y=[pick.time - nom_pick.time, pick.time - nom_pick.time], pen='b')


        #####################
        # Angle Calculation
        #####################

        u = np.array([self.setup.trajectory.vector.x,
                      self.setup.trajectory.vector.y,
                      self.setup.trajectory.vector.z])

        angle_off = []
        X = []
        for i in range(len(self.setup.fragmentation_point)):
            az = stn.times.fragmentation[i][0][1]
            tf = stn.times.fragmentation[i][0][2]

            az = np.radians(az)
            tf = np.radians(180 - tf)
            v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

            angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))
            X.append(self.setup.fragmentation_point[i].position.elev)
        angle_off = np.array(angle_off)
        try:
            best_indx = np.nanargmin(abs(angle_off - 90))
            print("Optimal Ballistic Height {:.2f} km with angle of {:.2f} deg".format(X[best_indx]/1000, angle_off[best_indx]))

            self.angle_canvas.plot(x=[X[best_indx], X[best_indx]], y=[np.nanmin(np.append(angle_off, 90)), np.nanmax(np.append(angle_off, 90))], pen=(0, 0, 255))
            self.height_canvas.plot(x=[X[best_indx], X[best_indx]], y=[np.nanmin(Y), np.nanmax(Y)], pen=(0, 0, 255))
            self.angle_canvas.scatterPlot(x=X, y=angle_off, pen=(255, 255, 255), symbol='o', brush=(255, 255, 255))

            # self.angle_canvas.setXRange(10000, 45000, padding=0)
            best_arr = []
            angle_arr = []
            
        except ValueError:
            best_indx = None

        angle_off = 0
        height = None
        for i in range(len(self.setup.fragmentation_point)):
            az, _ = np.array(self.obtainPerts(stn.times.fragmentation, i, pt=1))
            tf, _ = np.array(self.obtainPerts(stn.times.fragmentation, i, pt=2))
            az = np.radians(az)
            tf = np.radians(180 - tf)

            for j in range(len(az)):
                v = np.array([np.sin(az[j])*np.sin(tf[j]), np.cos(az[j])*np.sin(tf[j]), -np.cos(tf[j])])

                angle_off_new = np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v)))))
                self.angle_canvas.scatterPlot(x=[self.setup.fragmentation_point[i].position.elev], y=[angle_off_new], symbol='o')

                if abs(angle_off_new - 90) < abs(angle_off - 90) and not np.isnan(angle_off_new):
                    angle_off = angle_off_new

                    height = self.setup.fragmentation_point[i].position.elev

        if height is not None:
            self.angle_canvas.addItem(pg.InfiniteLine(pos=(height, 0), angle=90, pen=QColor(0, 0, 255)))
            self.height_canvas.addItem(pg.InfiniteLine(pos=(height, 0), angle=90, pen=QColor(0, 0, 255)))
        
        self.angle_canvas.addItem(pg.InfiniteLine(pos=(0, 90), angle=0, pen=QColor(255, 0, 0)))

        # if self.setup.perturb:
        #     data, remove = chauvenet(X[best_arr])
        #     data_angle, remove_angle = chauvenet(angle_arr)
        #     print("Ballistic Perturbation Range: {:.5f} - {:.5f} km, with angles of {:.2f} - {:.2f} deg"\
        #                     .format(np.nanmin(data)/1000, np.nanmax(data)/1000, np.nanmin(data_angle), np.nanmax(data_angle)))
        #     print("Removed Points: {:} km".format(remove))
        #     print("Removed Angles: {:} deg".format(remove_angle))



        phigh = pg.PlotCurveItem([np.nanmin(X), np.nanmax(X)], [65, 65], pen = 'g')           
        plow = pg.PlotCurveItem([np.nanmin(X), np.nanmax(X)], [115, 115], pen = 'g')                  
        pfill = pg.FillBetweenItem(phigh, plow, brush = (0, 0, 255, 150))
        self.angle_canvas.addItem(phigh)
        self.angle_canvas.addItem(plow)
        self.angle_canvas.addItem(pfill)

        # self.height_canvas.plot(x=[20500, 20500], y=[-40, 100], pen=(255, 255, 255))
        # self.height_canvas.plot(x=[21500, 21500], y=[-40, 100], pen=(255, 255, 255))
        # self.height_canvas.plot(x=[25500, 25500], y=[-40, 100], pen=(255, 255, 255))
        # self.height_canvas.plot(x=[30400, 30400], y=[-40, 100], pen=(255, 255, 255))


        self.height_canvas.setTitle('Fragmentation Height Prediction of Given Pick', color=(0, 0, 0))
        self.angle_canvas.setTitle('Angles of Initial Acoustic Wave Path', color=(0, 0, 0))
        self.height_canvas.setLabel('left', 'Difference in Time from {:.2f}'.format(nom_pick.time), units='s')
        self.angle_canvas.setLabel('left', 'Angle Away from Trajectory of Initial Acoustic Wave Path', units='deg')
        self.height_canvas.setLabel('bottom', 'Height of Solution', units='m')
        self.angle_canvas.setLabel('bottom', 'Height of Solution', units='m')

        self.height_canvas.setLimits(xMin=0, xMax=50000, yMin=-40, yMax=100, minXRange=1000, maxXRange=33000, minYRange=2, maxYRange=140)

        self.setLayout(layout)


    def export(self):

        # set export parameters if needed
        #exporter.parameters()['width'] = 1000   # (note this also affects height parameter)
        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = dlg[0]

        exporter = pg.exporters.SVGExporter(self.plot_view.scene())

        file_name = file_name + '.svg'
        exporter.export(file_name)

    def obtainPerts(self, data, frag, pt=0):
        data_new = []

        for i in range(len(data[frag][1])):
            data_new.append(data[frag][1][i][pt])
        data, remove = chauvenet(data_new)

        return data, remove