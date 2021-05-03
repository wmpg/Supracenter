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

        self.buildGUI()
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

        self.dots_x = []
        self.dots_y = []

        #######################
        # Plot nominal points
        #######################
        base_points = pg.ScatterPlotItem()
        for i in range(len(self.setup.fragmentation_point)):

            f_time = stn.times.fragmentation[i][0][0]

            X = self.setup.fragmentation_point[i].position.elev
            Y = f_time - nom_pick.time
            self.dots_x.append(X)
            self.dots_y.append(Y)

            base_points.addPoints(x=[X], y=[Y], pen=(255, 0, 238), brush=(255, 0, 238), symbol='o')


        ptb_colors = [(0, 255, 26, 150), (3, 252, 176, 150), (252, 3, 3, 150), (176, 252, 3, 150), (255, 133, 3, 150),
                      (149, 0, 255, 150), (76, 128, 4, 150), (82, 27, 27, 150), (101, 128, 125, 150), (5, 176, 249, 150)]
        
        
        base_points.setZValue(1)
        self.height_canvas.addItem(base_points, update=True)

        #########################
        # Plot Precursor Points
        #########################

        pre_points = pg.ScatterPlotItem()
        
        for i in range(len(self.setup.fragmentation_point)):
        
            X = self.setup.fragmentation_point[i].position.elev
            
            # Distance between frag point and the ground below it
            v_dist = X - stn.metadata.position.elev
            
            # Horizontal distance betweent the new ground point and the stn
            h_dist = self.setup.fragmentation_point[i].position.ground_distance(stn.metadata.position)
            
            # Speed of wave in air
            v_time = v_dist/310

            # Speed of wave in ground
            h_time = h_dist/3100

            # Total travel time
            Y = v_time + h_time - nom_pick.time

            pre_points.addPoints(x=[X], y=[Y], pen=(210, 235, 52), brush=(210, 235, 52), symbol='o')

        self.height_canvas.addItem(pre_points, update=True)

        #########################
        # Perturbation points
        #########################
        prt_points = pg.ScatterPlotItem()
        for i in range(len(self.setup.fragmentation_point)):
            data, remove = self.obtainPerts(stn.times.fragmentation, i)
            Y = []
            X = self.setup.fragmentation_point[i].position.elev
            for pt in data:
                
                Y = (pt - nom_pick.time)

                self.dots_x.append(X)
                self.dots_y.append(Y)
                prt_points.addPoints(x=[X], y=[Y], pen=(255, 0, 238, 150), brush=(255, 0, 238, 150), symbol='o')

            self.height_canvas.addItem(prt_points, update=True)

        for pick in self.pick_list:
            if pick.group == 0:
                self.height_canvas.addItem(pg.InfiniteLine(pos=(0, pick.time - nom_pick.time), angle=0, pen=QColor(0, 255, 0)))
            else:
                self.height_canvas.addItem(pg.InfiniteLine(pos=(0, pick.time - nom_pick.time), angle=0, pen=QColor(0, 0, 255)))

        self.dots = np.array([self.dots_x, self.dots_y])

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

        ###############################
        # Find optimal ballistic angle
        ###############################
        try:
            best_indx = np.nanargmin(abs(angle_off - 90))
            print("Optimal Ballistic Height {:.2f} km with angle of {:.2f} deg".format(X[best_indx]/1000, angle_off[best_indx]))
            self.angle_canvas.addItem(pg.InfiniteLine(pos=(X[best_indx], 0), angle=90, pen=QColor(0, 0, 255)))
            self.height_canvas.addItem(pg.InfiniteLine(pos=(X[best_indx], 0), angle=90, pen=QColor(0, 0, 255)))
            self.angle_canvas.scatterPlot(x=X, y=angle_off, pen=(255, 255, 255), symbol='o', brush=(255, 255, 255))

            best_arr = []
            angle_arr = []
            
        except ValueError:
            best_indx = None

        angle_off = 0
        height = None
        for i in range(len(self.setup.fragmentation_point)):
            for j in range(len(stn.times.fragmentation[i][1])):
                az = stn.times.fragmentation[i][1][j][1]
                tf = stn.times.fragmentation[i][1][j][2]
                az = np.radians(az)
                tf = np.radians(180 - tf)
                v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

                angle_off_new = np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v)))))
                self.angle_canvas.scatterPlot(x=[self.setup.fragmentation_point[i].position.elev], y=[angle_off_new], symbol='o')

                if abs(angle_off_new - 90) < abs(angle_off - 90) and not np.isnan(angle_off_new):
                    angle_off = angle_off_new

                    height = self.setup.fragmentation_point[i].position.elev

        if height is not None:
            self.angle_canvas.addItem(pg.InfiniteLine(pos=(height, 0), angle=90, pen=QColor(0, 0, 255)))
            self.height_canvas.addItem(pg.InfiniteLine(pos=(height, 0), angle=90, pen=QColor(0, 0, 255)))
        
        self.angle_canvas.addItem(pg.InfiniteLine(pos=(0, 90), angle=0, pen=QColor(255, 0, 0)))

        #####################
        # Build plot window
        #####################

        # 25 deg tolerance window
        phigh = pg.PlotCurveItem([np.nanmin(X), np.nanmax(X)], [65, 65], pen = 'g')           
        plow = pg.PlotCurveItem([np.nanmin(X), np.nanmax(X)], [115, 115], pen = 'g')                  
        pfill = pg.FillBetweenItem(phigh, plow, brush = (0, 0, 255, 100))
        self.angle_canvas.addItem(phigh)
        self.angle_canvas.addItem(plow)
        self.angle_canvas.addItem(pfill)

        # Build axes
        self.height_canvas.setTitle('Fragmentation Height Prediction of Given Pick', color=(0, 0, 0))
        self.angle_canvas.setTitle('Angles of Initial Acoustic Wave Path', color=(0, 0, 0))
        self.height_canvas.setLabel('left', 'Difference in Time from {:.2f}'.format(nom_pick.time), units='s')
        self.angle_canvas.setLabel('left', 'Angle Away from Trajectory of Initial Acoustic Wave Path', units='deg')
        self.height_canvas.setLabel('bottom', 'Height of Solution', units='m')
        self.angle_canvas.setLabel('bottom', 'Height of Solution', units='m')

        # self.height_canvas.setLimits(xMin=B.elev, xMax=A.elev, yMin=-40, yMax=100, minXRange=1000, maxXRange=33000, minYRange=2, maxYRange=140)

        # Fonts
        font= QFont()
        font.setPixelSize(20)
        self.height_canvas.getAxis("bottom").tickFont = font
        self.height_canvas.getAxis("left").tickFont = font
        self.height_canvas.getAxis('bottom').setPen((255, 255, 255)) 
        self.height_canvas.getAxis('left').setPen((255, 255, 255))
        self.angle_canvas.getAxis("bottom").tickFont = font
        self.angle_canvas.getAxis("left").tickFont = font
        self.angle_canvas.getAxis('bottom').setPen((255, 255, 255)) 
        self.angle_canvas.getAxis('left').setPen((255, 255, 255))
        self.setLayout(layout)

    def buildGUI(self):
        self.setWindowTitle('Fragmentation Staff')
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)
        
        theme(self)

    def export(self):

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