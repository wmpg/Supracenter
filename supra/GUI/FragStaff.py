import numpy as np

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.Utils.AngleConv import chauvenet

class FragmentationStaff(QWidget):

    def __init__(self, setup, pack):

        QWidget.__init__(self)
        self.setWindowTitle('Fragmentation Staff')
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)
        
        # Take important values from main window class
        self.arrTimes, self.current_station, self.pick_list = pack

        nom_pick = self.pick_list[0]

        layout = QVBoxLayout()

        # Pass setup value
        self.setup = setup

        # Create plot
        self.plot_view = pg.GraphicsLayoutWidget()
        self.height_canvas = self.plot_view.addPlot()
        self.plot_view.nextRow()
        self.angle_canvas = self.plot_view.addPlot()

        # Force x-axes to stay aligned
        self.angle_canvas.setXLink(self.height_canvas)   
        layout.addWidget(self.plot_view)
        self.plot_view.sizeHint = lambda: pg.QtCore.QSize(100, 100)

        export_button = QPushButton('Export')
        layout.addWidget(export_button)
        export_button.clicked.connect(self.export)

        guess = True
        X = []
        Y = []
        Y_M = []

        A = self.setup.trajectory.pos_i
        B = self.setup.trajectory.pos_f

        A.pos_loc(B)
        B.pos_loc(B)

        # Get prediction of time of the meteor, so the timing of each fragmentation can be known
        length_of_meteor = np.sqrt((A.x - B.x)**2 + (A.y - B.y)**2 + (A.z - B.z)**2)
        time_of_meteor = length_of_meteor/self.setup.trajectory.v

        # Plot nominal points
        for i in range(len(self.setup.fragmentation_point)):
            X.append(self.setup.fragmentation_point[i].position.elev)
            Y.append(self.arrTimes[0, self.current_station, 1, i] - nom_pick.time)
            Y_M.append(self.arrTimes[0, self.current_station, 1, i] - nom_pick.time + time_of_meteor)

        X = np.array(X)
        Y_M = np.array(Y_M)
        Y = np.array(Y)
        ptb_colors = [(0, 255, 26, 150), (3, 252, 219, 150), (252, 3, 3, 150), (223, 252, 3, 150), (255, 133, 3, 150),
                      (149, 0, 255, 150), (76, 128, 4, 150), (82, 27, 27, 150), (101, 128, 125, 150), (255, 230, 249, 150)]
        base_points = pg.ScatterPlotItem()
        base_points.addPoints(x=X, y=Y_M + X/self.setup.trajectory.pos_i.elev*(Y - Y_M), pen=(255, 0, 238), brush=(255, 0, 238), symbol='o')
        for i in range(len(self.setup.fragmentation_point)):
            try:
                txt = pg.TextItem("{:.2f}".format(self.arrTimes[0, self.current_station, 4, i]), color=(255, 0, 238))
                txt.setPos(X[i], Y[i])
                self.height_canvas.addItem(txt)

            except:
                pass
        
        err = pg.ErrorBarItem(x=X, y=(Y + Y_M)/2, top=abs(Y_M - (Y + Y_M)/2), bottom=abs(Y_M - (Y + Y_M)/2), beam=0.5, pen=(255, 0, 238))
        self.height_canvas.addItem(err) 
        
        base_points.setZValue(1)
        self.height_canvas.addItem(base_points, update=True)

        # Perturbation points
        Y_P = [[]]*self.setup.perturb_times
        Y_P_M = [[]]*self.setup.perturb_times
        opt_points = []
        for ptb in range(self.setup.perturb_times):
            y_p = [] 
            y_p_M = []
            for i in range(len(self.setup.fragmentation_point)):
                y_p.append(self.arrTimes[ptb, self.current_station, 1, i] - nom_pick.time)
                y_p_M.append(self.arrTimes[ptb, self.current_station, 1, i] - nom_pick.time + time_of_meteor)

            y_p = np.array(y_p)
            y_p_M = np.array(y_p_M)

            err = pg.ErrorBarItem(x=X, y=(y_p + y_p_M)/2, top=abs(y_p_M-(y_p + y_p_M)/2), bottom=abs(y_p_M-(y_p + y_p_M)/2), beam=0.5, pen=ptb_colors[ptb])
            self.height_canvas.addItem(err)
            for i in range(len(self.setup.fragmentation_point)):
                try:
                    txt = pg.TextItem("{:.2f}".format(self.arrTimes[ptb, self.current_station, 4, i]), color=ptb_colors[ptb])
                    txt.setPos(X[i], y_p[i])
                    self.height_canvas.addItem(txt)
                except:
                    pass

               
            Y_P[ptb].append(y_p)
            Y_P_M[ptb].append(y_p_M)

            y_p = np.array(y_p)
            y_p_M = np.array(y_p_M)

            perturbation_points = pg.ScatterPlotItem()
            perturbation_points.addPoints(x=X, y=y_p_M + X/self.setup.trajectory.pos_i.elev*(y_p - y_p_M), pen=ptb_colors[ptb], brush=ptb_colors[ptb], symbol='o')
            perturbation_points.setZValue(0)
            self.height_canvas.addItem(perturbation_points, update=True)

            P_p_y = y_p_M + X/self.setup.trajectory.pos_i.elev*(y_p - y_p_M)
            idx = np.isfinite(X) & np.isfinite(P_p_y)
            if len(P_p_y[idx]) > 0:
                P_p = np.polyfit(X[idx], P_p_y[idx], 1)
                pert_line_y = P_p[0]*X + P_p[1]
                #self.height_canvas.plot(x=X, y=pert_line_y, pen=(0, 255, 26, 50))
                pert_opt = -P_p[1]/P_p[0]
                opt_points.append(pert_opt)
                #self.height_canvas.scatterPlot(x=[pert_opt], y=[0], pen=(0, 255, 26), symbol='+')
    
        # Error Bars


        idx = np.isfinite(X) & np.isfinite(Y)
        if len(Y[idx]) > 0 and len(X) > 0:
            P = np.polyfit(X[idx], Y[idx], 1)
            y = P[0]*np.array(X) + P[1]

            #self.height_canvas.plot(x=X, y=y, pen='b')
        
            optimal_height = -P[1]/P[0]

        idx = np.isfinite(X) & np.isfinite(Y)
        if len(Y[idx]) > 0 and len(X) > 0:
            P_M = np.polyfit(X[idx], Y[idx] + time_of_meteor, 1)
            y_M = P_M[0]*np.array(X) + P_M[1]
            #self.height_canvas.plot(x=X, y=y_M, pen='b')
            optimal_height_M = -P_M[1]/P_M[0]

        P_o_y = Y_M + X/self.setup.trajectory.pos_i.elev*(Y - Y_M)
        idx = np.isfinite(X) & np.isfinite(P_o_y)
        if len(P_o_y[idx]) > 0 and len(X) > 0:
            P_o = np.polyfit(X[idx], P_o_y[idx], 1)
            y_o = P_o[0]*X + P_o[1]
        #self.height_canvas.plot(X, y_o, pen=(255, 0, 238))
        try:
            best_guess_frac = optimal_height_M*self.setup.trajectory.pos_i.elev/(self.setup.trajectory.pos_i.elev - optimal_height + optimal_height_M)
        except:
            guess = False

        if guess:
            opt_points.append(best_guess_frac)
            opt_points = np.array(opt_points)

        for pick in self.pick_list:
            if pick.group == 0:
                self.height_canvas.plot(x=[np.nanmin(X), np.nanmax(X)], y=[pick.time - nom_pick.time, pick.time - nom_pick.time], pen='g')
            else:
                self.height_canvas.plot(x=[np.nanmin(X), np.nanmax(X)], y=[pick.time - nom_pick.time, pick.time - nom_pick.time], pen='b')

        if guess:
            #self.height_canvas.scatterPlot(x=[-P_o[1]/P_o[0]], y=[0], pen=(255, 0, 238), symbol='+')
            print("Optimal Solution: {:.5f} km".format(opt_points[-1]/1000))
            if self.setup.perturb:
                data, remove = chauvenet(opt_points)
                print("Perturbation Range: {:.5f} - {:.5f} km".format(np.nanmin(data)/1000, np.nanmax(data)/1000))
                print("Removed Points: {:} km".format(remove))
        
        u = np.array([self.setup.trajectory.vector.x,
                      self.setup.trajectory.vector.y,
                      self.setup.trajectory.vector.z])

        angle_off = []
        for i in range(len(self.setup.fragmentation_point)):
            az = self.arrTimes[0, self.current_station, 2, i]
            tf = self.arrTimes[0, self.current_station, 3, i]
            # az = angle2NDE(az)
            az = np.radians(az)
            tf = np.radians(180 - tf)
            v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

            angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))

        angle_off = np.array(angle_off)
        try:
            best_indx = np.nanargmin(abs(angle_off - 90))
            print("Optimal Ballistic Height {:.2f} km with angle of {:.2f} deg".format(X[best_indx]/1000, angle_off[best_indx]))

            self.angle_canvas.plot(x=[X[best_indx], X[best_indx]], y=[np.nanmin(np.append(angle_off, 90)), np.nanmax(np.append(angle_off, 90))], pen=(0, 0, 255))
            self.height_canvas.plot(x=[X[best_indx], X[best_indx]], y=[np.nanmin(Y), np.nanmax(Y)], pen=(0, 0, 255))
            self.angle_canvas.scatterPlot(x=X, y=angle_off, pen=(255, 255, 255), symbol='o', brush=(255, 255, 255))
            self.angle_canvas.plot(x=[np.nanmin(X), np.nanmax(X)], y=[90, 90], pen='r')
            self.angle_canvas.setXRange(15000, 45000, padding=0)
            best_arr = []
            angle_arr = []
            for ptb in range(self.setup.perturb_times):
                angle_off = []
                for i in range(len(self.setup.fragmentation_point)):
                    az = self.arrTimes[ptb, self.current_station, 2, i]
                    tf = self.arrTimes[ptb, self.current_station, 3, i]
                    # az = angle2NDE(az)
                    az = np.radians(az)
                    tf = np.radians(180 - tf)
                    v = np.array([np.sin(az)*np.sin(tf), np.cos(az)*np.sin(tf), -np.cos(tf)])

                    angle_off.append(np.degrees(np.arccos(np.dot(u/np.sqrt(u.dot(u)), v/np.sqrt(v.dot(v))))))

                angle_off = np.array(angle_off)

                best_indx = (np.nanargmin(abs(angle_off - 90)))
                best_arr.append(best_indx)
                angle_arr.append(angle_off[best_indx])
                self.angle_canvas.plot(x=[X[best_indx], X[best_indx]], y=[np.nanmin(np.append(angle_off, 90)), np.nanmax(np.append(angle_off, 90))], pen=(0, 0, 255, 100))
                self.height_canvas.plot(x=[X[best_indx], X[best_indx]], y=[np.nanmin(Y), np.nanmax(Y)], pen=(0, 0, 255, 100))
                self.angle_canvas.scatterPlot(x=X, y=angle_off, pen=ptb_colors[ptb], symbol='o', brush=ptb_colors[ptb])
                self.angle_canvas.plot(x=[np.nanmin(X), np.nanmax(X)], y=[90, 90], pen='r')

            if self.setup.perturb:
                data, remove = chauvenet(X[best_arr])
                data_angle, remove_angle = chauvenet(angle_arr)
                print("Ballistic Perturbation Range: {:.5f} - {:.5f} km, with angles of {:.2f} - {:.2f} deg"\
                                .format(np.nanmin(data)/1000, np.nanmax(data)/1000, np.nanmin(data_angle), np.nanmax(data_angle)))
                print("Removed Points: {:} km".format(remove))
                print("Removed Angles: {:} deg".format(remove_angle))
        except ValueError:
            best_indx = None


        # self.height_canvas.plot(x=[20500, 20500], y=[-40, 100], pen=(255, 255, 255))
        # self.height_canvas.plot(x=[21500, 21500], y=[-40, 100], pen=(255, 255, 255))
        # self.height_canvas.plot(x=[25500, 25500], y=[-40, 100], pen=(255, 255, 255))
        # self.height_canvas.plot(x=[30400, 30400], y=[-40, 100], pen=(255, 255, 255))


        self.height_canvas.setTitle('Fragmentation Height Prediction of Given Pick')
        self.angle_canvas.setTitle('Angles of Initial Acoustic Wave Path')
        self.height_canvas.setLabel('left', 'Difference in Time from {:.2f}'.format(nom_pick.time), units='s')
        self.angle_canvas.setLabel('left', 'Angle Away from Trajectory of Initial Acoustic Wave Path', units='deg')
        self.height_canvas.setLabel('bottom', 'Height of Solution', units='m')
        self.angle_canvas.setLabel('bottom', 'Height of Solution', units='m')

        self.height_canvas.setLimits(xMin=17000, xMax=50000, yMin=-40, yMax=100, minXRange=1000, maxXRange=33000, minYRange=2, maxYRange=140)

        self.setLayout(layout)

    def export(self):

        # set export parameters if needed
        #exporter.parameters()['width'] = 1000   # (note this also affects height parameter)
        dlg = QFileDialog.getSaveFileName(self, 'Save File')

        file_name = dlg[0]

        exporter = pg.exporters.SVGExporter(self.plot_view.scene())

        file_name = file_name + '.svg'
        exporter.export(file_name)