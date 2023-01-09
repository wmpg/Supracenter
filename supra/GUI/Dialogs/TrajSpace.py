import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from functools import partial

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *


from supra.GUI.Tools.CustomWidgets import MatplotlibPyQT
from supra.Stations.ProcessStation import *
from supra.Utils.Formatting import *
from supra.Geminus.geminusSearch import periodSearch, presSearch




def hypfunc(x, a, b, h, k):
    return b*np.sqrt(1 + ((x - h)/a)**2) + k
def invhypfunc(x, a, b, h, k):
    result = np.abs(a*np.sqrt(((x - k)/b)**2 - 1) + h)
    return result



class TrajSpace(QWidget):

    def __init__(self, bam):

        QWidget.__init__(self)
        
        self.height_points = []
        self.bam = bam
        self.buildGUI()
        self.calculate()


    def clearax(self):
        self.height_points = []
        # self.pvh_graph.ax1.clear()
        self.pvh_graph.ax2.clear()
        self.pvh_graph.ax3.clear()
        # self.pvh_graph.ax4.clear()
        self.pvh_graph.ax5.clear()
        self.pvh_graph.ax6.clear()

    def binify(self):
        bin_size = float(self.bin_edits.text())
        h_min = float(self.min_height_edits.text())*1000
        h_max = float(self.max_height_edits.text())*1000

        bins = np.arange(h_min, h_max + bin_size, bin_size)
        bin_content = [0]*len(bins)

        pts = self.height_points
        
        for pt in pts:
            h = pt[0]
            a = pt[1]

            indx = [n for n,i in enumerate(bins) if i >= h ][0] - 1
            bin_content[indx] += a

        # self.pvh_graph.ax4.scatter(bins + bin_size/2, bin_content, alpha=1.0)

    def calculate(self):

        h_min = float(self.min_height_edits.text())*1000
        h_max = float(self.max_height_edits.text())*1000
        
        self.geminus_heights = []
        self.geminus_p = []
        self.geminus_t = []
        self.geminus_stat = []

        self.clearax()

        infra_list = self.getInfraStats()

        popt_list = self.genHyperbola(infra_list)
        
        trace_list, resp_list = self.getInfraTraces(infra_list)

        N = int(self.N_edits.text())
        l = float(self.l_edits.text())


        for tt, (trace, resp, popt, infra) in enumerate(zip(trace_list, resp_list, popt_list, infra_list)):



            if popt is None:
                continue

            stn_name = "{:}-{:}".format(infra.metadata.network, infra.metadata.code)

            a, t = procTrace(trace, ref_datetime=self.bam.setup.fireball_datetime, resp=resp, bandpass=None, backup=False)
            a_new = [item for sublist in a for item in sublist]
            t_new = [item for sublist in t for item in sublist]

            a = np.array(a_new)
            t = np.array(t_new)
            # p_list, h_list = self.convertTimes(a[0], t[0], popt)
        
            heights = invhypfunc(t, *popt)
            h_indicies = np.where(np.logical_and(heights>=h_min, heights<=h_max))

            divide = 0
            for i in range(len(h_indicies[0])):

                if h_indicies[0][i] - h_indicies[0][i-1] != 1:
                    divide = i

            branch_1 = h_indicies[0][:divide]
            branch_2 = h_indicies[0][divide:-1]


            FASes = []
            logs = []
            L = len(a)
            spacer = int(L/(l*(1-N) + N))
            shifter = int(spacer*(1-l))
            print("{:} Window Length = {:.2f} s".format(stn_name, spacer/trace.stats.sampling_rate))
            print("{:} Window Shift = {:.2f} s".format(stn_name, shifter/trace.stats.sampling_rate))
            for i in range((L - spacer)//shifter + 1):
                a_temp, t_temp = procTrace(trace, ref_datetime=self.bam.setup.fireball_datetime, resp=resp, bandpass=None, backup=False)
                # if self.h_space_tog.isChecked():
                #     self.pvh_graph.ax1.plot(invhypfunc(t_temp[0][i*shifter:int(i*shifter + spacer)], *popt), a_temp[0][i*shifter:int(i*shifter + spacer)], alpha=0.5)
                # else:
                #     self.pvh_graph.ax1.plot(t_temp[0][i*shifter:int(i*shifter + spacer)], a_temp[0][i*shifter:int(i*shifter + spacer)], alpha=0.5)

                freq, FAS = genFFT(a[i*shifter:int(i*shifter + spacer)], trace.stats.sampling_rate)
                if i == 0:
                    FAS_N = FAS

                
                FASes.append(np.sum(FAS))
                logs.append(FAS/FAS_N)

            best_fas = np.argmax(FASes)
            #second best fas
            FASes[best_fas] = 0
            bad_fas = np.argmax(FASes)
            self.pvh_graph.ax2.plot(freq, logs[best_fas], label=stn_name)

            if self.auto_gain.isChecked():
                j = np.where(logs[best_fas] >= logs[bad_fas])
            else:
                j = np.where(logs[best_fas] >= float(self.gain_edits.text()))

            # plt.loglog(freq[j], logs[np.argmax(FASes)][j])

            filtered_freq = freq[j]

            # print("Optimal Frequency Range {:.2f} - {:.2f} Hz".format(filtered_freq[0], filtered_freq[-1]))

            if self.stat_bandpass.isChecked():
                bnps = infra.bandpass
            else:
                bnps = [filtered_freq[0], filtered_freq[-1]]

            if bnps is not None:
                print("Optimal Frequency Range {:.2f} - {:.2f} Hz".format(bnps[0], bnps[-1]))

            a, t = procTrace(trace, ref_datetime=self.bam.setup.fireball_datetime, resp=resp, bandpass=bnps, backup=False)



            a_new = [item for sublist in a for item in sublist]
            t_new = [item for sublist in t for item in sublist]

            a = np.array(a_new)
            t = np.array(t_new)


            s2n = np.max(a)/np.median(np.abs(a))
            filtered_wave = a

            if tt == 0:
                total_vals = []
                total_elements = []
                total_h = []


            if self.h_space_tog.isChecked():
                h = invhypfunc(t, *popt)
                if self.branchselector.isChecked():
                    h = h[branch_2]
                    vals = a[branch_2]
                else:
                    h = h[branch_1]
                    vals = a[branch_1]
                self.pvh_graph.ax3.plot(h, vals, alpha=0.3, label="{:}".format(stn_name))
                
                # hil = reHilbert(vals)

                # self.height_points.append([h, vals])
                # # for ampl, heig in zip(vals, h):
                # #     self.height_points.append([heig, ampl])


                # # Last element
                # if tt == len(trace_list) - 1:


                #     # do the hilbert thing here

                #     self.binify()
                    
                # self.pvh_graph.ax3.plot(h, vals, alpha=0.3, label="{:}: Optimal Bandpass ({:.2f} - {:.2f} Hz) S/N {:.2f}".format(stn_name, filtered_freq[0], filtered_freq[-1], s2n))
            else:


                # hil = reHilbert(a)
                # total_vals.append(hil)
                # total_elements.append(len(hil))

                # # Last element

                # if tt == len(trace_list) - 1:
                #     min_trace = np.nanmin(total_elements)

                #     for vv, val_cut in enumerate(total_vals):
                #         if vv == 0:
                #             adjusted_cut = val_cut[:min_trace]
                #         else:
                #             adjusted_cut += val_cut[:min_trace]

                #     self.pvh_graph.ax4.plot(t[:len(adjusted_cut)], adjusted_cut[:len(t)], alpha=1.0)
                # self.pvh_graph.ax3.plot(t[0], a[0], alpha=0.3, label="{:}: Optimal Bandpass ({:.2f} - {:.2f} Hz) S/N {:.2f}".format(stn_name, filtered_freq[0], filtered_freq[-1], s2n))
                self.pvh_graph.ax3.plot(t, a, alpha=0.3, label="{:}".format(stn_name))


            # sig = a[0][j]
            # sig_times = t[0][j]
            # p, freq, FAS = findDominantPeriodPSD(sig, trace.stats.sampling_rate, normalize=False)

            # self.pvh_graph.ax4.plot(freq, FAS, label="Dominant Period of Signal: {:.2f} s".format(p))
            # print("Dominant Period of Signal: {:.2f} s".format(p))
            
            shortest_period = 1/trace.stats.sampling_rate
            longest_period = t[shifter] - t[0]
            ax5h = []
            ax5p = []
            ax6h = []
            ax6p = []
            for i in range((L - spacer)//shifter + 1):

                max_p = np.nanmax(a[i*shifter:int(i*shifter + spacer)])
                min_p = np.nanmin(a[i*shifter:int(i*shifter + spacer)])

                height_of_sol = invhypfunc(t[int(i*shifter + spacer//2)], *popt)

                if self.h_space_tog.isChecked():
                    h = invhypfunc(t[int(i*shifter + spacer//2)], *popt)

                    if h in heights[branch_1] and not self.branchselector.isChecked():
                        ax5h.append(h)
                        ax5p.append((max_p - min_p)/2)
                        
                    elif h in heights[branch_2] and self.branchselector.isChecked():
                        ax5h.append(h)
                        ax5p.append((max_p - min_p)/2)
                else:
                    ax5h.append(t[int(i*shifter + spacer//2)])
                    ax5p.append((max_p - min_p)/2)

                p, freq, FAS = findDominantPeriodPSD(a[i*shifter:int(i*shifter + spacer)], trace.stats.sampling_rate, normalize=False)
                
                # if h_min <= height_of_sol and height_of_sol <= h_max:


                if self.h_space_tog.isChecked():

                    h = invhypfunc(t[int(i*shifter + spacer//2)], *popt)

                    if h in heights[branch_1] and not self.branchselector.isChecked():
                        self.geminus_heights.append(height_of_sol)
                        self.geminus_p.append((max_p - min_p)/2)
                        self.geminus_stat.append(infra)
                        self.geminus_t.append(p)
                        ax6h.append(h)
                        ax6p.append(p)
                        
                    elif h in heights[branch_2] and self.branchselector.isChecked():
                        self.geminus_heights.append(height_of_sol)
                        self.geminus_p.append((max_p - min_p)/2)
                        self.geminus_stat.append(infra)
                        self.geminus_t.append(p) 
                        ax6h.append(h)
                        ax6p.append(p)
                else:
                    ax6h.append(t[int(i*shifter + spacer//2)])
                    ax6p.append(p)

                self.pvh_graph.ax6.axhline(y=shortest_period, linestyle='-')

            self.pvh_graph.ax5.scatter(ax5h, ax5p, label=stn_name)
            self.pvh_graph.ax6.scatter(ax6h, ax6p, label=stn_name)



            t_in_range = t[branch_1]
            try:
                t_min_range_1 = t_in_range[0]
                t_max_range_1 = t_in_range[-1]
            except IndexError:
                t_min_range_1, t_max_range_1 = np.nan, np.nan
            t_in_range = t[branch_2]
            try:
                t_min_range_2 = t_in_range[0]
                t_max_range_2 = t_in_range[-1]
            except IndexError:
                t_min_range_2, t_max_range_2 = np.nan, np.nan

        # plt.axhline(y=longest_period, color='k', linestyle='-')
        # spacer = 5*60//2
        # shifter = 100

        # period_list = []
        # period_h_list = []

        # for pp in range(len(waveform_list)):
        #     periods = []
        #     period_times = []
        #     for ii in range(len(waveform_list[pp][0])//shifter):


        #         if shifter*ii - spacer >= 0 and shifter*ii + spacer <= len(waveform_list[pp][0]):
        #             temp_waveform = waveform_list[pp][0][shifter*ii-spacer:shifter*ii+spacer]
        #             temp_time = time_list[pp][0][shifter*ii-spacer:shifter*ii+spacer]
        #             try:
        #                 st = infra_list[pp].stream.select(channel="*DF")[0]
        #             except IndexError:
        #                 sf = infra_list[pp].stream.select(channel="*HZ")[0]
        #             period, freq, FAS = findDominantPeriodPSD(st, resp=infra_list[pp].response)
        #             plt.semilogx(freq, FAS)

        #             periods.append(period)
        #             period_times.append(time_list[pp][0][shifter*ii])

        #     if popt_list[pp] is None:
        #         continue
        #     min_height = 17000
        #     max_height = 40000
        #     height_periods = invhypfunc(period_times, *popt_list[pp])
        #     h_indicies = np.where(np.logical_and(height_periods>=min_height, height_periods<=max_height))

        #     new_heights = []
        #     new_period = []
        #     for hh in h_indicies[0]:

        #         new_heights.append(height_periods[hh])
        #         new_period.append(periods[hh])


        #     period_list.append(new_period)
        #     period_h_list.append(new_heights)


        
        # plt.show()

        # # This gives the wrong station labels

        # for pp in range(len(p_list)):

        #     self.pvh_graph.ax1.plot(h_list[pp], np.abs(p_list[pp]), label="{:}".format(infra_list[pp].metadata.code), alpha=0.5)
        # #     self.pvh_graph.ax2.plot(period_h_list[pp], period_list[pp], label="{:}".format(infra_list[pp].metadata.code), alpha=0.5)
        
        try:
            t_min = float(self.min_time_edits.text())
            t_max = float(self.max_time_edits.text())
        except:
            t_min = t[0]
            t_max = t[-1]


        # if self.h_space_tog.isChecked():
        #     self.pvh_graph.ax1.set_xlabel("Height [m]")
        #     self.pvh_graph.ax1.set_xlim([h_min, h_max])
        # else:
        #     self.pvh_graph.ax1.set_xlabel("Time [s]")
        #     self.pvh_graph.ax1.set_xlim([t_min, t_max])

        # self.pvh_graph.ax1.set_ylabel("Overpressure [Pa]")

        self.pvh_graph.ax2.set_xlabel("Frequency [Hz]")
        self.pvh_graph.ax2.set_ylabel("Gain")
        self.pvh_graph.ax2.set_xscale('log')
        self.pvh_graph.ax2.set_yscale('log')
        self.pvh_graph.ax2.axvline(x=filtered_freq[0], linestyle='-')
        self.pvh_graph.ax2.axvline(x=filtered_freq[-1], linestyle='-')

        if self.h_space_tog.isChecked():
            self.pvh_graph.ax3.set_xlabel("Height [m]")
            self.pvh_graph.ax3.set_xlim([h_min, h_max])
        else:
            self.pvh_graph.ax3.set_xlabel("Time [s]")
            self.pvh_graph.ax3.set_xlim([t_min, t_max])
            if self.branchselector.isChecked():
                self.pvh_graph.ax3.axvline(x=t_min_range_2, linestyle='-')
                self.pvh_graph.ax3.axvline(x=t_max_range_2, linestyle='-')
            else:
                self.pvh_graph.ax3.axvline(x=t_min_range_1, linestyle='-')
                self.pvh_graph.ax3.axvline(x=t_max_range_1, linestyle='-')
        self.pvh_graph.ax3.set_ylabel("Overpressure [Pa]")
        
        # if self.h_space_tog.isChecked():
        #     self.pvh_graph.ax4.set_xlabel("Height [m]")
        #     self.pvh_graph.ax4.set_xlim([h_min, h_max])
        # else:
        #     self.pvh_graph.ax4.set_xlabel("Time [s]")
        #     self.pvh_graph.ax4.set_xlim([t_min, t_max])
        # self.pvh_graph.ax4.set_ylabel("Overpressure [Pa]")        

        # self.pvh_graph.ax4.set_xlabel("Frequency [Hz]")
        # self.pvh_graph.ax4.set_ylabel("Gain")
        # self.pvh_graph.ax4.set_xscale('log')

        if self.h_space_tog.isChecked():
            self.pvh_graph.ax5.set_xlabel("Height [m]")
            self.pvh_graph.ax5.set_xlim([h_min, h_max])
        else:
            self.pvh_graph.ax5.set_xlabel("Time [s]")
            self.pvh_graph.ax5.set_xlim([t_min, t_max])
            if self.branchselector.isChecked():
                self.pvh_graph.ax5.axvline(x=t_min_range_2, linestyle='-')
                self.pvh_graph.ax5.axvline(x=t_max_range_2, linestyle='-')
            else:
                self.pvh_graph.ax5.axvline(x=t_min_range_1, linestyle='-')
                self.pvh_graph.ax5.axvline(x=t_max_range_1, linestyle='-')
        self.pvh_graph.ax5.set_ylabel("Overpressure [Pa]")

        if self.h_space_tog.isChecked():
            self.pvh_graph.ax6.set_xlabel("Height [m]")
            self.pvh_graph.ax6.set_xlim([h_min, h_max])
        else:
            self.pvh_graph.ax6.set_xlabel("Time [s]")
            self.pvh_graph.ax6.set_xlim([t_min, t_max])
            if self.branchselector.isChecked():
                self.pvh_graph.ax6.axvline(x=t_min_range_2, linestyle='-')
                self.pvh_graph.ax6.axvline(x=t_max_range_2, linestyle='-')
            else:
                self.pvh_graph.ax6.axvline(x=t_min_range_1, linestyle='-')
                self.pvh_graph.ax6.axvline(x=t_max_range_1, linestyle='-')
        self.pvh_graph.ax6.set_ylabel("Dominant Period [s]")
        self.pvh_graph.ax2.legend()
        # self.pvh_graph.ax3.legend()
        # self.pvh_graph.ax5.legend()
        # self.pvh_graph.ax6.legend()
        self.pvh_graph.show()

    def buildGUI(self):
        self.setWindowTitle('Traj Space')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'BAM_no_wave.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)

        self.pvh_graph = MatplotlibPyQT()
        # self.pvh_graph.ax1 = self.pvh_graph.figure.add_subplot(611)
        self.pvh_graph.ax2 = self.pvh_graph.figure.add_subplot(411)
        self.pvh_graph.ax3 = self.pvh_graph.figure.add_subplot(412)
        # self.pvh_graph.ax4 = self.pvh_graph.figure.add_subplot(513)
        self.pvh_graph.ax5 = self.pvh_graph.figure.add_subplot(413)
        self.pvh_graph.ax6 = self.pvh_graph.figure.add_subplot(414)
        layout.addWidget(self.pvh_graph, 1, 1, 100, 1)
        export_raw = createButton("Export Overpressures and Periods", layout, 101, 1, self.exportraw, args=[])


        run_button = createButton("Make Plots", layout, 1, 3, self.calculate, args=[])
        _, self.N_edits = createLabelEditObj("Number of Windows", layout, 2, width=1, h_shift=1, tool_tip='', validate='int', default_txt='100')
        _, self.l_edits = createLabelEditObj("Percent Overlap", layout, 3, width=1, h_shift=1, tool_tip='', validate='float', default_txt='0.5')
        self.h_space_tog = createToggle("Height Space", layout, 4, width=1, h_shift=2, tool_tip='')
        self.h_space_tog.setChecked(True)
        self.auto_gain = createToggle("Auto Gain Limits", layout, 5, width=1, h_shift=2, tool_tip='')
        _, self.gain_edits = createLabelEditObj("Gain Cutoff", layout, 6, width=1, h_shift=1, tool_tip='', validate='float', default_txt='5')
        self.auto_gain.setChecked(True)
        ro_button = createButton("Calculate Relaxation Radii", layout, 7, 3, self.geminusify, args=[])
        _, self.min_height_edits = createLabelEditObj("Minimum Height [km]", layout, 8, width=1, h_shift=1, tool_tip='', validate='float', default_txt='17')
        _, self.max_height_edits = createLabelEditObj("Maximum Height [km]", layout, 9, width=1, h_shift=1, tool_tip='', validate='float', default_txt='40')
        self.branchselector = createToggle("Use Branch 2", layout, 10, width=1, h_shift=2, tool_tip='')
        self.branchselector.setChecked(True)
        _, self.min_time_edits = createLabelEditObj("Minimum Time [s]", layout, 11, width=1, h_shift=1, tool_tip='', validate='float', default_txt='300')
        _, self.max_time_edits = createLabelEditObj("Maximum Time [s]", layout, 12, width=1, h_shift=1, tool_tip='', validate='float', default_txt='600')
        self.stat_bandpass = createToggle("Use Station Bandpass", layout, 13, width=1, h_shift=2, tool_tip='')
        _, self.bin_edits = createLabelEditObj("Size of Bins [m]", layout, 14, width=1, h_shift=1, tool_tip='', validate='float', default_txt='100')
        self.bin_edits.editingFinished.connect(self.binify)

        self.ro_graph = MatplotlibPyQT()
        self.ro_graph.ax = self.ro_graph.figure.add_subplot(111)
        layout.addWidget(self.ro_graph, 15, 2, 90, 2)
        export_ro = createButton("Export Relaxation Radii Curve", layout, 105, 3, self.exportro, args=[])

    def exportraw(self):
            
        file_name = saveFile('csv')

        with open(file_name, 'w+') as f:
            f.write('Station, Height [m], Overpressure [Pa], Dominant Period [s] \n')
            for i in range(len(self.geminus_heights)):
                f.write('{:}, {:}, {:}, {:} \n'.format(self.geminus_stat[i].metadata.code, self.geminus_heights[i], self.geminus_p[i], self.geminus_t[i]))

        errorMessage('Raw Data Saved!', 0, info='Saved to File {:}'.format(file_name))
    
    def exportro(self):
            
        file_name = saveFile('csv')

        with open(file_name, 'w+') as f:
            f.write('Station, Height [m], Overpressure [Pa], Dominant Period [s], Ro (Weak-Shock, Overpressure) [m], Ro (Linear, Overpressure) [m], Ro (Weak-Shock, Dominant Period) [m], Ro (Linear, Dominant Period) [m]\n')
            for i in range(len(self.geminus_heights)):
                f.write('{:}, {:}, {:}, {:}, {:}, {:}, {:}, {:} \n'.format(\
                    self.geminus_stat[i].metadata.code, self.geminus_heights[i], self.geminus_p[i], self.geminus_t[i], \
                    self.ro_data[i, 1], self.ro_data[i, 2], self.ro_data[i, 3], self.ro_data[i, 4]))

        errorMessage('Data saved to CSV!', 0, info='Saved to File {:}'.format(file_name))    
            

    def geminusify(self):
        data = []
        max_steps = len(self.geminus_heights)
        for hh, h in enumerate(self.geminus_heights):
            traj = self.bam.setup.trajectory

            source = traj.findGeo(h)

            self.sounding_pres = None
            source_list = [source.lat, source.lon, source.elev/1000]


            stat = self.geminus_stat[hh]
            stat_pos = stat.metadata.position
            stat = [stat_pos.lat, stat_pos.lon, stat_pos.elev/1000]

            v = traj.v/1000

            theta = 90 - traj.zenith.deg

            dphi = np.degrees(np.arctan2(stat_pos.lon - source.lon, stat_pos.lat - source.lat)) - traj.azimuth.deg

            # Switch 3 doesn't do anything in this version of overpressure.py
            sw = [1, 0, 1]
            
            lat =   [source.lat, stat_pos.lat]
            lon =   [source.lon, stat_pos.lon]
            elev =  [source.elev, stat_pos.elev]

            sounding, _ = self.bam.atmos.getSounding(lat=lat, lon=lon, heights=elev, ref_time=self.bam.setup.fireball_datetime)

            pres = 10*101.325*np.exp(-0.00012*sounding[:, 0])
            
            sounding_pres = np.zeros((sounding.shape[0], sounding.shape[1]+1))
            sounding_pres[:, :-1] = sounding
            sounding_pres[:, -1] = pres
            sounding_pres[:, 1] -= 273.15

            sounding_pres = np.flip(sounding_pres, axis=0)

            gem_inputs = [source_list, stat, v, theta, dphi, sounding_pres, sw, True, True]
            Ro_ws_p, Ro_lin_p = presSearch(self.geminus_p[hh], gem_inputs, paths=False)
            Ro_ws_t, Ro_lin_t = periodSearch(self.geminus_t[hh], gem_inputs, paths=False)

            data.append([h, Ro_ws_p, Ro_lin_p, Ro_ws_t, Ro_lin_t])
            print("Complete Step {:} of {:}".format(hh + 1, max_steps))
        for row in data:
            self.ro_graph.ax.scatter(row[0], row[1], c='m', label="Weak-Shock Period")
            self.ro_graph.ax.scatter(row[0], row[2], c='c', label="Linear Period")
            self.ro_graph.ax.scatter(row[0], row[3], c='y', label="Weak-Shock Overpressure")
            self.ro_graph.ax.scatter(row[0], row[4], c='g', label="Linear Overpressure")

        self.ro_data = np.array(data)
        self.ro_graph.ax.set_xlabel("Height [m]")
        self.ro_graph.ax.set_ylabel("Relaxation Radius [m]")
        self.ro_graph.ax.legend()
        self.ro_graph.show()
    def getInfraStats(self):    
        
        infra_list = []

        for stn in self.bam.stn_list:
            if len(stn.stream.select(channel="*DF")) > 0:
                infra_list.append(stn)
            # elif len(stn.stream.select(channel="*HZ")) > 0:
            #     infra_list.append(stn)
        return infra_list
        
    def genHyperbola(self, infra_list):

        popt_list = []
        for stn in infra_list:
            points_x = []
            points_y = []


            for i in range(len(self.bam.setup.fragmentation_point)):

                f_time = stn.times.fragmentation[i][0][0]

                points_x.append(self.bam.setup.fragmentation_point[i].position.elev)
                points_y.append(f_time)

            ### ADD PERTURBATION POINTS HERE

            try:
                x_vals = []
                y_vals = []

                for pp in range(len(points_y)):
                    if not np.isnan(points_y[pp]):

                        x_vals.append(points_x[pp])
                        y_vals.append(points_y[pp])

                x_vals, y_vals = np.array(x_vals), np.array(y_vals)



                # Hyperbola in the form y = kx (since we invert the x values)
                from scipy.optimize import curve_fit

                popt, pcov = curve_fit(hypfunc, x_vals, y_vals)
                popt_list.append(popt)
            except TypeError:
                print("Could not generate hyperbola fit!")
                popt_list.append(None)
            except RuntimeError:
                print("Optimal hyperbola not found!")
                popt_list.append(None)
            except ValueError:
                print("No Arrivals!")
                popt_list.append(None)
        return popt_list

    def getInfraTraces(self, infra_list):
        resp_list = []
        trace_list = []
        for stn in infra_list:
            st, resp, gap_times = procStream(stn, ref_time=self.bam.setup.fireball_datetime)
            resp_list.append(resp)
            temp_st = findChn(st, "*DF")
            # if len(temp_st) == 0:
            #     temp_st = findChn(st, "*HZ")
            trace_list.append(temp_st)
        return trace_list, resp_list
            # waveform_data, time_data = procTrace(st, ref_datetime=self.bam.setup.fireball_datetime,\
            #         resp=resp, bandpass=None)

            # wave_d = []
            # time_d = []
            # for wave, time in zip(waveform_data, time_data):
            #     wave_d.append(wave)
            #     time_d.append(time)

            # wave_list.append(wave_d)
            # time_list.append(time_d)

        # return wave_list, time_list

    def convertTimes(self, wave, time, popt):

        press_list = []
        height_list = []

        if popt is None:
            return [], []
        
        height_peaks = invhypfunc(time, *popt)

        ##################
        # Get Bounds for Heights
        ##################
        h_min = float(self.min_height_edits.text())*1000
        h_max = float(self.max_height_edits.text())*1000
        h_indicies = np.where(np.logical_and(height_peaks>=h_min, height_peaks<=h_max))
        new_heights = []
        new_press = []

        for hh in h_indicies[0]:

            new_heights.append(height_peaks[hh])
            new_press.append(wave[hh])

        press_list.append(new_press)
        height_list.append(new_heights)

        return press_list, height_list

    def genPFFAS(self, wave_list, time_list):


        a = wave_list
        t = time_list

        ### optimal bandpass is when S/N is a maximum
        sampling_rate = 1/(t[1] - t[0])

        nyq = sampling_rate/2
        low_freq = 1/(t[-1] - t[0])


        p, freq, FAS = findDominantPeriodPSD(a, sampling_rate, normalize=True)

        return p, freq, FAS

    def getOptBandpass(self, trace_list, infra_list):

        b_list = [[1e-3, 9.9]]

        for bandpass in b_list:
            for trace, stn in zip(trace_list, infra_list):
                waveform_data, time_data = procTrace(trace, ref_datetime=self.bam.setup.fireball_datetime,\
                   resp=stn.response, bandpass=None)

                for ii in range(100):
                    L = len(waveform_data[0])//100

                    p, freq, FAS = self.genPFFAS(waveform_data[0][ii*L: (ii+1)*L], time_data[0][ii*L: (ii+1)*L])
                    if ii == 0:
                        FAS_n = FAS

                    f2 = interp1d(freq, FAS)
                    xnew = np.logspace(-3, np.log10(9))
                    plt.semilogx(xnew, f2(xnew))
                    plt.semilogx(freq, FAS)

        plt.show()