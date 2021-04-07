import os
import numpy as np
import pickle

from functools import partial

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import pyqtgraph as pg

from supra.GUI.Tools.Theme import theme
from supra.GUI.Tools.GUITools import *
from supra.GUI.Tools.ReportHelper import *

from supra.Files.SaveObjs import Prefs

class PreferenceWindow(QWidget):

    def __init__(self):

        QWidget.__init__(self)
        
        self.buildGUI()
        self.loadPref()

    def buildGUI(self):
        self.setWindowTitle('Program Wide Preferences')
        app_icon = QtGui.QIcon()
        app_icon.addFile(os.path.join('supra', 'GUI', 'Images', 'preferences.png'), QtCore.QSize(16, 16))
        self.setWindowIcon(app_icon)

        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)
        self.setPalette(p)

        theme(self)

        layout = QGridLayout()
        self.setLayout(layout)
     
        tab_widget = QTabWidget()
        tab_widget.blockSignals(True)
        tab_widget.blockSignals(False)
        layout.addWidget(tab_widget, 1, 1, 1, 2)

        general_tab, general_tab_layout = createTab(tab_widget, 'General')
        _, self.workdir = createLabelEditObj('Working Directory', general_tab_layout, 0, width=1, h_shift=0, tool_tip='')
        _, self.avg_sp_sound = createLabelEditObj('Average Speed of Sound', general_tab_layout, 1, width=1, h_shift=0, tool_tip='', validate='float', default_txt='310')
        _, self.contour_res = createLabelEditObj('Contour Resolution', general_tab_layout, 2, width=1, h_shift=0, tool_tip='', validate='float', default_txt='25')
        self.debug = createToggle('Debug', general_tab_layout, 3, width=1, h_shift=1, tool_tip='')
        self.ballistic_en = createToggle('Show Ballistic', general_tab_layout, 4, width=1, h_shift=1, tool_tip='')
        self.frag_en = createToggle('Show Fragmentation', general_tab_layout, 5, width=1, h_shift=1, tool_tip='')
        self.recalc_times = createToggle('Recalculate Times', general_tab_layout, 6, width=1, h_shift=1, tool_tip='')
        self.recalc_sigs = createToggle('Recalculate Signals', general_tab_layout, 7, width=1, h_shift=1, tool_tip='')

        browse_button = QPushButton('Browse')
        general_tab_layout.addWidget(browse_button, 0, 3)
        browse_button.clicked.connect(partial(folderSearch, self.workdir))

        atmos_tab, atmos_tab_layout = createTab(tab_widget, 'Atmospheric')
        self.wind_en = createToggle('Enable Winds', atmos_tab_layout, 0, width=1, h_shift=1, tool_tip='')
        self.pert_en = createToggle('Enable Perturbations', atmos_tab_layout, 1, width=1, h_shift=1, tool_tip='')
        _, self.atm_type = createComboBoxObj('Atmosphere Type', atmos_tab_layout, 2, items=['none', 'ecmwf'], width=1, h_shift=0, tool_tip='')
        _, self.pert_type = createComboBoxObj('Perturbation Type', atmos_tab_layout, 3, items=['none', 'spread'], width=1, h_shift=0, tool_tip='')
        _, self.pert_num = createLabelEditObj('Perturb Times', atmos_tab_layout, 4, width=1, h_shift=0, tool_tip='', validate='int', default_txt='1')

        pso_tab, pso_tab_layout = createTab(tab_widget, 'PSO Search')
        self.pso_debug = createToggle('PSO Debug', pso_tab_layout, 0, width=1, h_shift=1, tool_tip='')
        _, self.pso_theta = createLabelEditObj('Theta Resolution', pso_tab_layout, 1, width=1, h_shift=0, tool_tip='', validate='int', default_txt='90')
        _, self.pso_phi = createLabelEditObj('Phi Resolution', pso_tab_layout, 2, width=1, h_shift=0, tool_tip='', validate='int', default_txt='90')
        _, self.pso_min_ang = createLabelEditObj('Horizontal Tolerance', pso_tab_layout, 3, width=1, h_shift=0, tool_tip='', validate='float', default_txt='1e-5')
        _, self.pso_min_dist = createLabelEditObj('Vertical Tolerance', pso_tab_layout, 4, width=1, h_shift=0, tool_tip='', validate='float', default_txt='1000.0')
        _, self.pso_max_iter = createLabelEditObj('Maximum Iterations', pso_tab_layout, 5, width=1, h_shift=0, tool_tip='', validate='int', default_txt='100')
        _, self.pso_swarm_size = createLabelEditObj('Swarm Size', pso_tab_layout, 6, width=1, h_shift=0, tool_tip='', validate='int', default_txt='100')
        _, self.pso_run_times = createLabelEditObj('Run Times', pso_tab_layout, 7, width=1, h_shift=0, tool_tip='', validate='int', default_txt='1')
        _, self.pso_min_error = createLabelEditObj('Minimum Error', pso_tab_layout, 8, width=1, h_shift=0, tool_tip='', validate='float', default_txt='1e-8')
        _, self.pso_min_step = createLabelEditObj('Minimum Step', pso_tab_layout, 9, width=1, h_shift=0, tool_tip='', validate='float', default_txt='1e-8')
        _, self.pso_phi_p = createLabelEditObj('Phi_p', pso_tab_layout, 10, width=1, h_shift=0, tool_tip='', validate='float', default_txt='0.5')
        _, self.pso_phi_g = createLabelEditObj('Phi_g', pso_tab_layout, 11, width=1, h_shift=0, tool_tip='', validate='float', default_txt='0.5')
        _, self.pso_omega = createLabelEditObj('Omega', pso_tab_layout, 12, width=1, h_shift=0, tool_tip='', validate='float', default_txt='0.5')

        save_button = QPushButton('Save')
        layout.addWidget(save_button, 2, 2)
        save_button.clicked.connect(self.savePref)

    def savePref(self):

        pref_obj = Prefs()
        pref_obj.save(self)
        file_name = os.path.join('supra', 'Misc', 'BAMprefs.bam')
        errorMessage('Preferences Saved!', 0, title='Saved!', info='File saved to {:} ({:})'.format(file_name, byteify(os.stat(file_name).st_size)))
        self.close()

    def loadPref(self):
        try:
            with open(os.path.join('supra', 'Misc', 'BAMprefs.bam'), 'rb') as f:
                prefs = pickle.load(f)
        except EOFError as e:
            # Unable to load preferences
            # Read a default prefs file
            pass

        self.workdir.setText(prefs.workdir)       
        self.avg_sp_sound.setText(str(prefs.avg_sp_sound))  
        self.contour_res.setText(str(prefs.contour_res))   
        self.debug.setChecked(prefs.debug)         
        self.ballistic_en.setChecked(prefs.ballistic_en)  
        self.recalc_times.setChecked(prefs.recalc_times)
        try:
            self.recalc_sigs.setChecked(prefs.recalc_sigs)  
        except AttributeError:
            self.recalc_sigs.setChecked(False) 

        self.frag_en.setChecked(prefs.frag_en)       
        self.wind_en.setChecked(prefs.wind_en)       
        self.pert_en.setChecked(prefs.pert_en)       
        comboSet(self.atm_type, prefs.atm_type)      
        comboSet(self.pert_type, prefs.pert_type)     
        self.pert_num.setText(str(prefs.pert_num))      
        self.pso_debug.setChecked(prefs.pso_debug)     
        self.pso_theta.setText(str(prefs.pso_theta))     
        self.pso_phi.setText(str(prefs.pso_phi))       
        self.pso_min_ang.setText(str(prefs.pso_min_ang))   
        self.pso_min_dist.setText(str(prefs.pso_min_dist))  
        self.pso_max_iter.setText(str(prefs.pso_max_iter))  
        self.pso_swarm_size.setText(str(prefs.pso_swarm_size))
        self.pso_run_times.setText(str(prefs.pso_run_times)) 
        self.pso_min_error.setText(str(prefs.pso_min_error)) 
        self.pso_min_step.setText(str(prefs.pso_min_step))  
        self.pso_phi_p.setText(str(prefs.pso_phi_p))     
        self.pso_phi_g.setText(str(prefs.pso_phi_g))     
        self.pso_omega.setText(str(prefs.pso_omega))     