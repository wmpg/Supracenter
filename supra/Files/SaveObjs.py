import pickle
import os

from supra.Utils.TryObj import *

from supra.Atmosphere.AtmosClass import *
from supra.Utils.Classes import Config
import numpy as np
from supra.GUI.Tools.GUITools import *

class BAMFile:

    def __init__(self):
        
        self.file_name = None
        self.atmos = Atmos()
        self.stn_list = [] 
        self.setup = Config()

    def __str__(self):

        A = not (self.file_name is None)
        B = not (len(self.stn_list) == 0)
        C = not (self.atmos.isEmpty())
        D = not (self.setup.isEmpty())


        return "BAM Object: \n File_Name:    {:} \n Station_list: {:} \n Atmos:        {:} \n Setup:        {:} \n".format(\
                    str(A), str(B), str(C), str(D))



class Prefs:

    def __init__(self):
        self.workdir = os.path.expanduser("~")
        self.avg_sp_sound = 330
        self.contour_res = 25
        self.debug = True
        self.ballistic_en = False
        self.frag_en = False
        self.recalc_times = False
        self.wind_en = False
        self.pert_en = False
        self.atm_type = 'none'
        self.pert_type = 'none'
        self.pert_num = 0
        self.pso_debug = False
        self.pso_theta      = 90
        self.pso_phi        = 90
        self.pso_min_ang    = 330.0
        self.pso_min_dist   = 3000.0
        self.pso_max_iter   = 100
        self.pso_swarm_size = 100
        self.pso_run_times  = 1
        self.pso_min_error  = 1e-8
        self.pso_min_step   = 1e-8
        self.pso_phi_p      = 0.5
        self.pso_phi_g      = 0.5
        self.pso_omega      = 0.5

    def save(self, obj):

        self.workdir = obj.workdir.text()
        self.avg_sp_sound = tryFloat(obj.avg_sp_sound.text())
        self.contour_res = tryFloat(obj.contour_res.text())
        self.debug = obj.debug.isChecked()
        self.ballistic_en = obj.ballistic_en.isChecked()
        self.frag_en = obj.frag_en.isChecked()
        self.recalc_times = obj.recalc_times.isChecked()
        self.wind_en = obj.wind_en.isChecked()
        self.pert_en = obj.pert_en.isChecked()
        self.atm_type = obj.atm_type.currentText()
        self.pert_type = obj.pert_type.currentText()
        self.pert_num = tryInt(obj.pert_num.text())
        self.pso_debug = obj.pso_debug.isChecked()
        self.pso_theta      = tryInt(obj.pso_theta.text())
        self.pso_phi        = tryInt(obj.pso_phi.text())
        self.pso_min_ang    = tryFloat(obj.pso_min_ang.text())
        self.pso_min_dist   = tryFloat(obj.pso_min_dist.text())
        self.pso_max_iter   = tryInt(obj.pso_max_iter.text())
        self.pso_swarm_size = tryInt(obj.pso_swarm_size.text())
        self.pso_run_times  = tryInt(obj.pso_run_times.text())
        self.pso_min_error  = tryFloat(obj.pso_min_error.text())
        self.pso_min_step   = tryFloat(obj.pso_min_step.text())
        self.pso_phi_p      = tryFloat(obj.pso_phi_p.text())
        self.pso_phi_g      = tryFloat(obj.pso_phi_g.text())
        self.pso_omega      = tryFloat(obj.pso_omega.text())

        with open(os.path.join('supra', 'Misc', 'BAMprefs.bam'), 'wb') as f:
            pickle.dump(self, f)

        print('STATUS: Saved Preferences')

    def load(self):
        try:
            with open(os.path.join('supra', 'Misc', 'BAMprefs.bam'), 'rb') as f:
                return pickle.load(f)
        except EOFError as e:
            return None

class Stats:
    def __init__(self):
        pass

class Atmos:
    def __init__(self, avg_sp_sound=310):

        self.avg_sp_sound = avg_sp_sound
        self.default_weather = np.array([[    0.0, self.avg_sp_sound, 0.0, 0.0],
                                         [    0.0, self.avg_sp_sound, 0.0, 0.0],
                                         [99999.0, self.avg_sp_sound, 0.0, 0.0]])

    def __str__(self):
        
        return "Atmosphere class"

    def isEmpty(self):

        return not hasattr(self, 'ecmwf')

    def addECMWF(self, lat, lon, rng, time, file_name):
        
        self.ecmwf = ECMWF(lat, lon, rng, time, file_name)

    def addECMWFSpread(self, lat, lon, rng, time, file_name):
        
        self.ecmwf.spread_file = self.ecmwf.spread(lat, lon, rng, time, file_name)

    def addRadio(self):
        pass

    def loadSounding(self, file_name, weather_type, lat=0, lon=0, rng=0, time=0):
        
        if weather_type == 'ecmwf':
            self.addECMWF(lat, lon, rng, time, file_name)
        elif weather_type == 'spread':
            self.addECMWFSpread(lat, lon, rng, time, file_name)
        elif weather_type == 'radio':
            print('Not done yet')
        else:
            print('Unrecognized weather type')

    def getSounding(self, lat=0, lon=0, heights=[100000, 0], spline=100):
        ''' lat = [start lat, end lat], etc
        '''
        prefs = Prefs()
        prefs = prefs.load()
        self.avg_sp_sound = 310
        
        if prefs.atm_type == 'ecmwf':
            if hasattr(self, 'ecmwf'):
                return self.ecmwf.getProfile(lat, lon, heights, prefs, spline=spline)
            else:
                raise Exception("ECMWF is not loaded into the BAM file, trying to export default weather")
                self.default_weather = np.array([[heights[1], self.avg_sp_sound, 0.0, 0.0],
                                                 [heights[1], self.avg_sp_sound, 0.0, 0.0],
                                                 [heights[0], self.avg_sp_sound, 0.0, 0.0]])
                return self.default_weather, None
        elif prefs.atm_type == 'radio':
            self.default_weather = np.array([[heights[1], self.avg_sp_sound, 0.0, 0.0],
                                 [heights[1], self.avg_sp_sound, 0.0, 0.0],
                                 [heights[0], self.avg_sp_sound, 0.0, 0.0]])
            return self.default_weather, None
        elif prefs.atm_type == 'none':

            self.default_weather = np.array([[heights[1], self.avg_sp_sound, 0.0, 0.0],
                                             [heights[1], self.avg_sp_sound, 0.0, 0.0],
                                             [heights[0], self.avg_sp_sound, 0.0, 0.0]])

            return self.default_weather, None

        else:
            print('Unrecognized weather type')
            return self.default_weather, None


