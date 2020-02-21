import pickle

from PyQt5.QtCore import *

from supra.Utils.TryObj import *
from supra.Utils.Classes import Position, Station

from supra.GUI.GUITools import *

def saveGUI(self, write=True, autosave=False):

    """if write == True - used for saving data to setup obj, but not writing
    """
    if write and not autosave:
        dlg = QFileDialog.getSaveFileName(self, 'Save File')

    self.setup.fireball_name = self.fireball_name_edits.text()
    self.setup.get_data = tryBool(self.get_data_edits.currentText())
    self.setup.run_mode = self.run_mode_edits.currentText()
    self.setup.debug = tryBool(self.debug_edits.currentText())

    self.setup.working_directory = self.working_directory_edits.text()
    self.setup.arrival_times_file = self.arrival_times_edits.text()
    self.setup.sounding_file = self.sounding_file_edits.text()
    self.setup.perturbation_spread_file = self.perturbation_file_edits.text()
    self.setup.station_picks_file = self.station_picks_edits.text()
    self.setup.replot_points_file = self.points_name_edits.text()

    self.setup.lat_centre = tryFloat(self.lat_centre_edits.text())
    self.setup.lon_centre = tryFloat(self.lon_centre_edits.text())
    self.setup.deg_radius = tryFloat(self.deg_radius_edits.text())
    self.setup.fireball_datetime = self.fireball_datetime_edits.dateTime().toPyDateTime()
    self.setup.v_sound = tryFloat(self.v_sound_edits.text())

    self.setup.t0 = tryFloat(self.t0_edits.text())
    self.setup.v = tryFloat(self.v_edits.text())
    self.setup.azimuth = tryAngle(self.azim_edits.text())
    self.setup.zenith = tryAngle(self.zangle_edits.text())
    self.setup.lat_i = tryFloat(self.lat_i_edits.text())
    self.setup.lon_i = tryFloat(self.lon_i_edits.text())
    self.setup.elev_i = tryFloat(self.elev_i_edits.text())
    self.setup.lat_f = tryFloat(self.lat_f_edits.text())
    self.setup.lon_f = tryFloat(self.lon_f_edits.text())
    self.setup.elev_f = tryFloat(self.elev_f_edits.text())

    self.setup.pos_i = tryPosition(self.setup.lat_i, self.setup.lon_i, self.setup.elev_i)
    self.setup.pos_f = tryPosition(self.setup.lat_f, self.setup.lon_f, self.setup.elev_f)

    self.setup.trajectory = tryTrajectory(self.setup.t0, self.setup.v, self.setup.azimuth, self.setup.zenith, self.setup.pos_i, self.setup.pos_f)

    self.setup.show_ballistic_waveform = tryBool(self.show_ballistic_waveform_edits.currentText())

    self.setup.fragmentation_point = trySupracenter(fromTable(self.fragmentation_point))
    self.setup.show_fragmentation_waveform = tryBool(self.show_fragmentation_waveform_edits.currentText())

    frag_lat  = tryFloat(self.lat_frag_edits.text())
    frag_lon  = tryFloat(self.lon_frag_edits.text())
    frag_elev = tryFloat(self.elev_frag_edits.text())
    frag_time = tryFloat(self.time_frag_edits.text())
    self.setup.manual_fragmentation_search = trySupracenter(tryPosition(frag_lat, frag_lon, frag_elev), frag_time)

    self.setup.v_fixed = tryFloat(self.v_fixed_edits.text())
    self.setup.enable_restricted_time = self.restricted_time_check.isChecked()
    self.setup.restricted_time = self.restricted_time_edits.dateTime().toPyDateTime()

    self.setup.azimuth_min = tryAngle(self.azimuth_min_edits.text())
    self.setup.azimuth_max = tryAngle(self.azimuth_max_edits.text())
    self.setup.zenith_min = tryAngle(self.zangle_min_edits.text())
    self.setup.zenith_max = tryAngle(self.zangle_max_edits.text())
    self.setup.lat_min = tryFloat(self.lat_min_edits.text())
    self.setup.lat_max = tryFloat(self.lat_max_edits.text())
    self.setup.lon_min = tryFloat(self.lon_min_edits.text())
    self.setup.lon_max = tryFloat(self.lon_max_edits.text())
    self.setup.elev_min = tryFloat(self.elev_min_edits.text())
    self.setup.elev_max = tryFloat(self.elev_max_edits.text())
    self.setup.t_min = tryFloat(self.t_min_edits.text())
    self.setup.t_max = tryFloat(self.t_max_edits.text())
    self.setup.v_min = tryFloat(self.v_min_edits.text())
    self.setup.v_max = tryFloat(self.v_max_edits.text())
    self.setup.weight_distance_min = tryFloat(self.weight_distance_min_edits.text())
    self.setup.weight_distance_max = tryFloat(self.weight_distance_max_edits.text())

    self.setup.pos_min = tryPosition(self.setup.lat_min, self.setup.lon_min, self.setup.elev_min)
    self.setup.pos_max = tryPosition(self.setup.lat_max, self.setup.lon_max, self.setup.elev_max)

    self.setup.supra_min = trySupracenter(self.setup.pos_min, self.setup.t_min)
    self.setup.supra_max = trySupracenter(self.setup.pos_max, self.setup.t_max)

    self.setup.traj_min = tryTrajectory(self.setup.t_min, self.setup.v_min, self.setup.azimuth_min, \
                                     self.setup.zenith_min, self.setup.pos_min, self.setup.pos_min)
    self.setup.traj_max = tryTrajectory(self.setup.t_max, self.setup.v_max, self.setup.azimuth_max, \
                                     self.setup.zenith_max, self.setup.pos_max, self.setup.pos_max)

    self.setup.enable_winds = tryBool(self.enable_winds_edits.currentText())
    self.setup.weather_type = self.weather_type_edits.currentText()

    self.setup.perturb_times = tryInt(self.perturb_times_edits.text())
    self.setup.observe_frag_no = tryInt(self.frag_no_edits.text())
    self.setup.perturb = tryBool(self.perturb_edits.currentText())
    self.setup.perturb_method = self.perturb_method_edits.currentText()

    self.setup.n_theta = tryInt(self.n_theta_edits.text())
    self.setup.n_phi = tryInt(self.n_phi_edits.text())
    self.setup.h_tol = tryFloat(self.h_tol_edits.text())
    self.setup.v_tol = tryFloat(self.v_tol_edits.text())

    self.setup.maxiter = tryInt(self.maxiter_edits.text())
    self.setup.swarmsize = tryInt(self.swarmsize_edits.text())
    self.setup.run_times = tryInt(self.run_times_edits.text())
    self.setup.minfunc = tryFloat(self.minfunc_edits.text())
    self.setup.minstep = tryFloat(self.minstep_edits.text())
    self.setup.phip = tryFloat(self.phip_edits.text())
    self.setup.phig = tryFloat(self.phig_edits.text())
    self.setup.omega = tryFloat(self.omega_edits.text())
    self.setup.pso_debug = tryBool(self.pso_debug_edits.currentText())

    self.setup.contour_res = tryInt(self.contour_res_edits.text())
    self.setup.high_f = self.high_f_edits.text()
    self.setup.high_b = self.high_b_edits.text()
    self.setup.rm_stat = self.rm_stat_edits.text()

    self.setup.stations = []
    temp_stats = fromTable(self.extra_point)
    for stn in temp_stats:
        try:
            self.setup.stations.append(Station(stn[0], stn[1], Position(stn[2], stn[3], stn[4]), stn[5], stn[6], stn[7]))
        except TypeError:

            print("WARNING: Station: {:}-{:} could not be read".format(stn[0], stn[1]))
    
    # Check if it has the file name to autosave to
    if not hasattr(self, 'setup_file'):
        dlg = QFileDialog.getSaveFileName(self, 'Save File') 
        autosave = False

    if autosave:    
        with open(self.setup_file, 'wb') as f:
            pickle.dump(self.setup, f)

    else:
        if write:
            if '.pkl' not in dlg[0]:
                output = dlg[0] + '.pkl'
            else:
                output = dlg[0]

            with open(output, 'wb') as f:
                pickle.dump(self.setup, f)
                self.setup_file = dlg[0]

    self.setup = saveDefaults(self.setup)

    errorMessage("Setup saved!", 0, title='Saved!')


def loadGUI(self):
    dlg = QFileDialog()
    dlg.setFileMode(QFileDialog.AnyFile)
    dlg.setNameFilters(['Pickle File (*.pkl)'])
    #filenames = QStringList()

    dlg.exec_()

    filename = dlg.selectedFiles()

    if len(filename) == 0:
        errorMessage("No File Selected", 1)
        return None

    if '.pkl' not in filename[0]:
        errorMessage("No File Selected!", 1)
        return None
    
    self.setup_file = filename[0]

    try:
        with open(filename[0], 'rb') as f:
            self.setup = pickle.load(f)
    except EOFError as e:
        errorMessage('Unable to read file', 2, info='file: {:}'.format(filename[0]), detail='{:}'.format(e))
        return None

    self.setup = saveDefaults(self.setup)

    self.fireball_name_edits.setText(self.setup.fireball_name)
    comboSet(self.get_data_edits, self.setup.get_data)
    comboSet(self.run_mode_edits, self.setup.run_mode)
    comboSet(self.debug_edits, self.setup.debug)

    self.working_directory_edits.setText(self.setup.working_directory)
    self.arrival_times_edits.setText(self.setup.arrival_times_file)
    self.sounding_file_edits.setText(self.setup.sounding_file)
    self.perturbation_file_edits.setText(self.setup.perturbation_spread_file)
    self.station_picks_edits.setText(self.setup.station_picks_file)
    self.points_name_edits.setText(self.setup.replot_points_file)

    self.lat_centre_edits.setText(str(self.setup.lat_centre))
    self.lon_centre_edits.setText(str(self.setup.lon_centre))
    self.deg_radius_edits.setText(str(self.setup.deg_radius))
    self.fireball_datetime_edits.setDateTime(self.setup.fireball_datetime)
    self.v_sound_edits.setText(str(self.setup.v_sound))

    self.t0_edits.setText(str(self.setup.trajectory.t))
    self.v_edits.setText(str(self.setup.trajectory.v))
    self.azim_edits.setText(str(self.setup.trajectory.azimuth.deg))
    self.zangle_edits.setText(str(self.setup.trajectory.zenith.deg))
    self.lat_i_edits.setText(str(self.setup.trajectory.pos_i.lat))
    self.lon_i_edits.setText(str(self.setup.trajectory.pos_i.lon))
    self.elev_i_edits.setText(str(self.setup.trajectory.pos_i.elev))
    self.lat_f_edits.setText(str(self.setup.trajectory.pos_f.lat))
    self.lon_f_edits.setText(str(self.setup.trajectory.pos_f.lon))
    self.elev_f_edits.setText(str(self.setup.trajectory.pos_f.elev))

    comboSet(self.show_ballistic_waveform_edits, self.setup.show_ballistic_waveform)

    frag_list = []
    for element in self.setup.fragmentation_point:
        frag_list.append(element.toList())

    toTable(self.fragmentation_point, frag_list)
    comboSet(self.show_fragmentation_waveform_edits, self.setup.show_fragmentation_waveform)

    if self.setup.manual_fragmentation_search == None:
        self.lat_frag_edits.setText('')
        self.lon_frag_edits.setText('')
        self.elev_frag_edits.setText('')
        self.time_frag_edits.setText('')
    else: 
        self.lat_frag_edits.setText(str(self.setup.manual_fragmentation_search[0].position.lat))
        self.lon_frag_edits.setText(str(self.setup.manual_fragmentation_search[0].position.lon))
        self.elev_frag_edits.setText(str(self.setup.manual_fragmentation_search[0].position.elev))
        self.time_frag_edits.setText(str(self.setup.manual_fragmentation_search[0].time))

    self.v_fixed_edits.setText(str(self.setup.v_fixed))
    self.restricted_time_check.setChecked(self.setup.enable_restricted_time)
    self.restricted_time_edits.setDateTime(self.setup.restricted_time)

    self.azimuth_min_edits.setText(str(self.setup.traj_min.azimuth.deg))
    self.azimuth_max_edits.setText(str(self.setup.traj_max.azimuth.deg))
    self.zangle_min_edits.setText(str(self.setup.traj_min.zenith.deg))
    self.zangle_max_edits.setText(str(self.setup.traj_max.zenith.deg))
    self.lat_min_edits.setText(str(self.setup.traj_min.pos_f.lat))
    self.lat_max_edits.setText(str(self.setup.traj_max.pos_f.lat))
    self.lon_min_edits.setText(str(self.setup.traj_min.pos_f.lon))
    self.lon_max_edits.setText(str(self.setup.traj_max.pos_f.lon))
    self.elev_min_edits.setText(str(self.setup.traj_min.pos_f.elev))
    self.elev_max_edits.setText(str(self.setup.traj_max.pos_f.elev))
    self.t_min_edits.setText(str(self.setup.traj_min.t))
    self.t_max_edits.setText(str(self.setup.traj_max.t))
    self.v_min_edits.setText(str(self.setup.traj_min.v))
    self.v_max_edits.setText(str(self.setup.traj_max.v))
    self.weight_distance_min_edits.setText(str(self.setup.weight_distance_min))
    self.weight_distance_max_edits.setText(str(self.setup.weight_distance_max))

    comboSet(self.enable_winds_edits, self.setup.enable_winds)
    comboSet(self.weather_type_edits, self.setup.weather_type)

    self.perturb_times_edits.setText(str(self.setup.perturb_times))
    self.frag_no_edits.setText(str(self.setup.observe_frag_no))
    comboSet(self.perturb_edits, self.setup.perturb)
    comboSet(self.perturb_method_edits, self.setup.perturb_method)

    self.n_theta_edits.setText(str(self.setup.n_theta))
    self.n_phi_edits.setText(str(self.setup.n_phi))
    self.h_tol_edits.setText(str(self.setup.h_tol))
    self.v_tol_edits.setText(str(self.setup.v_tol))

    self.maxiter_edits.setText(str(self.setup.maxiter))
    self.swarmsize_edits.setText(str(self.setup.swarmsize))
    self.run_times_edits.setText(str(self.setup.run_times))
    self.minfunc_edits.setText(str(self.setup.minfunc))
    self.minstep_edits.setText(str(self.setup.minstep))
    self.phip_edits.setText(str(self.setup.phip))
    self.phig_edits.setText(str(self.setup.phig))
    self.omega_edits.setText(str(self.setup.omega))
    comboSet(self.pso_debug_edits, self.setup.pso_debug)

    self.contour_res_edits.setText(str(self.setup.contour_res))
    self.high_f_edits.setText(str(self.setup.high_f))
    self.high_b_edits.setText(str(self.setup.high_b))
    self.rm_stat_edits.setText(str(self.setup.rm_stat))

    toTableFromStn(self.extra_point, self.setup.stations)
