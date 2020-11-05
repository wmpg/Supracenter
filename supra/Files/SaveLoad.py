import pickle

from PyQt5.QtCore import *

from supra.Files.SaveObjs import BAMFile, Atmos
from supra.Utils.Formatting import *
from supra.Utils.TryObj import *
from supra.Utils.Classes import *
from supra.GUI.Tools.GUITools import *

def openPkl(filename):

    try:
        with open(filename, 'rb') as f:
            pkl = pickle.load(f)
        return pkl
    except EOFError as e:
        errorMessage('Unable to read file', 2, info='file: {:}'.format(filename[0]), detail='{:}'.format(e))
        return None

def saveSetup(obj):

    setup = Config()

    setup.fireball_name = obj.fireball_name_edits.text()
    setup.get_data = tryBool(obj.get_data_edits.currentText())
    setup.run_mode = obj.run_mode_edits.currentText()

    setup.arrival_times_file = obj.arrival_times_edits.text()
    setup.station_picks_file = obj.station_picks_edits.text()
    setup.replot_points_file = obj.points_name_edits.text()

    setup.lat_centre = tryFloat(obj.lat_centre_edits.text())
    setup.lon_centre = tryFloat(obj.lon_centre_edits.text())
    setup.deg_radius = tryFloat(obj.deg_radius_edits.text())
    setup.fireball_datetime = obj.fireball_datetime_edits.dateTime().toPyDateTime()
    obj.fatm_datetime_edits.setDateTime(setup.fireball_datetime)

    setup.t0 = tryFloat(obj.t0_edits.text())
    setup.v = tryFloat(obj.v_edits.text())
    setup.azimuth = tryAngle(obj.azim_edits.text())
    setup.zenith = tryAngle(obj.zangle_edits.text())
    setup.lat_i = tryFloat(obj.lat_i_edits.text())
    setup.lon_i = tryFloat(obj.lon_i_edits.text())
    setup.elev_i = tryFloat(obj.elev_i_edits.text())
    setup.lat_f = tryFloat(obj.lat_f_edits.text())
    setup.lon_f = tryFloat(obj.lon_f_edits.text())
    setup.elev_f = tryFloat(obj.elev_f_edits.text())

    setup.pos_i = tryPosition(setup.lat_i, setup.lon_i, setup.elev_i)
    setup.pos_f = tryPosition(setup.lat_f, setup.lon_f, setup.elev_f)

    setup.v_f = tryFloat(obj.vf_edits.text())

    # try:
    #     setup.trajectory = tryTrajectory(setup.t0, setup.v, setup.azimuth, setup.zenith, setup.pos_i, setup.pos_f, v_f=setup.v_f)
    # except AttributeError as e:
    #     errorMessage("Cannot build Trajectory!", 2)

    # setup.fragmentation_point = trySupracenter(fromTable(obj.fragmentation_point))

    frag_lat  = tryFloat(obj.lat_frag_edits.text())
    frag_lon  = tryFloat(obj.lon_frag_edits.text())
    frag_elev = tryFloat(obj.elev_frag_edits.text())
    frag_time = tryFloat(obj.time_frag_edits.text())
    setup.manual_fragmentation_search = trySupracenter(tryPosition(frag_lat, frag_lon, frag_elev), frag_time)

    setup.v_fixed = tryFloat(obj.v_fixed_edits.text())
    setup.enable_restricted_time = obj.restricted_time_check.isChecked()
    setup.restricted_time = obj.restricted_time_edits.dateTime().toPyDateTime()

    setup.azimuth_min = tryAngle(obj.azimuth_min_edits.text())
    setup.azimuth_max = tryAngle(obj.azimuth_max_edits.text())
    setup.zenith_min = tryAngle(obj.zangle_min_edits.text())
    setup.zenith_max = tryAngle(obj.zangle_max_edits.text())
    setup.lat_min = tryFloat(obj.lat_min_edits.text())
    setup.lat_max = tryFloat(obj.lat_max_edits.text())
    setup.lon_min = tryFloat(obj.lon_min_edits.text())
    setup.lon_max = tryFloat(obj.lon_max_edits.text())
    setup.elev_min = tryFloat(obj.elev_min_edits.text())
    setup.elev_max = tryFloat(obj.elev_max_edits.text())
    setup.t_min = tryFloat(obj.t_min_edits.text())
    setup.t_max = tryFloat(obj.t_max_edits.text())
    setup.v_min = tryFloat(obj.v_min_edits.text())
    setup.v_max = tryFloat(obj.v_max_edits.text())
    setup.weight_distance_min = tryFloat(obj.weight_distance_min_edits.text())
    setup.weight_distance_max = tryFloat(obj.weight_distance_max_edits.text())

    setup.pos_min = tryPosition(setup.lat_min, setup.lon_min, setup.elev_min)
    setup.pos_max = tryPosition(setup.lat_max, setup.lon_max, setup.elev_max)

    setup.supra_min = trySupracenter(setup.pos_min, setup.t_min)
    setup.supra_max = trySupracenter(setup.pos_max, setup.t_max)

    try:
        setup.traj_min = tryTrajectory(setup.t_min, setup.v_min, setup.azimuth_min, \
                                         setup.zenith_min, setup.pos_min, setup.pos_min)
        setup.traj_max = tryTrajectory(setup.t_max, setup.v_max, setup.azimuth_max, \
                                         setup.zenith_max, setup.pos_max, setup.pos_max)
    except AttributeError as e:
        errorMessage("Unable to build trajectory with given data", 0, info="Ignore if not using minimum or maximum trajectories", detail='{:}'.format(e))

    setup.observe_frag_no = tryInt(obj.frag_no_edits.text())

    setup.high_f = obj.high_f_edits.text()
    setup.high_b = obj.high_b_edits.text()
    setup.rm_stat = obj.rm_stat_edits.text()

    setup.stations = []
    temp_stats = fromTable(obj.extra_point)
    for stn in temp_stats:
        try:
            setup.stations.append(Station(stn[0], stn[1], Position(stn[2], stn[3], stn[4]), stn[5], stn[7]))
        except TypeError:
            print("WARNING: Station: {:}-{:} could not be read".format(stn[0], stn[1]))

    return setup

def loadDisplay(setup, obj):

    obj.fireball_name_edits.setText(setup.fireball_name)
    comboSet(obj.get_data_edits, setup.get_data)
    comboSet(obj.run_mode_edits, setup.run_mode)

    obj.arrival_times_edits.setText(setup.arrival_times_file)
    obj.station_picks_edits.setText(setup.station_picks_file)
    obj.points_name_edits.setText(setup.replot_points_file)

    obj.lat_centre_edits.setText(str(setup.lat_centre))
    obj.lon_centre_edits.setText(str(setup.lon_centre))
    obj.fatm_start_lat.setText(str(setup.lat_centre))
    obj.fatm_start_lon.setText(str(setup.lon_centre))
    obj.fatm_end_lat.setText(str(setup.lat_centre))
    obj.fatm_end_lon.setText(str(setup.lon_centre))
    obj.deg_radius_edits.setText(str(setup.deg_radius))
    obj.fireball_datetime_edits.setDateTime(setup.fireball_datetime)
    obj.fatm_datetime_edits.setDateTime(setup.fireball_datetime)

    # obj.t0_edits.setText(str(setup.trajectory.t))
    # obj.v_edits.setText(str(setup.trajectory.v))
    # obj.azim_edits.setText(str(setup.trajectory.azimuth.deg))
    # obj.zangle_edits.setText(str(setup.trajectory.zenith.deg))
    # obj.lat_i_edits.setText(str(setup.trajectory.pos_i.lat))
    # obj.lon_i_edits.setText(str(setup.trajectory.pos_i.lon))
    # obj.elev_i_edits.setText(str(setup.trajectory.pos_i.elev))
    # obj.lat_f_edits.setText(str(setup.trajectory.pos_f.lat))
    # obj.lon_f_edits.setText(str(setup.trajectory.pos_f.lon))
    # obj.elev_f_edits.setText(str(setup.trajectory.pos_f.elev))

    obj.vf_edits.setText(str(setup.v_f))

    # frag_list = []
    # for element in setup.fragmentation_point:
    #     frag_list.append(element.toList())

    # toTable(obj.fragmentation_point, frag_list)

    if setup.manual_fragmentation_search == None:
        obj.lat_frag_edits.setText('')
        obj.lon_frag_edits.setText('')
        obj.elev_frag_edits.setText('')
        obj.time_frag_edits.setText('')
    else: 
        obj.lat_frag_edits.setText(str(setup.manual_fragmentation_search[0].position.lat))
        obj.lon_frag_edits.setText(str(setup.manual_fragmentation_search[0].position.lon))
        obj.elev_frag_edits.setText(str(setup.manual_fragmentation_search[0].position.elev))
        obj.time_frag_edits.setText(str(setup.manual_fragmentation_search[0].time))

    obj.v_fixed_edits.setText(str(setup.v_fixed))
    obj.restricted_time_check.setChecked(setup.enable_restricted_time)
    obj.restricted_time_edits.setDateTime(setup.restricted_time)

    obj.azimuth_min_edits.setText(str(setup.traj_min.azimuth.deg))
    obj.azimuth_max_edits.setText(str(setup.traj_max.azimuth.deg))
    obj.zangle_min_edits.setText(str(setup.traj_min.zenith.deg))
    obj.zangle_max_edits.setText(str(setup.traj_max.zenith.deg))
    obj.lat_min_edits.setText(str(setup.traj_min.pos_f.lat))
    obj.lat_max_edits.setText(str(setup.traj_max.pos_f.lat))
    obj.lon_min_edits.setText(str(setup.traj_min.pos_f.lon))
    obj.lon_max_edits.setText(str(setup.traj_max.pos_f.lon))
    obj.elev_min_edits.setText(str(setup.traj_min.pos_f.elev))
    obj.elev_max_edits.setText(str(setup.traj_max.pos_f.elev))
    obj.t_min_edits.setText(str(setup.traj_min.t))
    obj.t_max_edits.setText(str(setup.traj_max.t))
    obj.v_min_edits.setText(str(setup.traj_min.v))
    obj.v_max_edits.setText(str(setup.traj_max.v))
    obj.weight_distance_min_edits.setText(str(setup.weight_distance_min))
    obj.weight_distance_max_edits.setText(str(setup.weight_distance_max))

    obj.frag_no_edits.setText(str(setup.observe_frag_no))

    obj.high_f_edits.setText(str(setup.high_f))
    obj.high_b_edits.setText(str(setup.high_b))
    obj.rm_stat_edits.setText(str(setup.rm_stat))

    # toTableFromStn(obj.extra_point, setup.stations)

def loadSourcesIntoBam(bam):
    
    bam.setup.fragmentation_point = []
    bam.setup.trajectory = None

    bam.setup.frag_metadata = []
    bam.setup.traj_metadata = []

    if not hasattr(bam, 'source_list'):
        return bam

    for src in bam.source_list:

        if src.source_type == "Fragmentation":

            bam.setup.fragmentation_point.append(src.source)
            bam.setup.frag_metadata.append(src)

        elif src.source_type == "Ballistic":

            # Only one trajectory supported currently
            bam.setup.trajectory = src.source
            bam.setup.traj_metadata.append(src)

    return bam

def loadAtmos(bam, obj):
    
    avg_sp_sound = obj.prefs.avg_sp_sound
    
    if not hasattr(bam, 'atmos'):
        bam.atmos = Atmos(avg_sp_sound=avg_sp_sound)


def save(obj):
    
    #save setup
    obj.bam.setup = saveSetup(obj)

    if obj.bam.file_name is None:
        obj.bam.file_name = saveFile('.bam')

    with open(obj.bam.file_name, 'wb') as f:
        pickle.dump(obj.bam, f)

    print('setup saved')


def load(obj):

    # Open file dialog for user to select bam file
    filename = fileSearch(['bam file (*.bam)'], None)

    # bam file to obj
    bam = openPkl(filename)

    # save filename for autosaving feature
    bam.file_name = filename
    

    loadDisplay(bam.setup, obj)
    bam = loadSourcesIntoBam(bam)
    # loadAtmos(bam, obj)
    # print(bam.stats)
    obj.bam = bam
    