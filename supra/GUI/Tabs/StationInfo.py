
import os

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

from supra.GUI.Tools.CustomWidgets import StationEx
from supra.GUI.Dialogs.StationSource import StationWindow
from supra.Utils.Classes import Position
from supra.Files.SaveLoad import save
from supra.Stations.FetchStations import getAllStations
from supra.Stations.StationObj import Metadata, Station

def addStationWidgets(obj, stn_list):
    obj.stat_widget_lst = []
    for ii, stn in enumerate(stn_list):

        stat_obj = StationEx()
        stat_obj.network.setText('{:}'.format(stn.metadata.network))
        stat_obj.code.setText('{:}'.format(stn.metadata.code))
        stat_obj.position.lat.setText('{:}'.format(stn.metadata.position.lat))
        stat_obj.position.lon.setText('{:}'.format(stn.metadata.position.lon))
        stat_obj.position.elev.setText('{:}'.format(stn.metadata.position.elev))
        stat_obj.name.setText('{:}'.format(stn.metadata.name))
        stat_obj.stream = stn.stream
        stat_obj.toggle.setState(stn.metadata.enabled)

        obj.stat_widget_lst.append(stat_obj) 

    drawStationWidgets(obj)

def clearStationWidgets(obj):
    for i in reversed(range(obj.station_table_layout.count())): 
        obj.station_table_layout.itemAt(i).widget().setParent(None)
        obj.stat_widget_lst = []

def drawStationWidgets(obj):

    for widget in obj.stat_widget_lst:
        obj.station_table_layout.addWidget(widget)

def addStation(obj):
    obj.p = StationWindow(obj.bam)
    obj.p.setGeometry(QRect(500, 400, 500, 400))
    obj.p.show()

def refreshStation(obj):
    clearStationWidgets(obj)
    addStationWidgets(obj, obj.bam.stn_list)

def getStations(obj):

    # Create fireball folder
    if not obj.checkForWorkDir():
        return None

    # #Build seismic data path
    obj.dir_path = os.path.join(obj.prefs.workdir, obj.bam.setup.fireball_name)


    if obj.bam.setup.get_data:
        
        #Download all waveform files which are within the given geographical and temporal range ###
        stn_list = getAllStations(obj.bam.setup.lat_centre, obj.bam.setup.lon_centre, \
                obj.bam.setup.deg_radius, obj.bam.setup.fireball_datetime, \
                obj.dir_path, obj=obj)

    else:
        print("WARNING: get_data is turned off, this is the data only on this machine!")
        
    # data_file_path = os.path.join(obj.dir_path, DATA_FILE)

    # if os.path.isfile(data_file_path):
        
    #     stn_list = readStationAndWaveformsListFile(data_file_path)

    # else:
    #     # print('Station and waveform data file not found! Download the waveform files first!')
    #     errorMessage("Station and waveform data file not found!", 1, detail='Download the waveform files first, make sure get data is on')
    #     return None

    addStationWidgets(obj, stn_list)
    # display all stations
    # but don't save

def loadStations(obj):

    clearStationWidgets(obj)
    addStationWidgets(obj, obj.bam.stn_list)
    drawStationWidgets(obj)

def saveStations(obj):

    stn_list = []

    for stn_ex in obj.stat_widget_lst:
        if stn_ex.toggle.getState():
            pos = Position(float(stn_ex.position.lat.text()), float(stn_ex.position.lon.text()), \
                            float(stn_ex.position.elev.text()))

            meta = Metadata(stn_ex.network.text(), stn_ex.code.text(), pos, stn_ex.name.text())
            meta.enabled = True
            stn = Station(meta, stn_ex.stream)

            stn_list.append(stn)
        

    obj.bam.stn_list = stn_list
    save(obj)