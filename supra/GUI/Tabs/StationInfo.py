
import os

from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *

from supra.GUI.Tools.CustomWidgets import StationEx
from supra.GUI.Dialogs.StationSource import StationWindow
from supra.GUI.Tools.GUITools import *
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
        stat_obj.response = stn.response
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

def outputStation(obj):

    print("Station List:\n")
    for stn in obj.bam.stn_list:
        meta = stn.metadata
        print("['{:}-{:}', {:}, {:}, {:}],".format(meta.network, meta.code, \
                    meta.position.lat, meta.position.lon, meta.position.elev))

def countStation(obj):
    
    print("Station Count:\n")
    print("##########################")
    
    seis_count = 0
    infra_count = 0

    for stn in obj.bam.stn_list:
        stream = stn.stream

        if len(stn.stream.select(channel="*HZ")) > 0:
            seis_count += 1

        if len(stn.stream.select(channel="*DF")) > 0:
            infra_count += 1

    print("Seismic Stations:    {:}".format(seis_count))
    print("Infrasound Stations: {:}".format(infra_count))

def getStations(obj):

    # Create fireball folder
    if not obj.checkForWorkDir():
        return None

    if not hasattr(obj.bam.setup, "fireball_name"):
        errorMessage("The .bam file has not been created, please save your event file before downloading stations", 1)
        return None

    # #Build seismic data path
    obj.dir_path = os.path.join(obj.prefs.workdir, obj.bam.setup.fireball_name)


    if obj.bam.setup.get_data:
        
        errorMessage("Stations are being downloaded", 0, info="Check terminal for download status", title="Station Download")

        #Download all waveform files which are within the given geographical and temporal range ###
        stn_list = getAllStations(obj.bam.setup.lat_centre, obj.bam.setup.lon_centre, \
                obj.bam.setup.deg_radius, obj.bam.setup.fireball_datetime, \
                obj.dir_path, obj=obj)

        addStationWidgets(obj, stn_list)

    else:
        print('Get data is turned off, so that data is not accidentally downloaded. If this was intentional, please turn it back on')


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
            stn = Station(meta, stn_ex.stream, response=stn_ex.response)

            stn_list.append(stn)

    obj.bam.stn_list = stn_list
    save(obj)
