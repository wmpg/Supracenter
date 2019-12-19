import PyQt5
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQt5.QtWebEngineWidgets import QWebEnginePage as QWebPage
from PyQt5.QtWebEngineWidgets import QWebEngineView as QWebView
from PyQt5.QtNetwork import *
from OpenGL import GL
import sys
from optparse import OptionParser
import os

import folium
from folium import features



# tooltip = 'Click me!'
# folium.Marker([45.3288, -121.6625], 
#     popup='<i>Mt. Hood Meadows</i>', 
#     tooltip=tooltip, 
#     icon=folium.Icon(color='red', icon='info-sign')).add_to(m)

# folium.Circle( #meters
#     radius=100,
#     location=[45.5244, -122.6699],
#     popup='The Waterfront',
#     color='crimson',
#     fill=False,
# ).add_to(m)

# folium.CircleMarker( #pixels
#     location=[45.5215, -122.6261],
#     radius=100,
#     popup='Laurelhurst Park',
#     color='#3186cc',
#     fill=True,
#     fill_color='#3186cc'
# ).add_to(m)
def addTrajectory(setup, m):
    folium.CircleMarker(
        location=[setup.lat_centre, setup.lon_centre],
        radius=10,
        popup='Geometric Landing Point',
        color='#3186cc',
        fill=True,
        fill_color='#3186cc'
    ).add_to(m)

    try:
        folium.PolyLine(
            locations=((setup.trajectory.pos_i.lat, setup.trajectory.pos_i.lon),
             (setup.trajectory.pos_f.lat, setup.trajectory.pos_f.lon)),
             color='blue',
            weight=2.5).add_to(m)
    except:
        pass
        
    return m

def addStations(setup, stn_list, m):
    for stn in stn_list:

        if stn.channel == 'BDF':
            c = 'green'
        else:
            c = 'red'

        folium.Marker([stn.position.lat, stn.position.lon], 
        popup='{:}-{:} \n Latitude: {:} \n Longitude: {:} \n Elevation: {:}'.format(stn.network, stn.code, stn.position.lat, stn.position.lon, stn.position.elev), 
        icon=folium.Icon(color=c, icon='info-sign')).add_to(m)

    return m

def htmlBuilder(setup, stn_list):
    
    m = folium.Map(location=[setup.lat_centre, setup.lon_centre])
    m.add_child(folium.LatLngPopup())

    m = addStations(setup, stn_list, m)
    m = addTrajectory(setup, m)

    filename = os.path.join(setup.working_directory, 'index.html')
    m.save(filename)

    return filename

 


