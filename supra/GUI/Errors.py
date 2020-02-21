from pyqtgraph.Qt import QtGui, QtCore
from supra.GUI.GUITools import *

def errorCodes(obj, attr, debug=True):

    if not debug:
        return False

    if hasattr(obj, attr):
        return False

    else:
        if attr == 'current_station_waveform':
            errorMessage("Current waveform does not exist!", 2, info='Please load the waveforms', detail='[ATTR ERROR {:}]'.format(attr))
        if attr == 'stn_list':
            errorMessage('Station list has not been generated yet', 1, detail='''[ATTR ERROR {:}] Either the stations are not downloaded yet (check "Stations" tab), or the stations have not been loaded into the "Make Picks" tab (Load Station Data button)'''.format(attr))
        if attr == 'current_station':
            errorMessage('Current station not defined', 1, info='Could not update GUI!', detail='[ATTR ERROR {:}]'.format(attr))
        if attr == 'trajectory':
            errorMessage('Trajectory not defined', 1, detail='[ATTR ERROR {:}] Define a trajectory in the setup toolbar.'.format(attr))
        if attr == 'fragmentation_point':
            errorMessage('Fragmentation point not defined', 1, detail='[ATTR ERROR {:}] Define a fragmentation point in the setup toolbar.'.format(attr))
        if attr == 'sounding_file':
            errorMessage('No sounding file defined', 2, detail='[ATTR ERROR {:}] Define a sounding file'.format(attr))
        return True