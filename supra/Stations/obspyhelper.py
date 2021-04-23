'''
Misc obspy helper functions taken from the examples provided by Obspy 
'''

import obspy


def loadStreamFromFile(file_name):
    '''
    Extracts stream object from mseed file

    Arguments:
    file_name [String] - file path to the mseed file

    Returns:
    stream [Stream Obj] - Obspy stream object of the mseed file
    '''

    try:
        
        stream = obspy.read(file_name)
    
    except FileNotFoundError as e:
        
        print("[ERROR] {:}".format(e))
        
        # use [] instead of None so that len(stream) can be compared
        return []
    
    return stream

def cleanWaveform(tr, kwargs={"freqmin": 2, "freqmax":8, "corners":4, "zerophase":False}):
    '''
    clean and filter waveform using obspy libraries

    Arguments:
    tr [Trace Obj] - trace to be processed

    Keyword Arguments:
    kwargs [Dictionary] - dictionary of obspy filter keyword arguments
    '''

    tr.detrend()
    tr.filter("bandpass", **kwargs)


''' 
testing not really needed for these because they are already done in the library functions of Obspy
The first function is just an Obspy function with a try-except statement around it
The second function is an example given by Obspy

If neither of these work, it is on Obspy
'''