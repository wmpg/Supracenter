""" Inputs a single-component waveform stream and 
    outputs the times of the signals
    Author: Luke McFadden
"""


def usage():
    ''' Function to print out the usage of sourcefinder.py
    '''

    print("""
            ######### usage ###########
             python -m sourcefinder.py
             arguments:
             file_name - [string]
             keyword arguments:
             channel - [string or None]
             verbose - [Boolean]
            ###########################
        """)

import numpy as np
import matplotlib.pyplot as plt
import more_itertools as mit
from scipy.stats import median_abs_deviation

from supra.Stations.Pick import Pick
from supra.Stations.obspyhelper import loadStreamFromFile, cleanWaveform
from supra.Stations.listhelper import rmBelow, cluster, degroupPicks


def statRemove(waveform_data, window_size=0, window_overlap=1.0, std_range=5):

    """ Finds the points in a list which are significantly outside the noise standard distribution.
    Uses a sliding window (which may overlap itself) to compare points locally to each other rather than
    over the entire dataset

    Usage:
    statRemove(data_list, window_size=100, window_overlap=0.1, std_range=10)
    >>> returns a list of points in data_list which are approximately 10 standard deviations away
    from the noise

    This function uses the median and median absolute deviation in place of the mean and standard deviation
    in order to better approximate the distribution of noise in the data

    Arguments:
    waveform_data [list] - a list of the data to find significant points from

    Keyword Arguments:
    window_size [int] - the length of the sliding window [in number of points]
    window_overlap [float] - amount of the window length to slide over at each step 
                            (0.01 - nearly complete overlap, 1.0 - no overlap, don't use 0)
    std_range [float] - number of standard deviations away from the mean to be considered significant

    Returns:
    An ordered list of the indicies from waveform_data which were found to be significantly different from
    the noise

    """

    # Number of points in trace
    stream_len = len(waveform_data)

    # A window_overlap of 0 will result in the window not moving, make it 0.01
    if window_overlap == 0.0:
        print("[WARNING] 0.0 used for window overlap, turning into 0.01")
        window_overlap = 0.01

    # If window_size is not given, make it the whole stream
    if window_size == 0:
        window_size = stream_len

    # Number of points in a window (as close as we can get)
    pts_per_window = int(window_size)

    # number of points to shift window each step
    amount_to_shift_window = int(pts_per_window*window_overlap)

    # List of indicies of waveform_data where the windows begin and end
    window_starts = np.arange(0, stream_len, amount_to_shift_window)
    window_ends = window_starts + pts_per_window

    timing_indicies = []

    # Sliding window
    for i in range(len(window_starts)):

        # Create window
        window_data = waveform_data[window_starts[i]:window_ends[i]]

        # Use median here because it's a better representation of the distribution of noise
        window_mean = np.median(waveform_data)
        window_std =  median_abs_deviation(waveform_data, scale='normal')

        cutoff = window_mean + std_range*window_std
        
        # Only select points above the cutoff 
        window_indicies = rmBelow(window_data, cutoff)
        
        window_indicies = list(np.array(window_indicies) + window_starts[i])
        timing_indicies += window_indicies

    ### Combine and order windows
    # remove duplicates which may come from window overlaps
    timing_indicies = list(set(timing_indicies))
    timing_indicies.sort()

    return timing_indicies



def groups2PickObjs(groups, conversion_factor):

    '''
    Converts the groups of indicies into their pick objects with a conversion factor (index -> time)

    Usage:
    groups2PickObjs(mit.consecutive_groups(my_list), 0.25 #seconds per point#)

    Arguments:
    groups [list] - a list of clusters of points to be put into pick objects, each cluster as a list of indicies
    conversion_factor [float] - number of seconds per index

    returns [list of Pick Obj] - a list of converted Pick Objects
    '''

    pick_list = []
    
    # Take consecutive points and group them into picks
    for grp in groups:
        
        grp = list(grp)
        
        # Get boundaries of groups
        first_pt, last_pt = grp[0], grp[-1]

        # convert from point order to a time
        first_pt *= conversion_factor
        last_pt *= conversion_factor

        pk = Pick(first_pt, last_pt)

        pick_list.append(pk)

    return pick_list



def sourceFinder(stream, channel=None, verbose=False, plot=False, length_of_signal=5, \
        bandpass={"freqmin": 2, "freqmax":8, "corners":4, "zerophase":False}, \
        window_size=15, window_overlap=1.0, std_range=10):
    
    '''
    Takes a trace from a given stream and finds possible signals, returns the times relative to the start of the stream where
    the picks are.

    Usage:
    sourceFinder('myFile.mseed', channel='HHZ', **kwargs)
    >>> a list of times where the picks are in the HHZ channel of the provided file

    How it works:
    The waveform is assumed to be made of two components: noise and signal. The noise is made of random values
    under a normal distribution. The signal is assumed not to follow the noise standard distribution. To estimate
    this noise standard distribution, the median and the medain absolute deviation are taken instead of the mean
    and the standard deviation to try and separate the noise from the signals. Any data above a certain number of
    standard deviations away from this standard distribution are said to be signal.

    The method done by looking at local sections of the waveform, and finding a signal has a large amplitude in
    relation to the data immediately around it, not the whole waveform. This is done by using a sliding window.
    A sliding window of a certain length is moved along the data, possibly overlaping itself depending on
    window_overlap, and the statistical method described above is done on each window. The signals are then
    returned in a human-readable form.

    Arguments:
    stream [Obspy Stream Obj] - stream

    Keyword Arguments:
    channel [String] - Channel of the trace to do analysis on, defaults to the first one in the Stream
    verbose [Boolean] - Switch to print program status
    plot [Boolean] - Switch to plot the waveform with the signals highlighted
    length_of_signal [float] - signals within this many seconds of each other will be combined into one signal
                                regardless of if a signal is detected between them 
                                (ideally on the order of a couple seconds, source-type dependant)
    bandpass [Dictionary] - keyword arguments to pass to the Obspy bandpass filter
                            Format: {"freqmin": 2, "freqmax":8, "corners":4, "zerophase":False}

                            Obspy's docs are shown here:
                            :param freqmin: Pass band low corner frequency.
                            :param freqmax: Pass band high corner frequency.
                            :param corners: Filter corners / order.
                            :param zerophase: If True, apply filter once forwards and once backwards.
                                This results in twice the filter order but zero phase shift in
                                the resulting filtered trace.
                            (defaults work fine here, they just provide a quick filtering, and are not 
                            meant to be optimal)
    window_size [float] - length of sliding windows in seconds which scans across the waveform
                            (A few seconds is optimal, if it is too small, longer signals will not be detected)
    window_overlap [float] - how much of the window length to slide over at each step (0.0, 1.0]
    std_range [float] - number of standard deviations away from the noise to consider a signal 
                        (lower -> more sensitive, higher -> less sensitive, ideally keep above 3-5)

    Returns:
    [list] - A list of the signals found, with each signal time in the form: [first_point, last_point].
             The signal time is given in seconds after the beginning of the trace selected.
    '''

    ###############
    # LOAD STREAM #
    ###############

    # if verbose:
    #     print("[INFO] Loading stream from file ({:})".format(file_name))

    # stream = loadStreamFromFile(file_name)
    
    # No traces in stream (empty file), or couldn't find file
    if len(stream) == 0:
        usage()
        return None

    ####################################
    # EXTRACT SINGLE TRACE FROM STREAM #
    ####################################

    if verbose:
        print("[INFO] Loading trace from stream")

    try:
        single_channel_stream = stream.select(channel=channel)
    
    except AttributeError as e:
        print("[ERROR] Channel must be of type string")
        print("... {:}".format(e))
        usage()

        single_channel_stream = []

    # If exception was raised, or if channel wasn't found
    if len(single_channel_stream) > 0:
        trace = single_channel_stream[0]

    else:
        print("[WARNING] Unable to find channel '{:}', using first trace".format(channel))
        trace = stream[0]


    if verbose:
        print("[INFO] Channel {:} selected".format(channel))

    ######################
    # WAVEFORM FILTERING #
    ######################

    if verbose:
        print("[INFO] cleaning waveform (bandpass and detrend)")

    # Clean waveform as shown by Obspy documentation
    cleanWaveform(trace, bandpass)

    ######################
    # EXTRACT TRACE DATA #
    ######################

    # Frequency
    freq = trace.stats.sampling_rate
    
    # Number of points
    npts = trace.stats.npts

    # Time between points
    delta = trace.stats.delta

    # Amplitude of waveform data
    waveform_data_raw = trace.data
    
    # Times of each point
    waveform_time = np.arange(0, npts/freq, delta)

    waveform_data = np.abs(waveform_data_raw)

    if verbose:
        print("[INFO] extracted data")

    ################
    # FIND SOURCES #
    ################
    
    # Take out significant points from noise - get indicies
    timing_indicies = statRemove(waveform_data, window_size=int(window_size*freq), \
                            window_overlap=window_overlap, std_range=std_range)

    if verbose:
        print("[INFO] selected indicies of picks")

    #########################
    # GENERATE PICK OBJECTS #
    #########################

    # group consecutive indicies together
    timing_groups = mit.consecutive_groups(timing_indicies)

    # Turn groups into pick objects
    pick_list = groups2PickObjs(timing_groups, delta)

    if verbose:
        print("[INFO] making pick objects")

    ######################
    # COMBINE LIKE PICKS #
    ######################

    # For points within a default length beside each other, combine them
    filtered_list = degroupPicks(cluster(pick_list, length_of_signal))
    
    if filtered_list is None:
        return None

    """Put as list for single function output
    Bolide Acoustic Modeling uses the pick objects, so the return would be here when added to that
    """

    # Turn from pick objects into human-readable form
    output_list = []
    for p in filtered_list:
        output_list.append([p.first_pt, p.last_pt])


    ########
    # PLOT #
    ########

    if plot:
        plt.plot(waveform_time, waveform_data_raw, c='k')
        
        for r in output_list:
            plt.axvspan(r[0], r[1], alpha=0.5, color='red')

        plt.show()


    return output_list
    

if __name__ == "__main__":
    pass

    ################################################################
    # USER CONFIG
    ################################################################

    ### Arguments

    CHANNEL = "HHZ"
    VERBOSE_MODE = True
    SHOW_PLOT = True
    LEGNTH_OF_SIGNAL = 5.00 # maximum length of a signal [s]
    BANDPASS = {"freqmin": 2, "freqmax":8, "corners":4, "zerophase":False}
    WINDOW_SIZE = 15 # length of sliding window [s]
    WINDOW_OVERLAP = 1.0 # frac
    STD_RANGE = 10.0 # num of std away from mean noise to call it a pick


    ##########################################################################################
    # Function call
    ##########################################################################################

    results = sourceFinder(file_name, channel=CHANNEL, verbose=VERBOSE_MODE, plot=SHOW_PLOT, \
        length_of_signal=LEGNTH_OF_SIGNAL, bandpass=BANDPASS, \
        window_size=WINDOW_SIZE, window_overlap=WINDOW_OVERLAP, std_range=STD_RANGE)
    print(results)

 