
from supra.Stations.sourcefinder import sourceFinder

CHANNEL = "HHZ"
VERBOSE_MODE = True
SHOW_PLOT = False
LEGNTH_OF_SIGNAL = 5.00 # maximum length of a signal [s]
BANDPASS = {"freqmin": 2, "freqmax":8, "corners":4, "zerophase":False}
WINDOW_SIZE = 15 # length of sliding window [s]
WINDOW_OVERLAP = 1.0 # frac
STD_RANGE = 10.0 # num of std away from mean noise to call it a pick

def calcAllSigs(bam, prefs):

    stn_list = bam.stn_list

    for stn in stn_list:

        results = sourceFinder(stn.stream, channel=CHANNEL, verbose=VERBOSE_MODE, plot=SHOW_PLOT, \
            length_of_signal=LEGNTH_OF_SIGNAL, bandpass=BANDPASS, \
            window_size=WINDOW_SIZE, window_overlap=WINDOW_OVERLAP, std_range=STD_RANGE)

        stn.signals = results

    return bam.stn_list