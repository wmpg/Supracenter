
import numpy as np

from supra.Stations.Filters import butterworthBandpassFilter
import scipy.signal
from scipy.fft import fft

def procTrace(trace, ref_datetime=None, resp=None, bandpass=[2, 8], backup=False):
    ''' procTrace filters a waveform given a response file and a reference time

    Arguments:
    trace - A single obspy trace to filter

    Keyword Arguments:
    ref_datetime - datetime object to reference to
    resp - response data to remove from the trace
    bandpass - [low, high] bandpass filters to use for the data
    

    ref_datetime = self.bam.setup.fireball_datetime
    '''

    raw_trace = trace.copy()


    # Obtain metadata from trace
    delta           = trace.stats.delta
    start_datetime  = trace.stats.starttime.datetime
    end_datetime    = trace.stats.endtime.datetime
    npts            = trace.stats.npts
    sampling_rate   = trace.stats.sampling_rate

    if ref_datetime is not None:
        offset = (start_datetime - ref_datetime).total_seconds()
    else:
        offset = 0

    # self.current_waveform_delta = delta
    time_data = np.arange(0,  npts/sampling_rate, delta)

    trace.detrend()

    # Remove sensitivity and remove response do the same thing - response is better
    if resp is not None:
        trace.remove_response(inventory=resp, output="DISP")
        # trace.remove_sensitivity(resp) 

    ampl_data = trace.data

    ampl_data = ampl_data[:len(time_data)]
    time_data = time_data[:len(ampl_data)] + offset

    # Init the butterworth bandpass filter
    if bandpass is not None:

        low =  bandpass[0]
        high = bandpass[1]

        butter_b, butter_a = butterworthBandpassFilter(low, high, 1.0/delta, order=2)

        # Filter the data
        ampl_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(ampl_data))


    if backup:
        return ampl_data, time_data, raw_trace

    return ampl_data, time_data      

def subTrace(trace, begin_time, end_time, ref_time, clean=None):
    """ Returns the trace between beginTime and endTime given as two times after reference"""
    

    delta = trace.stats.delta
    start_datetime = trace.stats.starttime.datetime
    end_datetime = trace.stats.endtime.datetime

    offset = (start_datetime - ref_time).total_seconds()

    if clean is not None:
        waveform_data, time_data = procTrace(trace, **clean)
    else:
        waveform_data = trace.data
        time_data = np.arange(0, trace.stats.npts / trace.stats.sampling_rate, \
                     delta)

        waveform_data = waveform_data[:len(time_data)]
        time_data = time_data[:len(waveform_data)] + offset

    number_of_pts_per_s = trace.stats.sampling_rate

    len_of_region = end_time - begin_time

    num_of_pts_in_roi = len_of_region*number_of_pts_per_s

    num_of_pts_in_offset = np.abs(number_of_pts_per_s*offset)

    num_of_pts_to_roi = begin_time*number_of_pts_per_s

    pt_0 = int(num_of_pts_in_offset + num_of_pts_to_roi)
    pt_1 = int(pt_0 + num_of_pts_in_roi)

    cut_waveform = waveform_data[pt_0:pt_1]
    cut_time = time_data[pt_0:pt_1]

    return cut_waveform, cut_time

def genFFT(waveform, time_data):

    sampling_rate = 1/(time_data[1] - time_data[0])
    sps = sampling_rate
    dt = 1/sampling_rate
    length = len(waveform)
    len_of_region = time_data[-1] - time_data[0]

    freq = np.linspace(1/len_of_region, (sps/2), length)*sps/length

    FAS = abs(fft(waveform))

    return freq, FAS