
import numpy as np

from Filters import butterworthBandpassFilter

def procTrace(trace, ref_datetime=None, resp=None, bandpass=[2, 8]):
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

    if resp is not None:
        trace.remove_response(inventory=resp, output="DISP")
        trace.remove_sensitivity(resp) 

    ampl_data = trace.data

    ampl_data = ampl_data[:len(time_data)]
    time_data = time_data[:len(ampl_data)] + offset

    # Init the butterworth bandpass filter
    butter_b, butter_a = butterworthBandpassFilter(low, high, 1.0/delta, order=2)

    # Filter the data
    ampl_data = scipy.signal.filtfilt(butter_b, butter_a, ampl_data)

    if backup:
        return ampl_data, time_data, raw_trace

    return ampl_data, time_data      
