import obspy
import numpy as np

from supra.Stations.Filters import butterworthBandpassFilter
import scipy.signal
from scipy.fft import fft
from obspy.signal.filter import bandpass as bandpassFunc

from scipy.signal import butter, lfilter

from scipy.signal import butter, sosfilt, sosfreqz

from scipy import signal

def butter_bandpass(lowcut, highcut, fs, order=6):
    """ Scipy suggested bandpass filter with fix found on stackoverflow for low frequencies

    Arguments:
    lowcut [float] - Lowcut frequency [Hz]
    highcut [float] - Highcut frequency [Hz]
    fs [float] - Sampling rate [s]
    order [integer] - Order to run the bandpass (6 is good, 9 is too much, 
                        2 is okay for very low frequencies where higher orders don't work)
    """

    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    sos = butter_bandpass(lowcut, highcut, fs, order=order)
    y = sosfilt(sos, data)
    return y


def bandpassFilter(bandpass, delta, ampl_data):
    low =  bandpass[0]
    high = bandpass[1]

    butter_b, butter_a = butterworthBandpassFilter(low, high, 1.0/delta, order=2)

    # Filter the data
    ampl_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(ampl_data))

    return ampl_data


def procTrace(trace, ref_datetime=None, resp=None, bandpass=[2, 8], backup=False):
    ''' procTrace filters a waveform given a response file and a reference time

    Arguments:
    trace - A single obspy trace to filter

    Keyword Arguments:
    ref_datetime - datetime object to reference to
    resp - response data to remove from the trace
    bandpass - [low, high] bandpass filters to use for the data
    backup [Boolean] - keep a raw trace
    '''


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

    # Split traces into groups because some traces may be missing a few seconds of data - these are marked by red lines
    partial_traces = trace.split()

    total_ampl = []
    total_time = []

    for tr in partial_traces:
        tr.detrend()

        # Remove sensitivity and remove response do the same thing - response is better (obspy)

        if resp is not None:
            tr.remove_response(inventory=resp, output="DISP")
            # trace.remove_sensitivity(resp) 

        ampl_data = tr.data

        time_data = tr.times(reftime=obspy.core.utcdatetime.UTCDateTime(ref_datetime))
        ampl_data = ampl_data[:len(time_data)]
        # time_data = time_data[:len(ampl_data)]

        raw_trace = ampl_data
        # Init the butterworth bandpass filter
        if bandpass is not None:

            # ampl_data = bandpassFilter(ampl_data, bandpass[0], bandpass[1], sampling_rate, order=5)
            ampl_data = butter_bandpass_filter(ampl_data, bandpass[0], bandpass[1], sampling_rate, order=6)

        total_ampl.append(ampl_data)
        total_time.append(time_data)

    if backup:
        return total_ampl, total_time, raw_trace

    return total_ampl, total_time   
       

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
        time_data = trace.times(reftime=obspy.core.utcdatetime.UTCDateTime(ref_time))

        waveform_data = waveform_data[:len(time_data)]
        time_data = time_data[:len(waveform_data)]

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

# def genFFT(waveform, time_data):

#     sampling_rate = 1/(time_data[1] - time_data[0])
#     sps = sampling_rate
#     dt = 1/sampling_rate
#     length = len(waveform)
#     len_of_region = time_data[-1] - time_data[0]

#     freq = np.linspace(1/len_of_region, (sps/2), length)*sps/length

#     FAS = abs(fft(waveform))

#     return freq, FAS

def genFFT(waveform, sampling_rate):


    freqs, psd = signal.welch(waveform)

    return freqs*sampling_rate, psd

def genSHM(freq, sampling_rate, time):

    # npts = time*sampling_rate

    delta = 1/sampling_rate

    time_data = np.arange(0, time, delta)

    ampl_data = np.sin(2*np.pi*freq*time_data)

    return np.array(ampl_data), np.array(time_data)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    sampling_rate = 60 #Hz
    time = 60 #seconds

    # generate a waveform as a superposition of other waves
    # a, t = genSHM(1, sampling_rate, time)
    a = np.zeros(sampling_rate*time)

    for f in [2, 4, 0.2, 0.5, 1]:
        temp_a, t = genSHM(f, sampling_rate, time)
        a = a + temp_a

    # plt.subplot(1, 2, 1)
    # plt.plot(t, a, label="Raw")

    # bandpass = [1/35, 1/25]

    # a_new = bandpassFilter(bandpass, 1/sampling_rate, a)

    # plt.plot(t, a_new, label="Scipy")

    # a_obs = butter_bandpass_filter(a, bandpass[0], bandpass[1], sampling_rate)

    # plt.plot(t, a_obs, label="Scipy Cookbook")

    # plt.legend()

    # plt.subplot(1, 2, 2)


    # freq, FAS = genFFT(a_new, t)
    # plt.semilogx(freq, FAS, label="Scipy")



    # plt.axvline(x=bandpass[0]/60, c='k')
    # plt.axvline(x=bandpass[1]/60, c='k')

    # plt.legend()
    # plt.show()


    # Sample rate and desired cutoff frequencies (in Hz).
    fs = sampling_rate
    lowcut = 0.3
    highcut = 0.6

    # Plot the frequency response for a few different orders.
    # plt.figure(1)
    # plt.clf()
    # for order in [3, 6, 9]:
    #     sos = butter_bandpass(lowcut, highcut, fs, order=order)
    #     w, h = sosfreqz(sos, worN=2000)
    #     plt.semilogx((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    # plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)], '--', label='sqrt(0.5)')
    # plt.xlabel('Frequency (Hz)')
    # plt.ylabel('Gain')
    # plt.grid(True)
    # plt.legend(loc='best')

    # Filter a noisy signal.
    plt.figure(1)
    plt.clf()
    plt.plot(t, a, label='Noisy signal')

    y = butter_bandpass_filter(a, lowcut, highcut, fs, order=6)
    plt.plot(t, y, label='Filtered signal')



    plt.xlabel('time (seconds)')
    plt.grid(True)
    plt.axis('tight')
    plt.legend(loc='upper left')

    plt.show() 


    # freq, FAS = genFFT(a, t)
    freq, FAS = genFFT(a, sampling_rate)
    plt.semilogx(freq, FAS, label="Raw")

    freq, FAS = genFFT(y, sampling_rate)
    plt.semilogx(freq, FAS, label="Scipy Cookbook")


    plt.axvline(x=lowcut, c='k')
    plt.axvline(x=highcut, c='k')

    plt.legend()
    plt.show()