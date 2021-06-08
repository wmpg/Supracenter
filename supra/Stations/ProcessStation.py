import obspy
import numpy as np

from supra.Stations.Filters import butterworthBandpassFilter
import scipy.signal
from scipy.fft import fft
from obspy.signal.filter import bandpass as bandpassFunc

from scipy.signal import butter, lfilter

from scipy.signal import butter, sosfiltfilt, sosfreqz
from scipy.interpolate import interp1d
from scipy import signal
from scipy.signal import hilbert, chirp

def butter_bandpass(lowcut, highcut, fs, order=6):
    """ Scipy suggested bandpass filter with fix found on stackoverflow for low frequencies

    Arguments:
    lowcut [float] - Lowcut frequency [Hz]
    highcut [float] - Highcut frequency [Hz]
    fs [float] - Sampling rate [s]
    order [integer] - Order to run the bandpass (6 is good, 9 is too much, 
                        2 is okay for very low frequencies where higher orders don't work)
    """

    ### TODO
    # Add forward and backward filter so that there are no artifacts at the beginning of a signal

    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    return sos

def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):

    """ sos adds another term, this works better for low frequency bandpasses
    """

    sos = butter_bandpass(lowcut, highcut, fs, order=order)
    y = sosfiltfilt(sos, data)
    return y


def bandpassFilter(bandpass, delta, ampl_data):

    """ Old filter method
    """

    low =  bandpass[0]
    high = bandpass[1]

    butter_b, butter_a = butterworthBandpassFilter(low, high, 1.0/delta, order=2)

    # Filter the data
    ampl_data = scipy.signal.filtfilt(butter_b, butter_a, np.copy(ampl_data))

    return ampl_data

def procStream(stn, ref_time=None):

    mseed = stn.stream.copy()
    resp = stn.response

    # Not sure why this error happens, float modulo Obspy error??
    try:
        mseed.merge()
    except ZeroDivisionError:
        pass
    except Exception:
        pass

    gaps = mseed.get_gaps()

    gap_times = []
    for gap in gaps:
        start = gap[4]
        end = gap[5]

        ref = obspy.core.utcdatetime.UTCDateTime(ref_time)

        start_ref = start - ref
        end_ref = end - ref

        gap_times.append([start_ref, end_ref])


    return mseed, resp, gap_times

def findChn(st, chn):
    try:
        st = st.select(channel=chn)[0]
    except IndexError:
        return []

    return st


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

        try:
            time_data = tr.times(reftime=obspy.core.utcdatetime.UTCDateTime(ref_datetime))
        except TypeError:
            time_data = np.arange(0,  npts/sampling_rate, delta)
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

def findDominantPeriod(wave, time):

    """ Estimate the dominant period using the zero-crossing method
        This should be redone with PSDs in the future, once the infrasound curve method actually works
    """

    # Number of crossings to check. Too many and you are no longer getting the dominant period at the time
    num_of_cross = 4

    # find zero-crossings
    zero_cross = []
    for ii in range(len(wave) - 1):
        if wave[ii+1]/wave[ii] < 0:
            zero_cross.append((time[ii+1] + time[ii])/2)

    if len(zero_cross) == 0:
        return np.nan

    if len(zero_cross) < num_of_cross:
        num_of_cross = len(zero_cross)

    estimates = []
    # Estimate using 1/2 wavelengths
    for i in range(num_of_cross - 1):
        estimates.append((zero_cross[i+1] - zero_cross[i])*2)

    # Estimate using full wavelengths
    for i in range(num_of_cross - 2):
        estimates.append(zero_cross[i+2] - zero_cross[i])


    return np.mean(estimates)

def findDominantPeriodPSD(wave, sf, normalize=False):
    
    freq, FAS = genFFT(wave, sf)
    ii = np.nanargmax(FAS)

    if normalize:
        FAS /= np.nanmax(FAS)


    dom_freq = freq[ii]

    if dom_freq == 0:
        return np.nan, freq, FAS

    dom_period = 1/dom_freq

    return dom_period, freq, FAS


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

    cut_waveform = waveform_data[0][pt_0:pt_1]
    cut_time = time_data[0][pt_0:pt_1]

    return cut_waveform, cut_time



def genFFT(waveform, sampling_rate):

    freqs, psd = signal.welch(waveform)

    func = interp1d(freqs, psd)#, kind="cubic")
    f_new = np.logspace(np.log10(freqs[1]), np.log10(freqs[-2]))

    psd = func(f_new)

    return f_new*sampling_rate, psd

def genSHM(freq, sampling_rate, time, phase=0):

    # npts = time*sampling_rate

    delta = 1/sampling_rate

    time_data = np.arange(0, time, delta)

    ampl_data = np.sin(2*np.pi*freq*time_data) + phase

    return np.array(ampl_data), np.array(time_data)

if __name__ == "__main__":

    import matplotlib.pyplot as plt

    sampling_rate = 20 #Hz
    time = 60 #seconds

    # generate a waveform as a superposition of other waves
    # a, t = genSHM(1, sampling_rate, time)
    a = np.zeros(sampling_rate*time)


    import obspy
    st = obspy.read("C:\\Users\\lmcfd\\Documents\\Fireballs\\Stubenberg\\GR_I26H4_0.mseed")
    resp = obspy.read_inventory("C:\\Users\\lmcfd\\Documents\\Fireballs\\Stubenberg\\GR_I26H4_0.xml")

    a, t = procTrace(st[0], ref_datetime=None, resp=resp, bandpass=None, backup=False)


    plt.subplot(6, 1, 1)
    plt.xlabel('Time [s]')
    plt.ylabel('Overpressure [Pa]')
    plt.subplot(6, 1, 2)

    N = 50
    FASes = []
    logs = []
    L = len(a[0])
    for i in range(N):
        plt.subplot(6, 1, 1)
        plt.plot(t[0][i*L//N:(i+1)*L//N], a[0][i*L//N:(i+1)*L//N])
        plt.subplot(6, 1, 2)
        freq, FAS = genFFT(a[0][i*L//N:(i+1)*L//N], sampling_rate)
        if i == 0:
            FAS_N = FAS
        plt.loglog(freq, FAS/FAS_N)
        FASes.append(np.sum(FAS/FAS_N))
        logs.append(FAS/FAS_N)

    best_fas = np.argmax(FASes)
    #second best fas
    FASes[best_fas] = 0
    bad_fas = np.argmax(FASes)
    print(best_fas)
    print(bad_fas)

    j = np.where(logs[best_fas] >= logs[bad_fas])

    plt.subplot(6, 1, 2)
    # plt.loglog(freq[j], logs[np.argmax(FASes)][j])

    filtered_freq = freq[j]
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Gain')
    print("Optimal Frequency Range {:.2f} - {:.2f} Hz".format(filtered_freq[0], filtered_freq[-1]))

    plt.subplot(6, 1, 3)
    a, t = procTrace(st[0], ref_datetime=None, resp=resp, bandpass=[filtered_freq[0], filtered_freq[-1]], backup=False)
    s2n = np.max(a[0])/np.median(np.abs(a[0]))
    filtered_wave = a[0]

    plt.plot(t[0], a[0], label="Optimal {:.2f}".format(s2n))
    plt.legend()

    max_p = np.nanmax(a[0])
    min_p = np.nanmin(a[0])

    print("Maximum Overpressure of Signal: {:.2f} Pa".format((max_p - min_p)/2))
    sig = a[0][j]
    sig_times = t[0][j]
    p, freq, FAS = findDominantPeriodPSD(sig, sampling_rate, normalize=False)
    plt.subplot(6, 1, 4)
    plt.semilogx(freq, FAS)
    print("Dominant Period of Signal: {:.2f} s".format(p))
    

    shortest_period = 1/sampling_rate
    longest_period = t[0][L//N] - t[0][0]

    for i in range(N):
        plt.subplot(6, 1, 5)
        max_p = np.nanmax(a[0][i*L//N:(i+1)*L//N])
        min_p = np.nanmin(a[0][i*L//N:(i+1)*L//N])

        plt.scatter(t[0][i*L//N], (max_p - min_p)/2, c='r')
        p, freq, FAS = findDominantPeriodPSD(a[0][i*L//N:(i+1)*L//N], sampling_rate, normalize=False)
        plt.subplot(6, 1, 6)
        plt.scatter(t[0][i*L//N], p, c='r')
        plt.axhline(y=shortest_period, color='k', linestyle='-')
        # plt.axhline(y=longest_period, color='k', linestyle='-')
    plt.show()
    exit()
    # for f in [0.5, 0.5, 2, 4, 8]:
    #     temp_a, t = genSHM(f, sampling_rate, time, phase=f*np.pi/2)
    #     a = a + temp_a

    # plt.subplot(1, 2, 1)
    # plt.plot(t, a, label="Raw")

    # bandpass = [1/35, 1/25]

    # a_new = bandpassFilter(bandpass, 1/sampling_rate, a)

    # plt.plot(t, a_new, label="Scipy")

    # a_obs = butter_bandpass_filter(a, bandpass[0], bandpass[1], sampling_rate)

    # plt.plot(t, a_obs, label="Scipy Cookbook")

    # plt.legend()

    # plt.subplot(1, 2, 2)


    freq, FAS = genFFT(a, sampling_rate)
    plt.semilogx()
    plt.plot(freq, FAS, label="Scipy")
    plt.show()


    # plt.axvline(x=bandpass[0]/60, c='k')
    # plt.axvline(x=bandpass[1]/60, c='k')

    # plt.legend()
    # plt.show()


    # # Sample rate and desired cutoff frequencies (in Hz).
    # fs = sampling_rate
    # lowcut = 0.3
    # highcut = 0.6

    # # Plot the frequency response for a few different orders.
    # # plt.figure(1)
    # # plt.clf()
    # # for order in [3, 6, 9]:
    # #     sos = butter_bandpass(lowcut, highcut, fs, order=order)
    # #     w, h = sosfreqz(sos, worN=2000)
    # #     plt.semilogx((fs * 0.5 / np.pi) * w, abs(h), label="order = %d" % order)

    # # plt.plot([0, 0.5 * fs], [np.sqrt(0.5), np.sqrt(0.5)], '--', label='sqrt(0.5)')
    # # plt.xlabel('Frequency (Hz)')
    # # plt.ylabel('Gain')
    # # plt.grid(True)
    # # plt.legend(loc='best')

    # # Filter a noisy signal.
    # plt.figure(1)
    # plt.clf()
    # plt.plot(t, a, label='Noisy signal')

    # # y = butter_bandpass_filter(a, lowcut, highcut, fs, order=6)
    # # plt.plot(t, y, label='Filtered signal')



    # plt.xlabel('time (seconds)')
    # plt.grid(True)
    # plt.axis('tight')
    # plt.legend(loc='upper left')

    # plt.show() 

    # fft = np.fft.rfft(a, norm="ortho")
    
    # def abs2(x):
    #     return x.real**2 + x.imag**2

    # selfconvol=np.fft.irfft(abs2(fft), norm="ortho")
    # selfconvol = selfconvol/selfconvol[0]
    # plt.plot(selfconvol)
    # plt.show()

    # # let's get a max, assuming a least 4 periods...
    # multipleofperiod=np.argmax(selfconvol[1:len(a)//4])
    # Ltrunk=a[0:(len(a)//multipleofperiod)*multipleofperiod]

    # fft = np.fft.rfft(Ltrunk, norm="ortho")
    # selfconvol=np.fft.irfft(abs2(fft), norm="ortho")
    # selfconvol=selfconvol/selfconvol[0]

    # plt.figure()
    # plt.plot(selfconvol)
    # plt.savefig('second.jpg')
    # plt.show()


    # #get ranges for first min, second max
    # fmax=np.max(selfconvol[1:len(Ltrunk)//4])
    # fmin=np.min(selfconvol[1:len(Ltrunk)//4])
    # xstartmin=1
    # while selfconvol[xstartmin]>fmin+0.2*(fmax-fmin) and xstartmin< len(Ltrunk)//4:
    #     xstartmin=xstartmin+1

    # xstartmax=xstartmin
    # while selfconvol[xstartmax]<fmin+0.7*(fmax-fmin) and xstartmax< len(Ltrunk)//4:
    #     xstartmax=xstartmax+1

    # xstartmin=xstartmax
    # while selfconvol[xstartmin]>fmin+0.2*(fmax-fmin) and xstartmin< len(Ltrunk)//4:
    #     xstartmin=xstartmin+1

    # period=np.argmax(selfconvol[xstartmax:xstartmin])+xstartmax

    # print("The period is ",period/fs)
    # # freq, FAS = genFFT(a, t)
    # freq, FAS = genFFT(a, sampling_rate)
    # plt.semilogx(freq, FAS, label="Raw")

    # freq, FAS = genFFT(y, sampling_rate)
    # plt.semilogx(freq, FAS, label="Scipy Cookbook")


    # plt.axvline(x=lowcut, c='k')
    # plt.axvline(x=highcut, c='k')

    # plt.legend()
    # plt.show()