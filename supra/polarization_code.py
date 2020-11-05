

import obspy
import numpy as np
import scipy

from obspy.signal.filter import bandpass


def jurkevicSum(z, n, e):

    N = np.min([len(z), len(n), len(e)])
    X = np.array([z[:N], n[:N], e[:N]])

    S = np.matmul(X, np.transpose(X))

    return S


def jurkevicPol(z, n, e):

    S = jurkevicSum(z, n, e)

    w, v = scipy.linalg.eig(S)
    
    u = np.array(v)[:,0]
    
    az = np.arctan2(u[1]*np.sign(u[0]), u[2]*np.sign(u[0]))
    
    return np.degrees(az)%360



if __name__ == "__main__":

    file_name = "C:\\Users\\lmcfd\\Desktop\\SL_CRES_14.mseed"
    st = obspy.read(file_name)
    z = st.select(channel="HHZ")[0]
    e = st.select(channel="HHE")[0]
    n = st.select(channel="HHN")[0]

    df = z.stats.sampling_rate
    low = 2
    high = 8

    z = bandpass(z, low, high, df)
    e = bandpass(e, low, high, df)
    n = bandpass(n, low, high, df)
    jurkevicPol(z, n, e)