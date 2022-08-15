import numpy as np


def pressureConv():
	ERA5_PRESSURE_LEVELS = np.array([1,2,3,
                5,7,10,
                20,30,50,
                70,100,125,
                150,175,200,
                225,250,300,
                350,400,450,
                500,550,600,
                650,700,750,
                775,800,825,
                850,875,900,
                925,950,975,
                1000])

	#Convert from hPa to Pa
	ERA5_PRESSURE_LEVELS *= 100

	return ERA5_PRESSURE_LEVELS

def estPressure(z):
    """ z in meters
    """

    p = 10*101.325*np.exp(-0.00012*z)*100
    # in Pa
    return p


def millibar2Pascal(mB):

    Pa = 100*mB

    return Pa


def pascal2Millibar(Pa):

    mB = Pa/100

    return mB
