
import numpy as np
from obspy.io.mseed.util import *



def shift():
	in_file = "C:\\Users\\lmcfd\\Documents\\Fireballs\\Romania\\RO_IPH2_11.mseed"
	out_file = "C:\\Users\\lmcfd\\Documents\\Fireballs\\Romania\\RO_IPH2_Shift_last.mseed"

	# shift in seconds
	shifter = -1.781270





	shift_time_of_file(in_file, out_file, -int(shifter*1e4))




if __name__ == "__main__":
	shift()