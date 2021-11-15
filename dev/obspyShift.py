
import numpy as np
from obspy.io.mseed.util import *



def shift():
	in_file = "F:\\Desktop\\Romania\\RO_IPH3_11.mseed"
	out_file = "F:\\Desktop\\Romania\\RO_IPH3_Shift_small.mseed"

	# shift in seconds
	shifter = -2.397




	shift_time_of_file(in_file, out_file, -int(shifter*1e4))




if __name__ == "__main__":
	shift()