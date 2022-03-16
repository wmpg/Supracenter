
import numpy as np
from obspy.io.mseed.util import *



def shift():
	in_file = "F:\\Documents\\Meteor Research\\Event\\CrawfordBay\\IM_I56H2_0.mseed"
	out_file = "F:\\Documents\\Meteor Research\\Event\\CrawfordBay\\IM_I56H2_shift.mseed"

	# shift in seconds
	shifter = -1.65




	shift_time_of_file(in_file, out_file, -int(shifter*1e4))




if __name__ == "__main__":
	shift()