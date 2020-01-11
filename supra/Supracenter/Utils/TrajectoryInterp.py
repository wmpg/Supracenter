import numpy as np

from supra.Fireballs.Program import position
# The difference between both of these methods is less than 1e-3 of a degree




if __name__ == "__main__":

	start = [48.05977, 13.10846, 85920]
	end = [48.3314, 13.0706, 0]

	# start = [48.2310, 13.0846, 31504.00]
	# end = [48.2337, 13.0843, 30644.80]

	# start = [48.2299, 13.0848, 31847.68]
	# end = [48.2357, 13.0840, 30014.72]

	# start = [48.2355, 13.0840, 30072.00]
	# end = [48.2655, 13.0798, 20620.80]
	# start = [48.2164, 13.0867, 36086.4]
	# end = [48.2219, 13.0859, 34368.0]

	az = 354.67
	ze = 19.69

	points = trajInterp(start, end, div=500)
	#points_angle = trajInterpAngle(start, az, ze, div=300)

	# for i in range(len(points)):
	# 	print(points[i].lat - points_angle[i].lat, \
	# 		  points[i].lon - points_angle[i].lon,\
	# 		  points[i].elev - points_angle[i].elev)
