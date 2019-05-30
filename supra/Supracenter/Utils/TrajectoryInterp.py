import numpy as np

from supra.Fireballs.Program import position
# The difference between both of these methods is less than 1e-3 of a degree
def trajInterp(start, end, div=10):

	A = position(*start)
	B = position(*end)

	A.pos_loc(B)
	B.pos_loc(B)

	v = np.array([B.x - A.x, B.y - A.y, B.z - A.z])

	dv = v/div

	P = [None]*div
	for i in range(div):

		P[i] = position(0, 0, 0)
		P[i].x = A.x + i*dv[0]
		P[i].y = A.y + i*dv[1]
		P[i].z = A.z + i*dv[2]

		P[i].pos_geo(B)

	P.append(B)
	
	for pt in P:
		print(pt)

	return P

def trajInterpAngle(start, az, ze, div=10):

	A = position(*start)
	A.pos_loc(A)

	azim = np.radians(az)
	zangle = np.radians(ze)

	v = np.array([np.sin(azim)*np.sin(zangle), np.cos(azim)*np.sin(zangle), -np.cos(zangle)])

	n = -A.elev/v[2]

	dv = n*v/div
	P = [None]*(div + 1)
	for i in range(div + 1):

		P[i] = position(0, 0, 0)
		P[i].x = A.x + i*dv[0]
		P[i].y = A.y + i*dv[1]
		P[i].z = A.z + i*dv[2]

		P[i].pos_geo(A)
	
	for pt in P:
		print(pt)

	return P


if __name__ == "__main__":

	start = [48.05977, 13.10846, 85920]
	end = [48.3314, 13.0706, 0]

	# start = [48.2310, 13.0846, 31504.00]
	# end = [48.2337, 13.0843, 30644.80]

	# start = [48.2299, 13.0848, 31847.68]
	# end = [48.2357, 13.0840, 30014.72]


	az = 354.67
	ze = 19.69

	points = trajInterp(start, end, div=100)
	#points_angle = trajInterpAngle(start, az, ze, div=300)

	# for i in range(len(points)):
	# 	print(points[i].lat - points_angle[i].lat, \
	# 		  points[i].lon - points_angle[i].lon,\
	# 		  points[i].elev - points_angle[i].elev)
