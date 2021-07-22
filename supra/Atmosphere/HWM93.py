import hwm93

def getHWM(dt, lat, lon, elev):
	""" Gets the zonal and meridional winds at a specifc time and space

	dt - datetime obj
	lat, lon, elev - floats of position given in deg, deg, km

	"""

	winds = hwm93.run(dt, altkm=elev, glat=lat, glon=lon, f107a=150, f107=150, ap=4)

	u, v = winds.zonal.values[0], winds.meridional.values[0]

	return u, v


if __name__ == "__main__":
	
	from datetime import datetime
	u, v = getHWM(datetime(2017, 11, 12, 8), 65, -148, 150)
	print(u, v)