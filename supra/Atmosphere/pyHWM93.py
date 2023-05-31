MODE = 93

if MODE == 93:
	import hwm93
else:
	from pyhwm2014 import HWM14

def getHWM(dt, lat, lon, elev):
	""" Gets the zonal and meridional winds at a specifc time and space

	dt - datetime obj
	lat, lon, elev - floats of position given in deg, deg, km

	"""
	try:
		if MODE == 93:
			winds = hwm93.run(dt, altkm=elev, glat=lat, glon=lon, f107a=150, f107=150, ap=4)
			u, v = winds.zonal.values[0], winds.meridional.values[0]

		else:
			dec_ut = dt.hour + dt.minute/60 + dt.second/3600
			doy = dt.timetuple().tm_yday

			winds = HWM14(alt=elev, glat=lat, glon=lon, ut=dec_ut, \
				year=dt.year, day=doy)
		return u, v
	except:
		print("[Error] Unable to get data from {:} - using no winds".format(MODE))
		return 0, 0

if __name__ == "__main__":
	
	from datetime import datetime
	u, v = getHWM(datetime(2017, 11, 12, 8), 65, -148, 150)
	print(u, v)