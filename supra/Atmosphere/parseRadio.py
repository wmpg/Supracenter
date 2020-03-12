import numpy as np

from supra.Utils.Classes import Constants

class Sounde():
	def __init__(self):
		pass

def openTXT(file_name, time):

	month_sounds = []
	with open(file_name, 'r') as f:
		for ii, line in enumerate(f):
			if '#' in line:
				line = line.split()
				if float(line[1]) == time[0] and float(line[2]) == time[1] and \
							np.abs(float(line[3]) - time[2]) <= 1:
					
					a = Sounde()
					a.id = line[0]
					a.year = float(line[1])
					a.month = float(line[2])
					a.day = float(line[3])
					a.hour = float(line[4])
					a.ms = float(line[5])
					a.levels = float(line[6])

					month_sounds.append(a)

	return month_sounds

def findClosestTime(soundes, time):

	ref = Sounde()
	ref.day = time[2]
	ref.hour = time[3]

	d = []
	
	for s in soundes:
		day_offset = (ref.day - s.day)*24
		hour_offset = ref.hour - s.hour

		d.append(np.abs(day_offset + hour_offset))

	return soundes[np.nanargmin(d)]

def grabClosestTime(file_name, key):
	sounding_data = []

	current_data = False

	with open(file_name, 'r') as f:
		for ii, line in enumerate(f):
			if '#' in line:
				line = line.split()
				if float(line[1]) == key.year and float(line[2]) == key.month and \
					float(line[3]) == key.day and float(line[4]) == key.hour:
					current_data = True
				else:
					current_data = False

			elif current_data:
				sounding_data.append(line)

	return sounding_data

def filterData(data, levels):

	consts = Constants()

	filt_data = np.empty((int(levels), 4))

	for ii, line in enumerate(data):
		
		filt_data[ii, 0] = float(line[16:21])
		a = float(line[22:27])/10 + 273.15
		
		try:
			filt_data[ii, 1] = (consts.GAMMA*consts.R/consts.M_0*a)**0.5
		except TypeError:
			filt_data[ii, 1] = np.nan

		filt_data[ii, 3] = (float(line[40:45]) + 180)%360
		filt_data[ii, 2] = float(line[46:51])/10

	filt_data[filt_data == -9999.0] = np.nan

	filt_data = filt_data[~np.isnan(filt_data).any(axis=1)]

	filt_data = np.flip(filt_data, axis=0)

	return filt_data

def parseRadio(file_name, time):

	soundes = openTXT(file_name, time)
	sounde = findClosestTime(soundes, time)
	data = grabClosestTime(file_name, sounde)
	filt_data = filterData(data, sounde.levels)

	return filt_data