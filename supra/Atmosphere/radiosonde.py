import datetime
import urllib.request

import numpy as np

from supra.Utils.Classes import Position

def parseStationList():

	file_name = 'supra/Atmosphere/igra2-station-list.txt'

	radio_stations = []

	now = datetime.datetime.now()

	with open(file_name, 'r') as f:
		for line in f:

			line = line.split()

			a = Station(line[0], line[0], Position(float(line[1]), float(line[2]), float(line[3])), 'RSNDE', line[4])

			if float(line[-2]) >= now.year:
				radio_stations.append(a)

	return radio_stations


def findNearestStation(lat, lon, stations):
	
	dist = []

	for stn in stations:
		dist.append(stn.position.ground_latlon_distance(Position(lat, lon, 0)))	

	return stations[np.nanargmin(dist)]

def buildURL(stn):
	url = 'ftp://ftp.ncdc.noaa.gov/pub/data/igra/data/data-y2d/{:}-data-beg2018.txt.zip'.format(stn.code)

	return url

def downloadRadio(lat, lon, download_loc, debug=False, year=None):

	r = parseStationList()
	closest_stat = findNearestStation(lat, lon, r)

	if debug:
		print("Closest Station: Lat: {:.4f} Lon: {:.4f}".format(closest_stat.position.lat, closest_stat.position.lon))

	url = buildURL(closest_stat)

	urllib.request.urlretrieve(url, download_loc)

