import os, sys
import urllib.request as urllibrary
import urllib
import datetime
import obspy
import pickle
import requests
from socket import timeout

from wmpl.Utils.OSTools import mkdirP

from supra.Stations.StationObj import Station, Metadata
from supra.Utils.Classes import Position

        # 'http://eida.gfz-potsdam.de/fdsnws/',
def urlList():
    with open(os.path.join('supra', 'Misc', 'BAMStationprefs.bam'), 'rb') as f:
        sources = pickle.load(f)

    urls = ['http://service.iris.edu/fdsnws/', 
        'http://service.ncedc.org/fdsnws/',
        'http://service.scedc.caltech.edu/fdsnws/',
        'http://rtserve.beg.utexas.edu/fdsnws/',
        'http://eida.bgr.de/fdsnws/',
        'http://eida-service.koeri.boun.edu.tr/fdsnws/',
        'http://eida.ethz.ch/fdsnws/',
        'http://geofon.gfz-potsdam.de/fdsnws/',
        'http://ws.icgc.cat/fdsnws/',
        'http://eida.ipgp.fr/fdsnws/',
        'http://webservices.ingv.it/fdsnws/',
        'http://erde.geophysik.uni-muenchen.de/fdsnws/',
        'http://eida-sc3.infp.ro/fdsnws/',
        'http://eida.gein.noa.gr/fdsnws/',
        'http://www.orfeus-eu.org/fdsnws/',
        'http://ws.resif.fr/fdsnws/',
        'http://seisrequest.iag.usp.br/fdsnws/',
        'http://fdsnws.raspberryshakedata.com/fdsnws/',
        'http://auspass.edu.au:8080/fdsnws/']

    refined_urls = []

    for i in range(len(sources)):
        if sources[i] == 1:
            refined_urls.append(urls[i])   


    return refined_urls

def filterStations(data_list):

    return data_list

def loadIntoStations(data_list):

    stn_list = []

    for station in data_list:
        meta = Metadata(station[0], station[1], \
                Position(float(station[2]), float(station[3]), float(station[4])), station[5], source=station[9])

        st = obspy.read(station[8])

        resp = obspy.read_inventory(station[10])

        stn = Station(meta, st, response=resp)

        stn_list.append(stn)

    return stn_list

def getAllStations(lat_centre, lon_centre, deg_radius, fireball_datetime, dir_path, obj=None):
    
    mkdirP(dir_path)

    urls = urlList()

    station_list = []

    sd = fireball_datetime
    start_date = "{:04d}-{:02d}-{:02d}".format(sd.year, sd.month, sd.day)
    ed = sd + datetime.timedelta(days=1)
    end_date = "{:04d}-{:02d}-{:02d}".format(ed.year, ed.month, ed.day)

    sd = fireball_datetime - datetime.timedelta(minutes=5)
    start_time = "{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}".format(sd.year, sd.month, sd.day, \
        sd.hour, sd.minute, sd.second, sd.microsecond//1000)

    ed = fireball_datetime + datetime.timedelta(minutes=55)
    end_time = "{:04d}-{:02d}-{:02d}T{:02d}:{:02d}:{:02d}.{:03d}".format(ed.year, ed.month, ed.day, \
        ed.hour, ed.minute, ed.second, ed.microsecond//1000)

    station_list = []

    print('STATIONS')

    for ii, u in enumerate(urls):
        print('Downloading from: {:}'.format(u))

        query ="{:}station/1/query?network=*&latitude={:.3f}&longitude={:.3f}" \
            "&maxradius={:.3f}&start={:s}&end={:s}&channel=*&format=text" \
            "&includerestricted=false&nodata=404".format(u, lat_centre, lon_centre, deg_radius, \
            start_date, end_date)

        try:
            txt = urllibrary.urlopen(query).read().decode('utf-8')
        except urllib.error.HTTPError as e:
            if e.code == 404:
                pass
            else:
                print(e)
            txt = ''

        if txt is not None:
            for entry in txt.split('\n')[1:]:

                entry = entry.split('|')

                # Skip empty rows
                if len(entry) != 8:
                    continue

                # Unpack the line
                network_new, station_code, lat, lon, elev, station_name, start_work, end_work = entry

                if entry not in station_list:
                    
                    # Construct a file name
                    mseed_file = network_new + '_' + station_code + '_' + str(ii) + '.mseed'
                    mseed_file_path = os.path.join(dir_path, mseed_file)

                    resp_file = network_new + '_' + station_code + '_' + str(ii) + '.xml'
                    resp_file_path = os.path.join(dir_path, resp_file)

                    stn_query = ("{:}dataselect/1/query?network={:s}&station={:s}" \
                        "&channel=*&start={:s}&end={:s}").format(u, network_new, station_code, start_time, \
                        end_time)


                    resp = ("{:}station/1/query?net={:s}&sta={:s}&starttime={:s}&"\
                                "format=xml&nodata=404&level=response".format(u, network_new, station_code, start_time))

                    try:
                        stn_query_txt = urllibrary.urlopen(stn_query)
                    except urllib.error.HTTPError as e:
                        if e.code == 404:
                            pass
                        else:
                            print(e)

                    try:
                        xml = requests.get(resp)
                        with open(resp_file_path, 'wb') as f:
                           f.write(xml.content)
                    except:
                        pass


                    print('{:2}-{:4}'.format(network_new, station_code))

                    with open(mseed_file_path,'wb') as f:
                        
                        try:
                            if stn_query_txt:
                                f.write(stn_query_txt.read())
                        except urllibrary.URLError:
                            print('Connection error! Could not download the waveform!')

                    if os.stat(mseed_file_path).st_size == 0:
                        os.remove(mseed_file_path)
                        continue

                    station_list.append([*entry, mseed_file_path, u, resp_file_path])


    station_list = filterStations(station_list)
    stn_list = loadIntoStations(station_list)

    return stn_list