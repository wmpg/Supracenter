
import csv
import numpy as np
import matplotlib.pyplot as plt
import datetime

# Input files
GLM_1 = "F:\\Documents\\Meteor_Research\\Event\\Stereo-12_07_2020\\GLM\\GLM-16.csv"
GLM_2 = "F:\\Documents\\Meteor_Research\\Event\\Stereo-12_07_2020\\GLM\\GLM-17.csv"

REF_DATETIME = datetime.datetime(2020, 12, 7, 9, 52, 47)

def readCSV(file_name):
    """ Reads in GLM csv file and returns four columns without headers
    """
    data_list = []

    with open(file_name, newline="") as f:
        data = csv.reader(f, delimiter=",")
        for row in data:
            if "#" not in row[0] and len(row) > 1:
                data_list.append(row)

    f.close()

    data_list = np.array(data_list)

    time = data_list[1:, 0].astype(float)
    lon = data_list[1:, 1].astype(float)
    lat = data_list[1:, 2].astype(float)
    energy = data_list[1:, 3].astype(float)

    return time, lon, lat, energy

def energyConverter(energy):
    # See Jenniskens et al 2018


    # # Source to GLM satellite distance
    # R = 35780*1000 #m
    # #R = 42170000

    # # 4 pi r^2 : r - radius of the effective apperature
    # r = 1.61e16/0.0095

    # geo_f = 4*np.pi*R**2/(np.pi*r**2)

    # blackbody_f = 1.018e3

    time_f = 1/0.002

    blackbody = 5.67e-8*6000**4

    E = np.array(energy)

    # E = np.array(energy)*geo_f*blackbody_f*time_f

    return E

def timeConverter(time):

    # time is seconds since 1970, online docs are wrong! 
    timestamp = datetime.datetime(year=1970, month=1, day=1, hour=0, minute=0, second=0, microsecond=0)

    time_list = []
    for t in time:
        time_list.append((timestamp + datetime.timedelta(seconds=t/1e3) - REF_DATETIME).total_seconds())

    return time_list

time1, lon1, lat1, energy1 = readCSV(GLM_1)
time2, lon2, lat2, energy2 = readCSV(GLM_2)

energy1 = energyConverter(energy1)
energy2 = energyConverter(energy2)

time1 = timeConverter(time1)
time2 = timeConverter(time2)


plt.plot(time1, energy1, c="m", label="GLM-16")
plt.plot(time2, energy2, c="k", label="GLM-17")
plt.xlabel("Time after {:} [s]".format(REF_DATETIME))
plt.ylabel("Energy [J]")
plt.legend()
plt.show()

# plt.scatter(time1, lon1, c="m", label="GLM-16")
# plt.scatter(time2, lon2, c="k", label="GLM-17")
# plt.xlabel("Time after {:} [s]".format(REF_DATETIME))
# plt.ylabel("Longitude [deg E]")
# plt.legend()
# plt.show()

# plt.scatter(time1, lat1, c="m", label="GLM-16")
# plt.scatter(time2, lat2, c="k", label="GLM-17")
# plt.xlabel("Time after {:} [s]".format(REF_DATETIME))
# plt.ylabel("Latitude [deg N]")
# plt.legend()
plt.show()
