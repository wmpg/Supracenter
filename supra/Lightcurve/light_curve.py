
import numpy as np

class LightCurve:

    def __init__(self, station, t, h, M):

        self.station = station
        self.t = t
        self.h = h
        self.M = M
        self.I = 10**(-np.array(self.M)/2.5)

def processLightCurve(light_curve):

    t = []
    h = []
    M = []

    current_station = ""

    light_curve_list = []




    for ii, line in enumerate(light_curve):

        try:
            float(line[0])
            t.append(line[0])
            h.append(line[1])
            M.append(line[2])
        except ValueError:
            if len(t) > 0:

                L = LightCurve(current_station, t, h, M)
                light_curve_list.append(L)

            current_station = line[0]

            t = []
            h = []
            M = []


        if ii == len(light_curve) - 1:
            if len(t) > 0:

                L = LightCurve(current_station, t, h, M)
                light_curve_list.append(L)

            current_station = line[0]

            t = []
            h = []
            M = []



    return light_curve_list

def readLightCurve(csv):

    light_curve = []
    with open(csv, 'r+') as f:
    
        for line in f:
            line = line.split(',')
            temp_line = []
    
            for item in line:
                if line[0][0] == "#":

                    station = item[11:]
                    
                    if station != "":
                        temp_line.append(station)

                    break
                else:
                    temp_line.append(float(item.strip()))

    
            light_curve.append(temp_line)

    light_curve = [x for x in light_curve if x != []]

    return light_curve

if __name__ == "__main__":


    light_curve = readLightCurve(file_name)

    light_curve_list = processLightCurve(light_curve)

