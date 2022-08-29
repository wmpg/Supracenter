
import numpy as np

class LightCurve:

    def __init__(self, station, t, h, M):

        self.station = station
        self.t = t
        self.h = h
        self.M = M
        self.I = 1500*10**(-np.array(self.M)/2.5)

    def genJoules(self):

        self.J = []

        for ii in range(len(self.t) - 1):

            dt = self.t[ii + 1] - self.t[ii]

            J = self.I[ii + 1]*dt

            self.J.append(J)

    def getDataList(self):

        return np.array(self.M), np.array(self.h)


    def estimateHeight(self, traj):

        h = []

        for t in self.t:
            h.append(traj.approxHeight(t)/1000)

        h = np.array(h)

        self.h = h

    def removeNAN(self):

        t = []
        h = []
        I = []
        M = []

        for ii in range(len(self.h)):
            if not np.isnan(self.M[ii]):
                t.append(self.t[ii])
                h.append(self.h[ii])
                I.append(self.I[ii])
                M.append(self.M[ii])

        self.t = np.array(t)
        self.h = np.array(h)
        self.I = np.array(I)
        self.M = np.array(M)

    def interpCurve(self, dh=100):

        self.removeNAN()

        dH = self.h[0] - self.h[-1]

        N = dh

        x = np.linspace(self.h[0], self.h[-1], N)

        x = x[::-1]

        xp = np.array(self.h[::-1])
        fp = np.array(self.M[::-1])


        f = np.interp(x, xp, fp)

        return x[::-1], f[::-1]


def processLightCurve(light_curve):

    t = []
    h = []
    M = []

    current_station = ""

    light_curve_list = []


    if light_curve is None:
        return None

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
                    try:
                        temp_line.append(float(item.strip()))
                    except:
                        return None

    
            light_curve.append(temp_line)

    light_curve = [x for x in light_curve if x != []]

    return light_curve

def readCNEOSlc(file_name):

    data_list = []

    t_list = [] 
    I_list = []

    latitude = None
    longitude = None


    # EXTRACT DATA
    with open(file_name, "r+") as f:
        lines = f.readlines()
        for line in lines:
            if line[0] == "(":
                line = line[1:-3]
                line = line.split(",")
                data_list.append(line)

           


    # MAKE POINTS
    for pair in data_list:
        t = float(pair[0])
        I = float(pair[1])

        t_list.append(t)
        I_list.append(I)


    t_list = np.array(t_list)
    I_list = np.array(I_list)

    # USING BROWN ET AL 1996 
    M_bol = 6 - 2.5*np.log10(I_list)

    L = LightCurve("CNEOS", t_list, None, M_bol)

    L.I = I_list
    L.M = M_bol

    return [L]

if __name__ == "__main__":


    light_curve = readLightCurve(file_name)

    light_curve_list = processLightCurve(light_curve)

