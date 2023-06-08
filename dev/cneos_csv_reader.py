from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np

from global_land_mask import globe

FILENAME = "F:\\Desktop\\cneos_fireball_data.csv"

def collisionDetection(lat, lon, bounding_box):
    #bounding box given as (left, top, right, bottom)

    left, top, right, bottom = bounding_box

    if lon <= right and lon >= left and lat >= bottom and lat <= top:
        return True

    return False

class CNEOS():

    def __init__(self, data):

        self.time = data[0].strip('"')
        self.lat = data[1].strip('"')
        self.lon = data[2].strip('"')
        self.alt = data[3].strip('"')
        self.E = data[9].strip('"')
        self.region = None

        if len(self.lat) == 0:
            self.lat = None
        else:

            if self.lat[-1] == "N":
                self.lat = float(self.lat[:-1])
            elif self.lat[-1] == "S":
                self.lat = -float(self.lat[:-1])
            else:
                self.lat = None

            if self.lon[-1] == "E":
                self.lon = float(self.lon[:-1])
            elif self.lon[-1] == "W":
                self.lon = -float(self.lon[:-1])
            else:
                self.lon = None

        try:
            self.alt = float(self.alt)
        except:
            self.alt = None

        try:
            self.E = float(self.E)
        except:
            self.E = None

    def __str__(self):

        return "{:}: ({:.2f} N, {:.2f} E) {:.2f} kT".format(self.time, self.lat, self.lon, self.E)


hydroacoustic_stations = [[115.2, -34.3], [-132.5, 53.3], [-78.8, -33.6], [51.9, -46.4], [-61.1, 16.3], [-110.9, 18.7], [-31.2, 39.4], [72.4, -7.3], [-12.3, -37.1], [-14.4, -8.0], [166.6, 19.3]]
#can_station 53°18'00.0"N 132°30'00.0"W
pacific_E_bb = [-180, 61, -116, 10]
pacific_W_bb = [149, 61, 180, 10]
atlantic_bb = [-67, 64, -16, 10]

obj_list = []




with open(FILENAME) as f:
    lines = f.readlines()
    total_events = len(lines) - 1
    for line in lines:
        line = line.strip()

        data = line.split(",")

        obj = CNEOS(data)

        if obj.lat is None:
            continue
        else:
            obj_list.append(obj)

        if globe.is_ocean(obj.lat, obj.lon):
            obj.region = "pacific"
        # if collisionDetection(obj.lat, obj.lon, pacific_E_bb):

        #     obj.region = "pacific"

        # elif collisionDetection(obj.lat, obj.lon, pacific_W_bb):

        #     obj.region = "pacific"

        # elif collisionDetection(obj.lat, obj.lon, atlantic_bb):

        #     obj.region = "atlantic"

E_cutoff = 1 #kT

p_list = []
for obj in obj_list:
    if obj.region == "pacific" and obj.E > E_cutoff:
        p_list.append(obj)

p_list_sorted = sorted(p_list, key=lambda x: x.E, reverse=True)

a_list = []
for obj in obj_list:
    if obj.region == "atlantic" and obj.E > E_cutoff:
        a_list.append(obj)

a_list_sorted = sorted(a_list, key=lambda x: x.E, reverse=True)

print("######################")
print("Total List - {:} Events > {:} kT".format(len(p_list) + len(a_list), E_cutoff))
print("######################")


# print("Pacific List - {:} Events".format(len(p_list)))
# print("######################")
for p in p_list_sorted:
    print(p)
print("######################")


# print("Atlantic List - {:} Events".format(len(a_list)))
# print("######################")
# for a in a_list_sorted:
#     print(a)
# print("######################")



m = Basemap(projection='merc', \
            llcrnrlat=-70,\
            urcrnrlat=70, \
            llcrnrlon=0, \
            urcrnrlon=360, \
            lat_ts=1, \
            resolution='l')

m.fillcontinents(color='grey', lake_color='aqua')
m.drawcountries(color='black')
m.drawlsmask(ocean_color='aqua')

m.drawparallels(np.arange(-70, 70, 10), labels=[1,0,0,1], textcolor="k", fmt="%.1f")
meridians = m.drawmeridians(np.arange(-180, 180, 15), labels=[1,0,0,1], textcolor="k", rotation="horizontal", fmt="%.1f")

for mmmmmm in meridians:
    try:
        meridians[mmmmmm][1][0].set_rotation(45)
    except:
        pass
for p in p_list_sorted + a_list_sorted:
    x, y = m(p.lon%360, p.lat)
    if p.E > 10:
        plt.scatter(x, y, c='r', s=20*np.log10(p.E))
    elif p.E > 1:
        plt.scatter(x, y, c='orange', s=20*np.log10(p.E))
    else:
        plt.scatter(x, y, c='g', s=20*np.log10(p.E))

for stat in hydroacoustic_stations:
    x, y = m(stat[0]%360, stat[1])

    plt.scatter(x, y, c='b', marker='^')

plt.scatter(0, 0, c='r', label="> 10 kT", alpha=0)
plt.scatter(0, 0, c='orange', label="> 1 kT", alpha=0)
plt.scatter(0, 0, c='g', label="> 0.1 kT", alpha=0)
plt.scatter(0, 0, c='b', marker='^', alpha=0, label="Hydroacoustic Station")

collected_events = len(p_list) + len(a_list)
plt.title("{:} Events (out of {:} {:.2f}%) > {:} kT".format(collected_events, total_events, collected_events/total_events*100, E_cutoff))
leg = plt.legend()
for lh in leg.legendHandles: 
    lh.set_alpha(1)
plt.show()