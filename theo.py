import numpy as np
import matplotlib.pyplot as plt

from supra.Utils.AngleConv import chauvenet

def diff(data):
	return np.nanmax(data) - np.nanmin(data)

a = [x[:] for x in [[None] * 5] * 5]

a[0][0] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=11-ze=5.npy')
a[1][0] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=16-ze=5.npy')
a[2][0] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=21-ze=5.npy')
a[3][0] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=26-ze=5.npy')
a[4][0] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=31-ze=5.npy')
a[0][1] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=11-ze=45.npy')
a[1][1] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=16-ze=45.npy')
a[2][1] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=21-ze=45.npy')
a[3][1] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=26-ze=45.npy')
a[4][1] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=31-ze=45.npy')
a[0][2] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=11-ze=65.npy')
a[1][2] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=16-ze=65.npy')
a[2][2] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=21-ze=65.npy')
a[3][2] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=26-ze=65.npy')
a[4][2] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=31-ze=65.npy')
a[0][3] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=11-ze=85.npy')
a[1][3] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=16-ze=85.npy')
a[2][3] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=21-ze=85.npy')
a[3][3] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=26-ze=85.npy')
a[4][3] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=31-ze=85.npy')
a[0][4] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=11-ze=25.npy')
a[1][4] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=16-ze=25.npy')
a[2][4] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=21-ze=25.npy')
a[3][4] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=26-ze=25.npy')
a[4][4] = np.load('/home/luke/Desktop/Seismic_data/2016-03-06 Stubenberg fireball/v=31-ze=25.npy')

stns = ['BW-KW1', 'CZ-CKRC', 'BW-RMOA', 'BW-RTBE', 'GR-I26H1', 'BW-WETR', 'CZ-KHC', 'CZ-KHC', 'OE-MOA', 'OE-MOA', 'BW-MGS04']

lats = [48.121941, 48.8199, 47.761658, 47.745098, 48.846135, 49.14502, 49.1309, 49.1309, 47.8495, 47.8495, 47.979698]

lons = [12.597187, 13.3095, 12.864466, 12.824082, 13.71793, 12.87571, 13.5782, 13.5782, 14.2659, 14.2659, 47.979698]
elev = [512.0, 573.0, 825.0, 1111.0, 1098.0, 609.0, 700.0, 700.0, 572.0, 572.0, 607.0]

v = [11, 16, 21, 26, 31]
for stn in range(11):
	g = []
	h = []
	j = []
	k = []
	l = []
	for ptb in range(10):
		b = []
		c = []
		d = []
		e = []
		f = [] 
		for i in range(5):
			b.append(a[i][0][ptb, stn, 0, 0])
			c.append(a[i][1][ptb, stn, 0, 0])
			d.append(a[i][2][ptb, stn, 0, 0])
			e.append(a[i][3][ptb, stn, 0, 0])
			f.append(a[i][4][ptb, stn, 0, 0])

		# b_m = b[2]
		# c_m = c[2]
		# d_m = d[2]
		# e_m = e[2]
		# f_m = f[2]

		if ptb == 0:
			b_m = np.array(b[:])
			c_m = np.array(c[:])
			d_m = np.array(d[:])
			e_m = np.array(e[:])
			f_m = np.array(f[:])

			# plt.plot(v, b - b_m, c='b', label='Zenith = 5')
			# plt.plot(v, f - f_m, c='c', label='Zenith = 25')
			# plt.plot(v, c - c_m, c='r', label='Zenith = 45')
			# plt.plot(v, d - d_m, c='g', label='Zenith = 65')
			# plt.plot(v, e - e_m, c='m', label='Zenith = 85')
			
		else:
			
			g.append((np.array(b) - b_m)[0])
			h.append((np.array(f) - f_m)[0])
			j.append((np.array(c) - c_m)[0])
			k.append((np.array(d) - d_m)[0])
			l.append((np.array(e) - e_m)[0])
			# plt.plot(v, b - b_m, c='b', alpha=0.1)
			# plt.plot(v, f - f_m, c='c', alpha=0.1)
			# plt.plot(v, c - c_m, c='r', alpha=0.1)
			# plt.plot(v, d - d_m, c='g', alpha=0.1)		
	print(stns[stn])		# plt.plot(v, e - e_m, c='m', alpha=0.1)
	print("Ze =  5, Mean: {:}, Filtered Mean: {:}, Range: {:}, Filtered Range: {:}".format(np.nanmean(g), np.nanmean(chauvenet(g)[0]), diff(g), diff(chauvenet(g)[0])))
	print("Ze = 25, Mean: {:}, Filtered Mean: {:}, Range: {:}, Filtered Range: {:}".format(np.nanmean(h), np.nanmean(chauvenet(h)[0]), diff(h), diff(chauvenet(h)[0])))
	print("Ze = 45, Mean: {:}, Filtered Mean: {:}, Range: {:}, Filtered Range: {:}".format(np.nanmean(j), np.nanmean(chauvenet(j)[0]), diff(j), diff(chauvenet(j)[0])))
	print("Ze = 65, Mean: {:}, Filtered Mean: {:}, Range: {:}, Filtered Range: {:}".format(np.nanmean(k), np.nanmean(chauvenet(k)[0]), diff(k), diff(chauvenet(k)[0])))
	print("Ze = 85, Mean: {:}, Filtered Mean: {:}, Range: {:}, Filtered Range: {:}".format(np.nanmean(l), np.nanmean(chauvenet(l)[0]), diff(l), diff(chauvenet(l)[0])))
	# x = [11, 31]
	# y = [0, 0]
	# plt.plot(x, y, 'k--')
	# plt.xlabel('Velocity [km/s]')
	# plt.ylabel('Relative Arrival Time [s]')
	# degree_sign= u'\N{DEGREE SIGN}'
	# plt.title('{:} \n Lat: {:8.4f}{:}N, Lon: {:8.4f}{:}E, Elev: {:10.2f} m'.format(stns[stn], lats[stn], degree_sign, lons[stn], degree_sign, elev[stn]))
	# plt.xticks(v) 
	# plt.legend()
	# plt.show()
