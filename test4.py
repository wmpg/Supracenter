import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

stat_1 = np.array([1, 1, 1, 0])
stat_2 = np.array([2, 2, 2, 5])

c = 3

t = 1

#t - current time
#stat_i - arrival time
#t_i - explosion time

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

xs, ys, zs, ts = [], [], [], []
xs2, ys2, zs2, ts2 = [], [], [], []

for x_i in np.linspace(-10, 10):
	for y_i in np.linspace(-10, 10):
		for z_i in np.linspace(-10, 10):
			for t_i in np.linspace(-10, 10):
				if (x_i - stat_1[0])**2 + (y_i - stat_1[1])**2 + (z_i - stat_1[2])**2 == c**2*(t_i - stat_1[3] + t):
					xs.append(x_i)
					ys.append(y_i)
					zs.append(z_i)
					ts.append(t_i)
				if (x_i - stat_2[0])**2 + (y_i - stat_2[1])**2 + (z_i - stat_2[2])**2 == c**2*(t_i - stat_2[3] + t):
					xs2.append(x_i)
					ys2.append(y_i)
					zs2.append(z_i)
					ts2.append(t_i)

ax.scatter(xs, ys, zs, c='b')
#ax.scatter(xs2, ys2, zs2, c='r')
plt.show()