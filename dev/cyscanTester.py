import numpy as np
import sys
import time 
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Supracenter.cyscan2 import cyscan
from supra.Supracenter.cyscan5 import cyscan as cyscan5

import warnings
warnings.filterwarnings("ignore")

N = 1
angle = 90

Z = np.array([[    0.0, 330.0,  0.0, 0.0], 
              [ 2000.0, 330.0,  0.0, 0.0], 
              [ 3000.0, 330.0,  0.0, 0.0], 
              [ 4000.0, 330.0,  0.0, 0.0], 
              [ 5000.0, 330.0,  0.0, 0.0], 
              [ 6000.0, 330.0,  0.0, 0.0], 
              [ 7000.0, 330.0,  0.0, 0.0], 
              [ 8000.0, 330.0,  0.0, 0.0], 
              [ 9000.0, 330.0,  0.0, 0.0], 
              [40000.0, 330.0,  0.0, 0.0]])


S = np.array([0.0, 0000.0, 40000.0])
D = np.array([500000.0, 500000.0, 0.0])

dx = D[0] - S[0]
dy = D[1] - S[1]
dz = D[2] - S[2]

dh = np.sqrt(dy**2 + dx**2)

dist = np.sqrt((D[0] - S[0])**2 + (D[1] - S[1])**2 + (D[2] - S[2])**2)
speed = 330
etime = dist/speed

print("Expected Time:     {:.2f} s".format(etime))
print("Expected Azimuth:  {:.2f} deg".format(np.degrees(np.arctan2(dx, dy))))
print("Expected Takeoff:  {:.2f} deg".format(90 + np.degrees(np.arctan2(-dz, dh))))

t1 = time.time()
for i in range(N):
	test_1_res = cyscan(S, D, Z, wind=True, n_theta=angle, n_phi=angle, h_tol=330, v_tol=3000, debug=False)
t2 = time.time()
time_res = (t2-t1)/N

print("Cyscan 2: {:.2f} ms".format(time_res*1000))

t1 = time.time()
for i in range(N):
	test_3_res = cyscan5(S, D, Z, wind=True, h_tol=330, v_tol=3000)
t2 = time.time()
time_res2 = (t2-t1)/N

print("Cyscan 5: {:.2f} ms".format(time_res2*1000))

print("Difference in Results:")
print("Time:    {:.2f} s ({:.2f} and {:.2f})".format(test_1_res[0] - test_3_res[0], test_1_res[0], test_3_res[0]))
print("AZ:      {:.2f} deg ({:.2f} and {:.2f})".format(test_1_res[1] - test_3_res[1], test_1_res[1], test_3_res[1]))
print("TF:      {:.2f} deg ({:.2f} and {:.2f})".format(test_1_res[2] - test_3_res[2], test_1_res[2], test_3_res[2]))
print("Error:   {:.2f} ({:.2f} and {:.2f})".format(test_1_res[3] - test_3_res[3], test_1_res[3], test_3_res[3]))



print("Speed Difference")
print("{:.2f}%".format((time_res-time_res2)/time_res*100))
if __name__ == '__main__':

    pass