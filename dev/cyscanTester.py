import numpy as np
import sys
import time 
import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})

from supra.Supracenter.cyscan2 import cyscan
from supra.Supracenter.cyscan4 import cyscan as cyscan4

import warnings
warnings.filterwarnings("ignore")

N = 100
angle = 90

Z = np.array([[    0.0, 300.0,  90.0, 20.0], 
              [ 2000.0, 310.0,  90.0, 30.0], 
              [ 3000.0, 300.0,  90.0, 20.0], 
              [ 4000.0, 300.0,  90.0, 18.0], 
              [ 5000.0, 290.0,  94.0, 20.0], 
              [ 6000.0, 300.0,  95.0, 20.0], 
              [ 7000.0, 300.0,  90.0, 27.0], 
              [ 8000.0, 320.0, 270.0, 20.0], 
              [ 9000.0, 330.0, 180.0, 20.0], 
              [10000.0, 300.0,  90.0, 20.0]])


S = np.array([0.0, 0.0, 10000.0])
D = np.array([5000.0, 5000.0, 0.0])

t1 = time.time()
for i in range(N):
	test_1_res = cyscan(S, D, Z, wind=True, n_theta=angle, n_phi=angle, h_tol=330, v_tol=3000, debug=False)
t2 = time.time()
time_res = (t2-t1)/N

print("Cyscan 2: {:.2f} ms".format(time_res*1000))

t1 = time.time()
for i in range(N):
	test_3_res = cyscan4(S, D, Z, wind=True, n_theta=angle, n_phi=angle, h_tol=330, v_tol=3000)
t2 = time.time()
time_res2 = (t2-t1)/N

print("Cyscan 4: {:.2f} ms".format(time_res2*1000))

print("Difference in Results:")
print(test_1_res - test_3_res)

print("Speed Difference")
print("{:.2f}%".format((time_res-time_res2)/time_res*100))