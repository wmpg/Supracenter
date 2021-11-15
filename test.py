
import numpy as np
from supra.Supracenter.anglescan2 import anglescan
from supra.Supracenter.cyscan5 import cyscan
from supra.Supracenter.cyscan2 import cyscan as cy2

sounding = np.array([[		   0.0,  330.0,  	0.0, 0.0], 
					 [		2000.0,  330.0,  	0.0, 0.0], 
					 [		3000.0,  300.0,  	0.0, 0.0], 
					 [		4000.0,  300.0,  	0.0, 0.0], 
					 [		5000.0,  300.0,  	0.0, 0.0], 
					 [		6000.0,  300.0,  	0.0, 0.0], 
					 [		7000.0,  300.0,  	0.0, 0.0], 
					 [		8000.0,  300.0, 	0.0, 0.0], 
					 [		9000.0,  300.0, 	0.0, 0.0], 
					 [      10000.0, 300.0,  	0.0, 0.0]])

S = np.array([0.0, 0.0, 10000.0])
D = np.array([6.41420602e-13, 1.04751934e+04, 0.00000000e+00])
theta = 135
phi = 0

a2 = anglescan(S, phi, theta, sounding, wind=True, debug=True, trace=False, plot=False)

print("A2: ", a2)

c5 = cyscan(S, D, sounding, wind=True, h_tol=330, v_tol=330)

print("C5: ", c5)

c2 = cy2(S, D, sounding, wind=True, h_tol=330, v_tol=330, n_theta=1080, n_phi=1080)

print("C2: ", c2)

tt = np.sqrt(2)*7000/300 + np.sqrt(2)*3000/330

print("Expected time: ", tt)
