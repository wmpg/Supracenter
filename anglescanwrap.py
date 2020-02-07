import numpy as np
from supra.Supracenter.anglescan import anglescan
from supra.Utils.Classes import Position

S = np.array([1069.86, -11475.30, 32400])
#az, tf
phi, theta = 32.1531715393, 110.390098572 

z_profile = np.array(
[[  1.09800000e+03,   3.33249621e+02,   6.39769967e+00,   8.15753052e-01],
 [  1.32870000e+03,   3.32245814e+02,   7.46880547e+00,   8.21841037e-01],
 [  1.45991000e+03,   3.31259057e+02,   7.75688260e+00,   8.19828570e-01],
 [  1.75063000e+03,   3.30386460e+02,   8.42901661e+00,   8.28831141e-01],
 [  2.08109000e+03,   3.29449692e+02,   8.86068035e+00,   8.06714981e-01],
 [  2.26180000e+03,   3.28453268e+02,   8.97530606e+00,   7.34630784e-01],
 [  2.65469000e+03,   3.27410620e+02,   9.00045710e+00,   6.62316170e-01],
 [  3.08925000e+03,   3.25090182e+02,   8.78484173e+00,   5.17458770e-01],
 [  3.81482000e+03,   3.22613843e+02,   9.23816624e+00,   4.71988932e-01],
 [  4.34173000e+03,   3.19867468e+02,   1.02598161e+01,   4.81650572e-01],
 [  4.89602000e+03,   3.16732426e+02,   1.13200487e+01,   4.82942404e-01],
 [  5.75930000e+03,   3.13210326e+02,   1.17077738e+01,   5.04343430e-01],
 [  6.31500000e+03,   3.13210326e+02,   1.17077738e+01,   5.04343430e-01],
 [  6.63166000e+03,   3.09223385e+02,   1.23709315e+01,   5.72768257e-01],
 [  7.21409000e+03,   3.04648244e+02,   1.19807062e+01,   5.62036780e-01],
 [  8.38036000e+03,   2.99650957e+02,   1.19626229e+01,   3.42363272e-01],
 [  9.25570000e+03,   2.96384905e+02,   1.34617191e+01,   4.15566865e-01],
 [  1.04226400e+04,   2.97597801e+02,   1.72666631e+01,   6.36915559e-01],
 [  1.12979300e+04,   2.98152179e+02,   1.85458983e+01,   8.13275332e-01],
 [  1.15320000e+04,   2.98152179e+02,   1.85458983e+01,   8.13275332e-01],
 [  1.18902400e+04,   2.98841966e+02,   1.63087678e+01,   9.13307592e-01],
 [  1.27973000e+04,   2.99356526e+02,   1.61640549e+01,   8.97905113e-01],
 [  1.37271800e+04,   2.99292602e+02,   1.68634748e+01,   9.74438712e-01],
 [  1.50035000e+04,   2.98562063e+02,   1.65094500e+01,   1.08703972e+00],
 [  1.63228300e+04,   2.96805638e+02,   1.48915123e+01,   1.39067664e+00],
 [  1.67490000e+04,   2.96727136e+02,   1.51893465e+01,   1.39105321e+00],
 [  1.87563400e+04,   2.94340057e+02,   1.16511671e+01,   1.33977629e+00],
 [  2.06949000e+04,   2.93534851e+02,   1.68547503e+01,   1.62653654e+00],
 [  2.19660000e+04,   2.93534851e+02,   1.68547503e+01,   1.62653654e+00],
 [  2.39000200e+04,   2.90542342e+02,   1.85821413e+01,   1.74761806e+00],
 [  2.66355600e+04,   2.90029354e+02,   3.13705767e+01,   1.47533555e+00],
 [  2.71830000e+04,   2.90029354e+02,   3.13705767e+01,   1.47533555e+00],
 [  3.13309600e+04,   2.96627600e+02,   4.80634166e+01,   1.57133120e+00],
 [  3.24000000e+04,   2.96627600e+02,   4.80634166e+01,   1.57133120e+00]]
)
z_profile[:, 1] = 330
D = anglescan(S, phi, theta, z_profile, wind=False)
print(D)
ref_pos = Position(48.3314, 13.0706, 0)
a = Position(0, 0, 0)
a.x = D[0]
a.y = D[1]
a.z = D[2]

a.pos_geo(ref_pos)
print(a)
