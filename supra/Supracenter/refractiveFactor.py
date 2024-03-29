import numpy as np
from supra.Supracenter.anglescan2 import anglescan
from supra.Supracenter.cyscan5 import cyscan
from supra.Utils.Classes import *

import matplotlib.pyplot as plt

def addAngleComp(ang1, ang2, deg=False):
    ''' Adds perpendicular angles together to a combined angle
    '''

    if np.isnan(ang1) or np.isnan(ang2):
        return np.nan

    if deg:
        return np.degrees(np.arccos(np.cos(np.radians(ang1))*np.cos(np.radians(ang2))))
    else:
        return np.arccos(np.cos(ang1)*np.cos(ang2))

def refractiveFactor(S, D, zProfile, D_ANGLE=2):
    """ Compare the difference between launch angles of rays in an ideal and a realistic atmosphere going
    to the same points on the ground (D_ANGLE away from each other in the realistic atmosphere)
    """

    # put into local coordinates
    D.pos_loc(S)
    S.pos_loc(S)

    # Find launch angles in the ideal atmosphere (straight lines)
    dx, dy, dz = D.xyz - S.xyz 
    tf_ideal_n = np.degrees(np.arctan2(dz, np.sqrt((dy)**2 + (dx)**2))) + 90
    az_ideal_n = np.degrees(np.arctan2(dx, dy))%360


    # Find actual angles
    _, az_n, tf_n, _ = cyscan(S.xyz, D.xyz, zProfile, wind=True,\
            h_tol=330, v_tol=330)


    d_angle_ideal = [np.nan]

    # Number of angles to try
    RANGE = 4

    # Go through <RANGE> number of angles <D_ANGLE> away from the actual angles
    for i in range(RANGE):
        d_az = np.degrees(np.arctan(np.cos(i*2*np.pi/RANGE)*np.tan(np.radians(D_ANGLE))))
        d_tf = np.degrees(np.arctan(np.sin(i*2*np.pi/RANGE)*np.tan(np.radians(D_ANGLE))))

        D_n = anglescan(S.xyz, az_n + d_az, tf_n + d_tf, zProfile, wind=True)

        dx, dy, dz = D_n[0:3] - S.xyz

        # Find ideal atmosphere angles to get to new launch points
        tf = np.degrees(np.arctan2(dz, np.sqrt((dy)**2 + (dx)**2))) + 90
        az = (np.degrees(np.arctan2(dx, dy)))%360

        tf = (tf - tf_ideal_n)
        az = (az - az_ideal_n) 


        d_angle_ideal.append(addAngleComp(tf, az, deg=True))



    if np.isnan(d_angle_ideal).all():
        return np.nan


    d_angle_ideal = np.nanmean(d_angle_ideal)

    rf = np.sqrt(D_ANGLE/d_angle_ideal)
    
    return rf

if __name__ == "__main__":


    z_profile = np.array([[  7.06000000e+02,   3.24486958e+02,   6.36118647e+00,  -3.13443622e+00],
 [  8.24062114e+02,   3.19690862e+02,   1.57042734e+01,  -2.70366529e+00],
 [  1.64812423e+03,   3.19598265e+02,   1.63140898e+01,  -2.65343600e+00],
 [  2.47218634e+03,   3.17416190e+02,   1.61052447e+01,  -2.74364338e+00],
 [  3.29624846e+03,   3.15115425e+02,   1.59108914e+01,  -2.67985235e+00],
 [  4.12031057e+03,   3.12890966e+02,   1.48770415e+01,  -2.65940597e+00],
 [  4.94437268e+03,   3.09459895e+02,   1.26108906e+01,  -2.67571289e+00],
 [  5.76843480e+03,   3.06016466e+02,   9.09817466e+00,  -2.64769082e+00],
 [  6.59249691e+03,   3.03425325e+02,   1.02152550e+01,  -2.86508287e+00],
 [  7.41655902e+03,   3.01635624e+02,   1.56124774e+01,  -2.94619533e+00],
 [  8.24062114e+03,   3.00744302e+02,   2.18552900e+01,  -2.65513311e+00],
 [  9.06468325e+03,   3.00268735e+02,   2.06919232e+01,  -4.27408742e+00],
 [  9.88874537e+03,   3.00031743e+02,   1.88214704e+01,  -1.36819633e+00],
 [  1.07128075e+04,   3.00252516e+02,   1.93972913e+01,   3.83339583e+00],
 [  1.15368696e+04,   3.00174710e+02,   1.73513363e+01,   2.81939763e+00],
 [  1.23609317e+04,   2.99266812e+02,   1.71769913e+01,   3.00861678e+00],
 [  1.31849938e+04,   2.98389289e+02,   1.59544666e+01,   2.74511812e+00],
 [  1.40090559e+04,   2.97635750e+02,   1.57910274e+01,   2.52657135e+00],
 [  1.48331180e+04,   2.96687997e+02,   1.37572928e+01,   2.32235797e+00],
 [  1.56571802e+04,   2.95498878e+02,   1.13667841e+01,   2.04839164e+00],
 [  1.64812423e+04,   2.94208395e+02,   1.18130845e+01,   1.68327796e+00],
 [  1.73053044e+04,   2.93024814e+02,   1.42480264e+01,   1.33837977e+00],
 [  1.81293665e+04,   2.92170176e+02,   1.68676364e+01,   1.15532294e+00],
 [  1.89534286e+04,   2.91691806e+02,   1.85930806e+01,   1.17667563e+00],
 [  1.97774907e+04,   2.91333170e+02,   1.96066418e+01,   1.27273098e+00],
 [  2.06015528e+04,   2.90836766e+02,   2.02133674e+01,   1.31079682e+00],
 [  2.14256150e+04,   2.90195352e+02,   2.06885690e+01,   1.27936199e+00],
 [  2.22496771e+04,   2.89527993e+02,   2.12925499e+01,   1.22807758e+00],
 [  2.30737392e+04,   2.88954728e+02,   2.22854982e+01,   1.20706447e+00],
 [  2.38978013e+04,   2.88543635e+02,   2.38739144e+01,   1.25081116e+00],
 [  2.47218634e+04,   2.88128452e+02,   2.60221654e+01,   1.32330305e+00],
 [  2.55459255e+04,   2.87476885e+02,   2.86263878e+01,   1.36865860e+00],
 [  2.63699876e+04,   2.86456286e+02,   3.15283880e+01,   1.34605268e+00],
 [  2.71940498e+04,   2.85358836e+02,   3.43383663e+01,   1.27885404e+00],
 [  2.80181119e+04,   2.84589910e+02,   3.66048145e+01,   1.20753516e+00],
 [  2.88421740e+04,   2.84554881e+02,   3.78762236e+01,   1.17256869e+00],
 [  2.96662361e+04,   2.85659126e+02,   3.77010847e+01,   1.21442726e+00],
 [  3.04902982e+04,   2.88206622e+02,   3.57866391e+01,   1.35047805e+00],
 [  3.13143603e+04,   2.92063865e+02,   3.25250840e+01,   1.49839570e+00],
 [  3.21384224e+04,   2.96974158e+02,   2.84911777e+01,   1.55030282e+00],
 [  3.29624846e+04,   3.02356901e+02,   2.41274713e+01,   1.48072490e+00],
 [  3.37865467e+04,   3.07123406e+02,   1.96691229e+01,   1.39345015e+00],
 [  3.46106088e+04,   3.10218507e+02,   1.53486181e+01,   1.39632079e+00],
 [  3.54346709e+04,   3.11701099e+02,   1.16194159e+01,   1.49359076e+00],
 [  3.62587330e+04,   3.12461776e+02,   9.09994144e+00,   1.61218090e+00],
 [  3.70827951e+04,   3.13411082e+02,   8.41257663e+00,   1.67715711e+00],
 [  3.79068572e+04,   3.15457597e+02,   1.01768733e+01,   1.61377318e+00],
 [  3.87309193e+04,   3.19173832e+02,   1.45280557e+01,   1.37943963e+00],
 [  3.95549815e+04,   3.24418530e+02,   2.05726945e+01,   9.99864040e-01],
 [  4.03790436e+04,   3.30959210e+02,   2.72858867e+01,   5.09483119e-01],
 [  4.12031057e+04,   3.38550709e+02,   3.36582454e+01,  -5.65824567e-02],
 [  4.20271678e+04,   3.46513919e+02,   3.92114643e+01,  -6.39800750e-01],
 [  4.28512299e+04,   3.53637912e+02,   4.41180950e+01,  -1.15294849e+00],
 [  4.36752920e+04,   3.58679257e+02,   4.85904719e+01,  -1.50704870e+00],
 [  4.44993541e+04,   3.60394519e+02,   5.28409293e+01,  -1.61312441e+00],
 [  4.53234163e+04,   3.57540267e+02,   5.70818014e+01,  -1.38219864e+00],
 [  4.61474784e+04,   3.48873067e+02,   6.15254225e+01,  -7.25294411e-01],
 [  4.69715405e+04,   3.33230808e+02,   6.63696905e+01,   4.40790654e-01],
 [  4.77956026e+04,   3.18463918e+02,   7.02125991e+01,   1.55928936e+00],
 [  4.86196647e+04,   3.19491970e+02,   7.05341245e+01,   1.56760500e+00],
 [  4.94437268e+04,   3.21332068e+02,   7.05140434e+01,   1.51507895e+00],
 [  5.02677889e+04,   3.21914073e+02,   7.05769363e+01,   1.54475497e+00],
 [  5.10918511e+04,   3.22949087e+02,   7.03716495e+01,   1.53359211e+00],
 [  5.19159132e+04,   3.23637528e+02,   7.00617654e+01,   1.53510952e+00],
 [  5.27399753e+04,   3.24178122e+02,   6.95981961e+01,   1.53260540e+00],
 [  5.35640374e+04,   3.24481088e+02,   6.89850908e+01,   1.52971763e+00],
 [  5.43880995e+04,   3.24511770e+02,   6.82585399e+01,   1.52628016e+00],
 [  5.52121616e+04,   3.24236841e+02,   6.74463883e+01,   1.52211410e+00],
 [  5.60362237e+04,   3.23646137e+02,   6.65720516e+01,   1.51763512e+00],
 [  5.68602858e+04,   3.22769158e+02,   6.56612814e+01,   1.51295793e+00],
 [  5.76843480e+04,   3.21645129e+02,   6.47383532e+01,   1.50832888e+00],
 [  5.85084101e+04,   3.20312187e+02,   6.38266006e+01,   1.50398952e+00],
 [  5.93324722e+04,   3.18808332e+02,   6.29571463e+01,   1.50016531e+00],
 [  6.01565343e+04,   3.17170651e+02,   6.21238249e+01,   1.49728114e+00],
 [  6.09805964e+04,   3.15434746e+02,   6.13269043e+01,   1.49557301e+00],
 [  6.18046585e+04,   3.13634486e+02,   6.05802187e+01,   1.49489658e+00],
 [  6.26287206e+04,   3.11801776e+02,   5.98633181e+01,   1.49501921e+00],
 [  6.34527828e+04,   3.09966458e+02,   5.91518444e+01,   1.49573906e+00],
 [  6.42768449e+04,   3.08156332e+02,   5.84215612e+01,   1.49684306e+00],
 [  6.51009070e+04,   3.06397251e+02,   5.76479625e+01,   1.49810373e+00],
 [  6.59249691e+04,   3.04713285e+02,   5.68071593e+01,   1.49929495e+00],
 [  6.67490312e+04,   3.03126907e+02,   5.58740374e+01,   1.50016256e+00],
 [  6.75730933e+04,   3.01659241e+02,   5.48218028e+01,   1.50024923e+00],
 [  6.83971554e+04,   3.00330411e+02,   5.36387396e+01,   1.49877364e+00],
 [  6.92212176e+04,   2.99159441e+02,   5.23298999e+01,   1.49638487e+00],
 [  7.00452797e+04,   2.98166062e+02,   5.08952341e+01,   1.49380559e+00],
 [  7.08693418e+04,   2.97367246e+02,   4.93271564e+01,   1.49133443e+00],
 [  7.16934039e+04,   2.96783782e+02,   4.76203357e+01,   1.48948978e+00],
 [  7.25174660e+04,   2.96439674e+02,   4.57682088e+01,   1.48882171e+00],
 [  7.33415281e+04,   2.96257535e+02,   4.37683100e+01,   1.49003490e+00],
 [  7.41655902e+04,   2.96173139e+02,   4.16024384e+01,   1.49386135e+00],
 [  7.49896524e+04,   2.96163771e+02,   3.93034862e+01,   1.50169934e+00],
 [  7.58137145e+04,   2.96218851e+02,   3.69246469e+01,   1.51488428e+00],
 [  7.66377766e+04,   2.96324266e+02,   3.44518749e+01,   1.53339596e+00],
 [  7.74618387e+04,   2.96466816e+02,   3.18911767e+01,   1.55710182e+00],
 [  7.82859008e+04,   2.96633219e+02,   2.92618786e+01,   1.58602086e+00],
 [  7.91099629e+04,   2.96810097e+02,   2.65830318e+01,   1.62022235e+00],
 [  7.99340250e+04,   2.96984038e+02,   2.38692400e+01,   1.65979753e+00],
 [  7.99800000e+04,   2.96984038e+02,   2.38692400e+01,   1.65979753e+00]])
    
    D = Position(45.852, 26.6646, 706)

    
    A = Trajectory(-2.47, 27760, pos_i=Position(45.4560, 26.45, 85500), azimuth=Angle(52.20), zenith=Angle(47.00))

    points = A.trajInterp2(div=50, min_p=17000, max_p=80000, xyz=False)

    x = []
    y = []


    for pt in points:
        S = Position(pt[0], pt[1], pt[2])

        r = refractiveFactor(S, D, z_profile)
        print(pt[2], r)
        x.append(pt[2])
        y.append(r)


    plt.scatter(x, y)
    plt.show()