import sys
import time 

import numpy as np

import pyximport
pyximport.install(setup_args={'include_dirs':[np.get_include()]})
from supra.Supracenter.cyscan2 import cyscan
from supra.Supracenter.anglescan import anglescan
from supra.Supracenter.anglescanrev import anglescanrev
import warnings
warnings.filterwarnings("ignore")

### CONSTANTS
len_of_test = 4
htol = 330
vtol = 1000
n_angle = 360
module_list = ["cyscan2.pyx", "anglescan.py", "anglescanrev.py"]

def result(result):

    """ Prints results in a specific color
    """

    if result == "conditional":
        return "[ " + "(!) WARNING (!)" + " ]"
    else:
        if result:
            return "[ " + "(#) PASS (#)" + " ]"
        else:
            return "[ " + "<!!> FAILURE <!!>" + " ]"

def timeres(time):
    if time >= 0.01:
        return " TIME: " + "{:.4f} s".format(time)
    elif time >= 1e-5:
        return " TIME: " + "{:.4f} ms".format(time*1000)
    elif time >= 1e-8:
        return " TIME: " + "{:.4f} µs".format(time*1e6)
    elif time >= 1e-11:
        return " TIME: " + "{:.4f} ns".format(time*1e9)
    else:
        return " TIME: " + "~{:.0f} s".format(time)

################
# CYSCAN TESTS #
################
def cyscanTests():

    # TEST 1: STRAIGHT LINE TESTS

    Z = np.array([[    0.0, 300.0, 0.0, 0.0],
                  [ 5000.0, 300.0, 0.0, 0.0],
                  [10000.0, 300.0, 0.0, 0.0]])

    count = 0
    test_1 = True
    time_res = 0
    for i in range(len_of_test):
        S = np.array([0.0, 0.0, 10000.0])
        x, y = np.random.uniform(low=-10000.0, high=10000.0, size=2)
        D = np.array([x, y, 0.0])

        t1 = time.time()
        test_1_res = cyscan(S, D, Z, wind=False, n_theta=n_angle, n_phi=n_angle, h_tol=htol, v_tol=vtol, debug=False)
        t2 = time.time()
        time_res += (t2-t1)

        v = D - S
        
        if np.isnan(test_1_res.all()):
            test_1 = "conditional"
            break
        
        # Temporal Check
        if not test_1_res[0] - np.sqrt(np.dot(v, v))/300 < np.sqrt(htol**2 + vtol**2)/300:
            test_1 = False

        # Azimuth Check
        max_allowed_az_err = np.max([\
                            np.abs(np.degrees(np.arctan2(v[0], v[1]))%360 - np.degrees(np.arctan2(v[0]+htol, v[1]+htol))%360), \
                            np.abs(np.degrees(np.arctan2(v[0], v[1]))%360 - np.degrees(np.arctan2(v[0]-htol, v[1]+htol))%360), \
                            np.abs(np.degrees(np.arctan2(v[0], v[1]))%360 - np.degrees(np.arctan2(v[0]+htol, v[1]-htol))%360), \
                            np.abs(np.degrees(np.arctan2(v[0], v[1]))%360 - np.degrees(np.arctan2(v[0]-htol, v[1]-htol))%360)])
        if not np.abs(np.degrees(np.arctan2(v[0], v[1]))%360 - test_1_res[1]%360) < max_allowed_az_err:
            test_1 = False 

        # Zenith Check
        v_h = np.sqrt(v[0]**2 + v[1]**2)
        max_allowed_ze_err = np.max([\
                             np.abs(np.degrees(np.arctan2(v[2], v_h))%360 - np.degrees(np.arctan2(v[2]+vtol, v_h+htol))%360),\
                             np.abs(np.degrees(np.arctan2(v[2], v_h))%360 - np.degrees(np.arctan2(v[2]-vtol, v_h+htol))%360),\
                             np.abs(np.degrees(np.arctan2(v[2], v_h))%360 - np.degrees(np.arctan2(v[2]+vtol, v_h-htol))%360),\
                             np.abs(np.degrees(np.arctan2(v[2], v_h))%360 - np.degrees(np.arctan2(v[2]-vtol, v_h-htol))%360)])
        
        if not np.abs(np.degrees(np.arctan(-v[2]/v_h)) + 90 - test_1_res[2]) < max_allowed_ze_err:
            test_1 = False 

    print("STRAIGHT LINE CHECK: " + result(test_1) + " " + timeres(time_res/len_of_test))

    # TEST 2: Ray Check

    S = np.array([0.0, 0.0, 10000.0])
    D = np.array([5000.0, 5000.0, 0.0])

    Z = np.array([[    0.0, 300.0,  90.0, 20.0], 
                  [ 2000.0, 300.0,  90.0, 20.0], 
                  [ 3000.0, 300.0,  90.0, 20.0], 
                  [ 4000.0, 300.0,  90.0, 20.0], 
                  [ 5000.0, 300.0,  90.0, 20.0], 
                  [ 6000.0, 300.0,  90.0, 20.0], 
                  [ 7000.0, 300.0,  90.0, 20.0], 
                  [ 8000.0, 300.0, 270.0, 20.0], 
                  [ 9000.0, 300.0, 180.0, 20.0], 
                  [10000.0, 300.0,  90.0, 20.0]])

    t1 = time.time()
    test_2_res = cyscan(S, D, Z, wind=True, n_theta=n_angle, n_phi=n_angle, h_tol=htol, v_tol=vtol, debug=True)
    t2 = time.time()

    good_res = [35.4540329,   20.01325035, 160.94708252,  25.60055176]
    if (np.abs(test_2_res[0] - good_res[0]) < (np.sqrt(htol**2 + vtol**2)/330)) and \
       (np.abs(test_2_res[1] - good_res[1]) < 2) and \
       (np.abs(test_2_res[2] - good_res[2]) < 2):

        test_2 = True
    else: 
        test_2 = False

    print("RAY-TRACE CHECK: " + result(test_2) + " " + timeres(t2-t1))

    print("SANITY CHECK: " + result(True))

###################
# ANGLESCAN TESTS #
###################
def anglescanTests():
    len_of_test = 1
    htol = 330
    vtol = 1000
    n_angle = 360

    test_1 = True
    S = np.array([0.0, 0.0, 10000.0])
    phi, theta = 20.01325035, 160.94708252
    Z = np.array([[    0.0, 300.0, 20.0,  90.0], 
                  [ 2000.0, 300.0, 20.0,  90.0], 
                  [ 3000.0, 300.0, 20.0,  90.0], 
                  [ 4000.0, 300.0, 20.0,  90.0], 
                  [ 5000.0, 300.0, 20.0,  90.0], 
                  [ 6000.0, 300.0, 20.0,  90.0], 
                  [ 7000.0, 300.0, 20.0,  90.0], 
                  [ 8000.0, 300.0, 20.0, 270.0], 
                  [ 9000.0, 300.0, 20.0, 180.0], 
                  [10000.0, 300.0, 20.0,  90.0]])

    t1 = time.time()
    res = anglescan(S, phi, theta, Z, wind=True, h_tol=htol, v_tol=vtol, target=None, debug=False)
    t2 = time.time()

    count = 0

    #horizontal check
    if not np.abs(res[0] - 5000) < htol and np.abs(res[1] - 5000) < vtol:
        test_1 = False

    #vertical check
    if not res[2] == 0:
        test_1 = False

    #time check:
    if not np.abs(res[3] - 35.4540329) < (np.sqrt(htol**2 + vtol**2)/330):
        test_1 = False


    print("RAY-TRACE CHECK: " + result(test_1) + " " + timeres(t2-t1))

    S = []
    phi = []
    theta = []

    for i in range(len_of_test):
        z = np.random.uniform(low=1000.0, high=10000.0, size=1)[0]
        S.append(np.array([0.0, 0.0, z]))
        phi.append(np.random.uniform(low=0.0, high=360.0, size=1)[0])
        theta.append(np.random.uniform(low=90.1, high=179.9, size=1)[0])

    count = 0
    time_res = 0 
    test_2 = 0
    for ii in range(len(S)):

        t1 = time.time()
        D = anglescan(S[ii], phi[ii], theta[ii], Z, wind=True, h_tol=htol, v_tol=vtol, target=None, debug=False)
        
        if np.isnan(D.all()):
            count += 2
            continue

        res = cyscan(S[ii], D[:3], Z, wind=True, n_theta=n_angle, n_phi=n_angle, h_tol=htol, v_tol=vtol, debug=False)
        t2 = time.time()
        time_res += (t2-t1)

        if np.isnan(res[0]):
            test_2 = "conditional"

        #angle check
        if np.abs(phi - res[1]) < 2 and np.abs(theta - res[2]) < 2:
            count += 1

        #time check
        if np.abs(D[3] - res[0]) < (np.sqrt(htol**2 + vtol**2)/330):
            count += 1

    if test_2 != "conditional":
        test_2 = (count == 2*len(S))
    
    print("CYSCAN COMPATABILITY CHECK: " + result(test_2) + " " + timeres(time_res/len(S)))

    print("SANITY CHECK: " + result(True))


######################
# ANGLESCANREV TESTS #
######################
def anglescanrevTests():
    
    len_of_test = 1
    htol = 330
    vtol = 1000
    n_angle = 360

    Z = np.array([[     0.0, 300.0, 0.0, 0.0],
                  [  5000.0, 300.0, 0.0, 0.0],
                  [ 10000.0, 300.0, 0.0, 0.0]])

    counter = 0
    time_res = 0

    for i in range(len_of_test):
        
        x, y = np.random.uniform(low=-10000.0, high=10000.0, size=2)
        D = np.array([x, y, 0.0])

        phi   = (np.random.uniform(low=0.0, high=360.0, size=1)[0])
        theta = (np.random.uniform(low=0.0, high=90.0, size=1)[0])

        t1 = time.time()
        S = anglescanrev(D, phi, theta, Z, wind=False, trace=False)
        t2 = time.time()
        time_res += (t2-t1)
        
        dx = S[0] - D[0]
        dy = S[1] - D[1]
        dz = S[2] - D[2]
        dh = np.sqrt(dx**2 + dy**2)

        if np.abs(np.degrees(np.arctan2(dx, dy))%360 - phi) < 0.01 and \
                    np.abs(np.degrees(np.arctan2(dz, dh))%180 - theta) < 0.01:
            counter += 1

    test_1 = (counter == len_of_test)

    print("STRAIGHT LINE TEST: " + result(test_1) + " " + timeres(time_res/len_of_test))

    # test 2
    Z = np.array([[    0.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [ 1000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [ 2000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [ 3000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [ 4000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [ 5000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [ 6000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [ 7000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [ 8000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [ 9000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]], 
                  [10000.0, np.random.uniform(low=270.0, high=360.0, size=1)[0], np.random.uniform(low=0.0, high=40.0, size=1)[0], np.random.uniform(low=0.0, high=2*np.pi, size=1)[0]]])
    time_res = 0
    test_2 = True

    for i in range(len_of_test):
        
        x, y = np.random.uniform(low=-10000.0, high=10000.0, size=2)
        D = np.array([x, y, 0.0])

        phi   = (np.random.uniform(low=0.0, high=360.0, size=1)[0])
        theta = (np.random.uniform(low=0.0, high=90.0, size=1)[0])

        t1 = time.time()
        S = anglescanrev(D, phi, theta, np.copy(Z), wind=True, trace=False)

        if np.isnan(S[0].all()) or len(S) == 1:
            test_2 = 'conditional'
            break


        t2 = time.time()
        time_res += (t2-t1)

        test_2_res = cyscan(S, D, Z, wind=True, n_theta=n_angle, n_phi=n_angle, h_tol=htol, v_tol=vtol, debug=False)

        if not np.abs(test_2_res[0] - S[-1]) < 2:
            test_2 = False
            break
    
    print("CYSCAN TIME TEST: " + result(test_2) + " " + timeres(time_res/len_of_test))

    print("SANITY CHECK: " + result(True))

if __name__ == "__main__":

    print("####################################")
    print("    SUPRACENTER DEV TESTING KIT     ")
    print("####################################")
    print("")

    print("Which module would you like to test?")

    selection_idx = None

    while selection_idx is None:
        for mm, module in enumerate(module_list):
            print('{:}. : {:}'.format(mm+1, module))

        a = input()

        try:
            selection_idx = int(a)
            if selection_idx > len(module_list):
                selection_idx = None
        except ValueError:
            print("Looking for the integer value here...")
            selection_idx = None



    selection = module_list[int(a) - 1]

    print("MODULE SELECTION: " + "{:}".format(selection))
    print("####################################")
    print("RUNNING TEST SEQUENCE...")

    if selection == module_list[0]:
        cyscanTests()
    elif selection == module_list[1]:
        anglescanTests()
    elif selection == module_list[2]:
        anglescanrevTests()
if __name__ == '__main__':

    pass