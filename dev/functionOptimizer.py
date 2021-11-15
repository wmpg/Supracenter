import time
import numpy as np

from supra.Supracenter.cyscan5 import cyscan as cyscan5   
from supra.Supracenter.cyscan6 import cyscan as cyscan6


def optFunc(old_func, new_func, N, args=[], kwargs={}):


    t1 = time.time()
    for i in range(N):
        old_results = old_func(*args, **kwargs)
    t2 = time.time()
    t_old = t2 - t1

    t1 = time.time()
    for i in range(N):
        new_results = new_func(*args, **kwargs)
    t2 = time.time()
    t_new = t2 - t1

    print("Function Run (old):  {:.4f} s".format(t_old/N))
    print("Function Run (new):  {:.4f} s".format(t_new/N))

    del_results = np.array(new_results) - np.array(old_results)

    if del_results.all() == 0:
        print("Results: Unchanged!")
    else:
        print("Old Results: {:}".format(old_results))
        print("New Results: {:}".format(new_results))
        print("Results: {:}".format(del_results))

    print("Fractional Difference: {:.2f}%".format((t_old - t_new)/t_old*100))

if __name__ == "__main__":
    
    S = np.array([0, 0, 1000])


    z_profile = np.array([[    0, 330, 0, 0],
                          [500, 330, 0, 0],
                          [1000, 330, 0, 0]])
    D = ([1000, 10000, 0])
    print("Actual Answers: []")


    optFunc(cyscan5, cyscan6, 100, args=[S, D, z_profile], kwargs={"wind": True})
