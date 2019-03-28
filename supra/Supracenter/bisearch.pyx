"""Finds closest value to a given value from a list"""

import cython
cimport cython

@cython.wraparound(False) # This disables negative indexing, e.g. array[-1]
@cython.boundscheck(False) # This disables checking that the indices of an array are valid
@cython.cdivision(True) # This disables checks for divisions by zero
@cython.nonecheck(False)
cpdef int bisearch(x, float val):
    #Author: Wayne Edwards

    """ Finds the index of the value in an ordered vector x that is closest to, but not greater than, value, val, given

    Arguments:
        x: [ndarray] vector of ordered values 
        val: [float] value to be searched for

    Returns:
        i: [int] position of the nearest value, without going over, to x
    """

    cdef:
        int i = 0
        int j = len(x)
        int found = 0
        int k = 0
    
    while (found == 0):
        k = (i + j)//2
        if (val < x[k]):
            j = k
        else:
            i = k
        if (j <= i + 1):
            found = 1

    return i