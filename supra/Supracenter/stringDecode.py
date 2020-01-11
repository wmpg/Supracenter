import random
import sys

import numpy as np

def minIndex(a):
    return list(a).index(min(a))

def createString(leng, N):
    # Original Author: Wayne Edwards

    """ Purpose: Using the uniform distribution random number generator the createstring function constructs a
        series of N random binary strings of 1's and 0's with length long.

    Arguments:
        leng: [int] length of desired binary strings
        N: [int] number of desired binary strings    

    Returned: 
        Bstring(N, long): [string] matrix of N binary strings with length long.     

    """ 
    #initialize all strings
    all_string = []

    #Loop for each string
    for i in range(N):

        # create random floats of 'leng' long
        a = np.random.uniform(0, 1, leng)
        
        # a => 0.5 -> 0 | a < 0.5 -> 1
        b = a < 0.5

        # Turn b from Boolean to integer
        b = b*1

        all_string.append(b)

    return all_string


def crossover(X):
    # Original Author: Wayne Edwards 

    """ Purpose: The crossover function takes two Binary strings randomly chosen and if their probability of 
        crossig over is less then Pc then the pair copies the information from a random location along one 
        string to the same location in the other. And vise versa. The new pair is then placed in Y. Otherwise 
        the pair remain the same and are placed into Y. This is done for n/2 pairs to create a new set of 
        binary strings

    Argument:
        X(n, m): [string] set of n binary strings of length m

    Returns:
        Y(n, m): [string] newly crossed set of binary strings
    """
    
    # Probability of crossover
    Pc = 0.90  

    n, m = X.shape
    Y = np.zeros_like(X)

    for i in range(int(n/2)):

        # select first candidate
        ii = int(np.floor(random.uniform(0, 1)*(n - 1) + 1))
        m1 = X[ii, 0:m]

        # select second candidate
        ii = int(np.floor(random.uniform(0, 1)*(n - 1) + 1))
        m2 = X[ii, 0:m]

        # See if the Pair are suitable for crossover
        cross = random.uniform(0, 1)
        if (cross <= Pc):

            # Choose location along string to crossover
            p = int(np.floor(random.uniform(0, 1)*(m - 1) + 1))

            # Crossover points
            m1[0:p], m2[0:p] = m2[0:p], m1[0:p]

        # Place pair (crossed over or not) into new set of strings
        Y[i, :] = m1
        Y[n - i - 1, :] = m2

    return Y


def mutatestrings(X):
    # Original Author: Wayne Edwards 

    """ Purpose: Using uniformly distributed random number generator the mutatestrings function takes the set 
        of binary strings, X, and tests each bit for mutation. If the random number is less than the probability 
        of mutation then the polarity of the bit is changed. Otherwise the bit is left unaltered. User has a 
        choice of three probability of mutation functions Linear, Exponential & constant. Linear & exponential 
        functions starting bits in strings higher mutation rates while constant rate is the same for all bit 
        and on the order of 1/m.

    Arguments:
        X(n, m): [string] a set of n binary strings with length m

    Returned:
        Mstring(n, m): [string] Mutated set of strings
    """
    # size of the binary strings
    n, m = X.shape
    Mstring = X

    Pm = 0.01
    Pmutate = np.random.uniform(0, 1, size=(n, m))

    for i in range(0, m):
    
        Pmutate[:, i] = (Pmutate[:, i] <= Pm)
        here = [l for l, x in enumerate(Pmutate[:, i]) if x]
    
        if len(here) != 0:
            h = len(here)
            for j in range(h):
                if (Mstring[here[j], i] == 1): 
                    Mstring[here[j], i] = 0
                else:
                    Mstring[here[j], i] = 1
    return Mstring


def reproduce(X, misfit):
    # Original Author: Wayne Edwards 
    """ Purpose: The reproduce function is designed to produce a new set of n Binary strings from an initially
        given population. The members of this new population are chosen randomly from a pool of reproductions
        of the initial population. The number of reproductions (and members) in the pool, is based on the size
        of their misfits via the probability of reproduction, Pr

    Arguments:  
        X(n, m): [string] Initial set of n Binary Strings of length m               
        misfit(n): [float] misfit for each Binary String                             
          
    Returns:
        Y(n, m): [string] New binary string set made up of copies of members in X
    """

    n, m = X.shape
    misfit = np.array(misfit)

    # Calculate Statistics of the misfits
    # Find the mean and the maximum misfit
    MaxM = np.max(misfit)
    AvgM = np.mean(misfit)

    # Linear Model for Probability of reproduction
    B = 1/(n*(MaxM - AvgM))
    A = B*MaxM
    Pr = A - B*misfit
    Pr = Pr*1000

    # Begin reproducing initial population (Breed) according to model
    try:
        breed = np.zeros((int(np.floor(max(Pr)))*n + 1))
    except:
        print("ERROR: NaN value encountered. Try increasing swarm_size.")
        sys.exit()

    P = 0

    for j in range(n):
        i = int(np.floor(Pr[j]))
        if ((i != 0) and (np.isnan(i) == False)):
            breed[P + 1 : P + i] = j
            P += i

    # By default the string with the smallest misfit gets at least one copy of itself
    # placed in the output population
    i = minIndex(misfit)
    Y = np.zeros_like(X)
    Y[0, 0:m] = X[i, 0:m]

    # Out of the Breeding Pond randomly pick output population
    P = len(breed)
    pick = np.random.uniform(0, 1, size=n - 1)
    pick = np.floor(pick*(P - 1) + 1)
    pick = pick.astype(int)

    for j in range(1, n - 1):
        i = int(breed[pick[j]])
        Y[j, 0:m] = X[i, 0:m]

    return Y


def stringDecode(code, pSize, np):
    # Original Author: Wayne Edwards 
    """ The function stringDecode takes a binary string and converts it to a series of integer positions

    Arguments:
        code: [string] binary string of given position
        pSize: [int] greatest possible size of integer position (in powers of 2)
        np: [int] number of parameters to be decoded    

    Returns:
        x: [list] array of np integers giving a position
    """
    
    # Alternating Binary String Decoding
    start = len(code)
    n = int(pSize)
    number = [1]*np

    for i in range(1, n + 1):

        for j in range(np):

            number[j] += (2**(n - i))*code[start - np + j]
        start -= np
    x = number
    
    return x
    
