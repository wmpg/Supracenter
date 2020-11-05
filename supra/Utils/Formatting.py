import os
import sys
import random

def splashMessage():
    print('#########################################')
    print('#     Western Meteor Python Library     #')
    print('#       Bolide Acoustic Modelling       #')
    print('#            Luke McFadden,             #')
    print('#              Denis Vida,              #') 
    print('#              Peter Brown              #')
    print('#              2018 - 2020              #')
    print('#########################################')

def printWidth(char='#'):

    ts = os.get_terminal_size().columns
    ts = int(ts/len(char))
    print(char*ts)

def printSym(l, char='#'):
    print(char*l)

def loadingBar(message, step, max_steps):

    load_len = os.get_terminal_size().columns - len(message) - 12
    percent = step/max_steps
    load = round(percent*load_len) - 5
    pattern = ['¯-_0}', '-_¯0}', '_¯-0}', '-¯_0}', '¯_-0}', '_-¯0}']
    a = random.randint(0, len(pattern)-1)

    sys.stdout.write("\r{:} : {:}{:} {:7.2%} ".format(message, '*'*(load), pattern[a], percent))
    sys.stdout.flush()


def randospinny(message):
    sym = ['\\', '/', '|', '_', '-', '=', '+', '#', 'o']
    a = random.randint(0, len(sym)-1)

    sys.stdout.write('\r' + message + '   ' + sym[a])
    sys.stdout.flush()

def goodspinny(message, count, pattern=['+', 'x']):
    count = count%len(pattern)

    sys.stdout.write('\r' + message + '   ' + pattern[count])
    sys.stdout.flush()

def meteorspinny(message, count, pattern=['¯-_0}', '-_¯0}', '_¯-0}', '-¯_0}', '¯_-0}', '_-¯0}']):
    # count is the percent done)

    load_len = os.get_terminal_size().columns - len(message)
    a = random.randint(0, len(pattern)-1)
    count = count%load_len + len(pattern[0])
    sys.stdout.write(('\r' + message + ' '*count + pattern[a] + ' '*(load_len - count))[:load_len + len(message) + 1])
    sys.stdout.flush()

def checkExt(file, ext):
    
    if ext not in file:
        return file + ext
    else:
        return file