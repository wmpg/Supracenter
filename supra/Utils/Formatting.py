import os
import sys
import random

import numpy as np

try:
    from termcolor import colored
    termclr = True
except ModuleNotFoundError:
    termclr = False

def splashMessage():
    print('#########################################')
    print('#     Western Meteor Python Library     #')
    print('#       Bolide Acoustic Modelling       #')
    print('#            Luke McFadden,             #')
    print('#              Denis Vida,              #') 
    print('#              Peter Brown              #')
    print('#              2018 - 2020              #')
    print('#########################################')

def printTrue(result):
    if result:
        if termclr:
            return "< " + colored('TRUE', 'green') + " >"
        else:
            return "< TRUE >"
    else:
        if termclr:
            return "< " + colored('FALSE', 'red') + " >"
        else:
            return "< FALSE >"

def termchkr(t, color='white', rm_brace=False):

    if rm_brace:
        if termclr:
            return colored(t, color)
        else:
            return t
    else:
        if termclr:
            return "[" + colored(t.upper(), color) + "] "
        else:
            return "[" + t.upper() + "] "


def printMessage(t):
    
    if t.lower() == 'debug':
        return termchkr(t, color='cyan')

    elif t.lower() == 'status':
        return termchkr(t, color='blue')

    elif t.lower() == 'warning':
        return termchkr(t, color='yellow')

    elif t.lower() == 'error':
        return termchkr(t, color='red')

    elif t.lower() == 'fragmentation':
        return termchkr(t, color='green')

    elif t.lower() == 'ballistic':
        return termchkr(t, color='blue')

    else:
        return termchkr(t, color='white')

def printPercent(num, pas):

    if pas > 4:       
        return termchkr("{:.2f}%".format(num), color='green')
    elif pas == 4:
        return termchkr("{:.2f}%".format(num), color='yellow')
    else:
        return termchkr("{:.2f}%".format(num), color='red')


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

    if percent == 1.00:
        sys.stdout.write("\r{:} : {:}{:} {:7.2%}\n".format(message, '*'*(load - 1), pattern[a], percent))
    else:
        sys.stdout.write("\r{:} : {:}{:} {:7.2%} ".format(message, '*'*(load - 1), pattern[a], percent))
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

def stationFormat(stn, setup, ref_pos, more=True):

    lines = "#"*20
    cyan_tag = termchkr("#", color="cyan", rm_brace=True)
    
    try:
        b_time = stn.times.ballistic[0][0][0]
    except:
        b_time = None


    try:
        b_prts = stn.times.ballistic[0][1][0]
    except:
        b_prts = None


    print(termchkr(lines, color="cyan", rm_brace=True), termchkr(lines, color="cyan", rm_brace=True))
    print("{:} Station: {:2}-{:5}".format(cyan_tag, stn.metadata.network, stn.metadata.code))
    print(termchkr(lines, color="cyan", rm_brace=True))
    print("{:} {:}".format(cyan_tag, stn.metadata.name))
    print("{:} Latitude  {:.4f} °N".format(cyan_tag, stn.metadata.position.lat))
    print("{:} Longitude {:.4f} °E".format(cyan_tag, stn.metadata.position.lon))
    print("{:} Elevation {:.2f}  m".format(cyan_tag, stn.metadata.position.elev))
    print(termchkr(lines, color="cyan", rm_brace=True))
    print("{:} Filtering".format(cyan_tag))
    print(termchkr(lines, color="cyan", rm_brace=True))
    print("{:} Response Attached:        {:}".format(cyan_tag, printTrue(stn.hasResponse())))
    print("{:} Seismic Available:        {:}".format(cyan_tag, printTrue(stn.hasSeismic())))
    print("{:} Infrasound Available:     {:}".format(cyan_tag, printTrue(stn.hasInfrasound())))
    print(termchkr(lines, color="cyan", rm_brace=True))
    print("{:} Sources".format(cyan_tag))
    print(termchkr(lines, color="cyan", rm_brace=True))
    print("{:} Ground Distance From Reference: {:.2f} km".format(cyan_tag, stn.metadata.position.ground_distance(ref_pos)/1000))
    # if b_prts is not None:
    #     print("{:} Ballistic Arrival {:.2f} ({:+.2f} / {:+.2f}) s".format(cyan_tag, b_time, np.nanmax(b_prts) - b_time), np.nanmax(b_prts) - b_time)
    # else:
    if b_time is None:
        print("{:} No Ballistic Arrival".format(cyan_tag))
    else:
        print("{:} Ballistic Arrival {:.2f} s".format(cyan_tag, b_time))
    

    print(termchkr(lines, color="cyan", rm_brace=True))
    print("{:} Annotations".format(cyan_tag))
    print(termchkr(lines, color="cyan", rm_brace=True))
    if len(stn.annotation.annotation_list) == 0:
        print("{:} No Annotations".format(cyan_tag))
    for an in stn.annotation.annotation_list:
        print("{:} {:} at {:.2f} s, from source: {:}".format(cyan_tag, an.title, an.time, an.source))
    print(termchkr(lines, color="cyan", rm_brace=True))

def validate(v, name):



    if isinstance(v, float):
        if v is None:
            raise TypeError("Invalid input argument: {:} (None)".format(name))
            return None

        if np.isnan(v):
            raise TypeError("Invalid input argument: {:} (np.nan)".format(name))
            return None

    return v