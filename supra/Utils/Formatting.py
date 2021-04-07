import os
import sys
import random

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

def termchkr(t, color='white'):

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

def printPercent(num):

    if num >= 67:       
        return termchkr("{:.2f}%".format(num), color='green')
    elif num >= 33:
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

def stationFormat(network, code, available_channels, ground_distance):

    a = (" # # # {:}-{:} # # # ").format(network, code)
    print(a)
    print("#"*len(a))
    for chn in available_channels:
        print(chn)

    print("Ground Distance: {:7.3f} km".format(ground_distance/1000))
    print("")
