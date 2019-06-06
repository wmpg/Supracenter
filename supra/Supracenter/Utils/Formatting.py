import os
import sys

def printWidth(char='#'):

    ts = os.get_terminal_size().columns
    ts = int(ts/len(char))
    print(char*ts)

def printSym(l, char='#'):
    print(char*l)

def loadingBar(message, step, max_steps):

    load_len = os.get_terminal_size().columns - len(message) - 12
    percent = step/max_steps
    load = round(percent*load_len) + 1

    sys.stdout.write("\r{:} : {:}{:} {:7.2%}{:}".format(message, '#'*(load), ' '*(load_len -load), percent, ' '))
    sys.stdout.flush()
    if percent*100 == load_len*100:
        print('')
