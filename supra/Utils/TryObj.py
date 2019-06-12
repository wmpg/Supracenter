import datetime

from supra.Utils.Classes import Angle, Position, Trajectory

def strToBool(my_str):

    return (my_str.lower() == 'true')

def tryFloat(statement):
    
    try:
        return float(statement)
    except:
        return None

def tryInt(statement):

    try:
        return int(statement)
    except:
        return None

def tryBool(statement):

    statement = statement.lower()

    return strToBool(statement)

def tryStr(statement):

    return statement

def tryDateTime(statement):

    try:
        try:
            return datetime.datetime.strptime(statement, "%Y-%m-%d %H:%M:%S.%f")
        except:
            return datetime.datetime.strptime(statement, "%Y-%m-%d %H:%M:%S")
    except:
        return None

def tryEval(statement, try_float=False):
    
    try:
        statement = eval(statement)
    except:
        statement = None

    if try_float:
        statement = [tryFloat(i) for i in statement]

    return statement

def tryPosition(lat, lon, elev):

    try:
        result = Position(lat, lon, elev)
    except:
        result = None

    return result

def tryAngle(angle):

    try:
        result = Angle(angle)
    except:
        result = None

    return result

def tryTrajectory(t, v, az, ze, i, f):

    try:
        result = Trajectory(t, v, azimuth=az, zenith=ze, pos_i=i, pos_f=f)
    except:
        result = None

    return result

def trySupracenter(statement):

    supra_list = []
    
    for event in statement:
        supra_list.append(Supracenter(Position(event[0], event[1], event[2]), event[3]))

    return supra_list