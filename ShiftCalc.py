
from supra.Utils.Classes import Position
import numpy as np

def calculateShifts(array, source, speed_of_sound=330):

    shift_list = []

    for element in array:
        element_pos = element#.metadata.position
        source_pos  = source

        distance = source_pos.pos_distance(element_pos)

        time_shift = distance/speed_of_sound

        shift_list.append(time_shift)

    shift_list = shift_list - shift_list[0]

    return shift_list

if __name__ == "__main__":

    stats = [Position(45.8502, 26.6437, 701.0), Position(45.8539, 26.6454, 722.0), Position(45.8455, 26.6634, 650.0), Position(45.8419, 26.6415, 720.0)]



    source1 = Position(45.6786, 26.8569, 47140)

    source2 = Position(45.7042, 26.9041, 42760)

    source3 = Position(45.7280, 26.9480, 38690)

    source4 = Position(45.63604289860163, 26.778652254769103, 54438.40420266334)

    shift1 = calculateShifts(stats, source1)

    shift2 = calculateShifts(stats, source2)

    shift3 = calculateShifts(stats, source3)

    shift4 = calculateShifts(stats, source4)

    print(shift1)

    print(shift2)

    print(shift3)

    print(shift4)


[ 0.          0.37669839 -1.25779788 -1.06317511]