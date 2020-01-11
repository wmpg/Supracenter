"""Converts ECMWF levels into heights"""

import numpy as np

def readLevels(file_name='supra/Supracenter/level_conversion_ECMWF_37.txt', header=2):
    """ Gets the conversion of heights from a .txt file, for use with convLevels().

    Arguments:
        file_name: [string] name of conversion file, .txt
        header: [int] number of headers in the .txt file to ignore

    Returns:
        data: [ndarray] contents of the level conversion file to convert with
    """

    with open(file_name) as f:

        # Skip the header
        for i in range(header):
            next(f)

        data = np.array([0, 0, 0])

        # Parse file contents
        for line in f:

            # Remove the newline char
            line = line.replace('\n', '').replace('\r', '')

            # Split the line by the delimiter
            line = line.split()

            # Strip whitespaces from individual entries in the line
            for i, entry in enumerate(line):
                line[i] = float(entry.strip())

            # Add the contents of the line to the data list
            data = np.vstack((data, line))

        # First row was all zeroes
        data = np.delete(data, 0, 0)

        return data

def convLevels(typ=1):
    """ HELPER FUCNTION: Converts levels from ECMWF data into geopotential or geometeric heights.
        see https://www.ecmwf.int/en/forecasts/documentation-and-support/137-model-levels
        The conversion is done with a .txt file with the contents of that table.

    Arguments:
        typ: [int] 0 - levels to geopotential heights, 1 - levels to geometric heights 

    Returns:
        data: [list] list of converted heights
    """

    data = readLevels()
    
    if typ == 0:
        # Geopotential Heights
        return data[:, 1]

    else:

        # Geometric Heights
        return data[:, 2]

